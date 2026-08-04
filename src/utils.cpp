
#include <array>
#include <algorithm>
#include <stdexcept>

#include "utils.h"

#if defined(__x86_64__)
  #include <immintrin.h>
#endif

namespace bgen {

#if defined(__x86_64__)

// sum ploidy states with AVX2. Compiled for AVX2 in isolation, so that this is
// only ever entered after a runtime check for AVX2 support.
BGEN_TARGET_AVX2
static std::uint64_t ploidy_sum_avx2(std::uint8_t * x, std::uint32_t & size, std::uint32_t & i) {
  std::uint64_t arr[4];
  __m256i _sum = _mm256_setzero_si256();
  const __m256i zero = _mm256_setzero_si256();
  for (; i + 32 <= size; i += 32) {
    // sum each half of every 64-bit lane against zero, which accumulates the bytes
    // straight into 64-bit counters, so nothing can overflow part way through
    __m256i values = _mm256_loadu_si256((const __m256i*) &x[i]);
    _sum = _mm256_add_epi64(_sum, _mm256_sad_epu8(values, zero));
  }
  _mm256_storeu_si256((__m256i*) &arr[0], _sum);
  return arr[0] + arr[1] + arr[2] + arr[3];
}

// this code is solely here to avoid a bug on macosx x86_64 when AVX2 is not
// available. If not present, the final stage to clean up the remainder
// segfaults. It's a mystery why.
BGEN_TARGET_SSE4
static std::uint64_t ploidy_sum_sse4(std::uint8_t * x, std::uint32_t & size, std::uint32_t & i) {
  std::uint64_t arr[2];
  __m128i _sum = _mm_setzero_si128();
  const __m128i zero = _mm_setzero_si128();
  for (; i + 16 <= size; i += 16) {
    __m128i values = _mm_loadu_si128((const __m128i*) &x[i]);
    _sum = _mm_add_epi64(_sum, _mm_sad_epu8(values, zero));
  }
  _mm_storeu_si128((__m128i*) &arr[0], _sum);
  return arr[0] + arr[1];
}

// get min and max of ploidy values with AVX2
BGEN_TARGET_AVX2
static void range_avx2(std::uint8_t * x, std::uint32_t & size, size_t & i,
                       std::uint8_t & min_val, std::uint8_t & max_val) {
  std::array<std::uint8_t, 32> arr;
  __m256i values;
  __m256i _mins = _mm256_set_epi8(-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                                  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                                  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
  __m256i _maxs = _mm256_set_epi8(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                  0, 0, 0, 0);
  for (; i + 32 < size; i += 32) {
    values = _mm256_loadu_si256((const __m256i*) &x[i]);
    _mins = _mm256_min_epu8(_mins, values);
    _maxs = _mm256_max_epu8(_maxs, values);
  }
  _mm256_storeu_si256((__m256i*) &arr[0], _mins);
  for (auto v : arr) {
    min_val = std::min(min_val, v);
  }
  _mm256_storeu_si256((__m256i*) &arr[0], _maxs);
  for (auto v : arr) {
    max_val = std::max(max_val, v);
  }
}

// this code is solely here to avoid a bug on macosx x86_64 when AVX2 is not
// available. If not present, the final stage to clean up the remainder
// segfaults. It's a mystery why.
BGEN_TARGET_SSE4
static void range_sse4(std::uint8_t * x, std::uint32_t & size, size_t & i,
                       std::uint8_t & min_val, std::uint8_t & max_val) {
  std::array<std::uint8_t, 16> arr;
  __m128i values;
  __m128i _mins = _mm_set_epi8(-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                               -1, -1, -1, -1, -1);
  __m128i _maxs = _mm_set_epi8(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
  for (; i + 16 < size; i += 16) {
    values = _mm_loadu_si128((const __m128i*) &x[i]);
    _mins = _mm_min_epu8(_mins, values);
    _maxs = _mm_max_epu8(_maxs, values);
  }
  _mm_storeu_si128((__m128i*) &arr[0], _mins);
  for (auto v : arr) {
    min_val = std::min(min_val, v);
  }
  _mm_storeu_si128((__m128i*) &arr[0], _maxs);
  for (auto v : arr) {
    max_val = std::max(max_val, v);
  }
}

#endif

// Returns value of Binomial Coefficient C(n, k)
std::uint32_t n_choose_k(int n, int k) {
  std::uint32_t res = 1;

  // Since C(n, k) = C(n, n-k)
  if ( k > n - k ) {
    k = n - k;
  }

  // Calculate value of
  // [n * (n-1) *---* (n-k+1)] / [k * (k-1) *----* 1]
  for (std::uint32_t i = 0; i < (std::uint32_t) k; ++i) {
    res *= (n - i);
    res /= (i + 1);
  }

  return res;
}

/// check if the minor allele is certain (to 99.9999999999999& confidence)
///
///  Take the frequency, and number of individuals checked so far, and see if the
///  99.99..(fifteen nines) confidence interval overlaps 0.5. If not, then we can
///  be sure we've identified the minor allele, even without checking the full
///  population.
///
///  @param freq estimated minor allele frequency
///  @param n_checked number of individsuals checked so far
///  @param z standard normal deviate (eg 1.96 for 95% CI, here we use 10.0 for
///    stronger confidence, and the fact the normal approximation for confidence
///    intervals isn't perfect)
///  @return True/False for whether to halt the permuations
bool minor_certain(double freq, int n_checked, double z) {
    double delta = (z * std::sqrt((freq * (1 - freq)) / n_checked));
    
    // check if the confidence interval overlaps 0.5
    return !((freq - delta < 0.5) & (freq + delta > 0.5));
}

// sum array of 8-bit uints with vectorised operations
//
// Summing the ploidy array via numpy is more expensive than it should be, this
// is just a quick vectorized sum to speed that up.
//
// The bytes are accumulated straight into 64-bit counters, so no intermediate can
// overflow. A 32-bit accumulator could not manage that: the ploidy is at most 255
// per person, so the lanes and the sum of the lanes both wrap partway through a
// large cohort, and the size at which that happens depends on the ploidy values.
//
/// @param x array of floats
/// @param size size of array
/// @returns sum of array
std::uint64_t fast_ploidy_sum(std::uint8_t * x, std::uint32_t & size) {
  std::uint32_t i = 0;
  std::uint64_t total = 0;

#if defined(__x86_64__)
  if (__builtin_cpu_supports("avx2")) {
    total += ploidy_sum_avx2(x, size, i);
  } else if (__builtin_cpu_supports("sse4.1")) {
    total += ploidy_sum_sse4(x, size, i);
  }
#endif

  // include the remainder not used during vectorised operations
  for ( ; i < size; i++) {
    total += x[i];
  }
  return total;
}

// get min and max of the ploidy values in one pass
Range fast_range(std::uint8_t * x, std::uint32_t & size) {
  std::uint8_t min_val = 255;
  std::uint8_t max_val = 0;
  size_t i = 0;

#if defined(__x86_64__)
  if (__builtin_cpu_supports("avx2")) {
    range_avx2(x, size, i, min_val, max_val);
  } else if (__builtin_cpu_supports("sse4.1")) {
    range_sse4(x, size, i, min_val, max_val);
  }
#endif

  // include the remainder not used during vectorised operations
  for ( ; i < size; i++) {
    min_val = std::min(min_val, x[i]);
    max_val = std::max(max_val, x[i]);
  }
  return {min_val, max_val};
}

} // namespace bgen
