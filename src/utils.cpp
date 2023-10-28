
#include <array>
#include <algorithm>
#include <iostream>

#include "utils.h"

#if defined(__x86_64__)
  #include <immintrin.h>
#endif

namespace bgen {

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
// We need to load the data into vector registers, but we need to convert the
// 8-bit vectors to 32-bit to avoid overflow. 32-bit uints should be sufficient,
// since this sums ploidy states, which can be at most 255 per person, so this
// allows at least 134 million individuals (((2 ** 31) * 16) / 255).
//
/// @param x array of floats
/// @param size size of array
/// @returns sum of array
std::uint64_t fast_ploidy_sum(std::uint8_t * x, std::uint32_t & size) {
  std::cout << "summing ploidy" << std::endl;
  size_t i = 0;
  std::uint64_t total = 0;

#if defined(__x86_64__)
  if (__builtin_cpu_supports("avx2")) {
    std::cout << " - avx2 ploidy" << std::endl;
    std::uint32_t arr[8];
    __m128i initial;
    __m256i _vals1, _vals2;
    __m256i _sum1 = _mm256_set_epi32(0, 0, 0, 0, 0, 0, 0, 0);
    __m256i _sum2 = _mm256_set_epi32(0, 0, 0, 0, 0, 0, 0, 0);
    for (; i + 16 < size; i += 16) {
      // load data and convert to 32-bit uints
      initial = _mm_loadu_si128((const __m128i*) &x[i]);
      _vals1 = _mm256_cvtepu8_epi32(initial);
      _vals2 = _mm256_cvtepu8_epi32(_mm_bsrli_si128(initial, 8));

      _sum1 = _mm256_add_epi32(_sum1, _vals1);
      _sum2 = _mm256_add_epi32(_sum2, _vals2);
    }
    _mm256_storeu_si256((__m256i*) &arr[0], _sum1);
    total += arr[0] + arr[1] + arr[2] + arr[3] + arr[4] + arr[5] + arr[6] + arr[7];
    _mm256_storeu_si256((__m256i*) &arr[0], _sum2);
    total += arr[0] + arr[1] + arr[2] + arr[3] + arr[4] + arr[5] + arr[6] + arr[7];
  }
#endif

  std::cout << " - completing ploidy, n=" << i << std::endl;
  // include the remainder not used during vectorised sum
  for ( ; i < size; i++) {
    total += x[i];
  }
  std::cout << " - ploidy complete, n=" << i << std::endl;
  return total;
}

// get min and max of the ploidy values in one pass
Range fast_range(std::uint8_t * x, std::uint32_t & size) {
  std::uint8_t min_val = 255;
  std::uint8_t max_val = 0;
  size_t i = 0;
  std::cout << "ploidy range" << std::endl;

#if defined(__x86_64__)
  if (__builtin_cpu_supports("avx2")) {
    std::cout << " - avx2 ploidy range" << std::endl;
    std::array<std::uint8_t, 32> arr;
    __m256i values;
    __m256i _mins = _mm256_set_epi8(-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                                    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                                    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
    __m256i _maxs = _mm256_set_epi8(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0);
    for (; i + 32 < size; i += 32) {
      // load data and convert to 32-bit uints
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
#endif

  std::cout << " - completing ploidy range, n=" << i << std::endl;
  // include the remainder not used during vectorised operations
  for ( ; i < size; i++) {
    min_val = std::min(min_val, x[i]);
    max_val = std::max(max_val, x[i]);
  }
  std::cout << " - ploidy range complete, n=" << i << std::endl;
  return {min_val, max_val};
}

} // namespace bgen
