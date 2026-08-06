
#include <array>
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

#include "utils.h"

#if defined(__x86_64__)
  #include <immintrin.h>
#elif defined(__aarch64__)
  #include <arm_neon.h>
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

// scan 32 bytes at a time for the missingness flag. Compiled for AVX2 in isolation, so
// this is only ever entered after a runtime check for AVX2 support.
BGEN_TARGET_AVX2
static void missing_scan_avx2(const char * data, std::uint32_t & size, std::uint32_t & x,
                              std::vector<std::uint32_t> & missing) {
  for (; x + 32 <= size; x += 32) {
    __m256i values = _mm256_loadu_si256((const __m256i*) &data[x]);
    // the missingness flag is the top bit of each byte, which is what movemask gathers,
    // giving one bit per sample and zero when none of the 32 are missing
    std::uint32_t bits = (std::uint32_t) _mm256_movemask_epi8(values);
    while (bits) {
      missing.push_back(x + (std::uint32_t) __builtin_ctz(bits));
      bits &= bits - 1;
    }
  }
}

// sum dosages into four 64-bit lanes, skipping missing samples. Compiled for AVX2 in
// isolation, so this is only ever entered after a runtime check for AVX2 support.
BGEN_TARGET_AVX2
static double dosage_sum_avx2(float * dose, std::uint32_t & size, std::uint32_t & i,
                              std::uint64_t & called) {
  __m256d acc = _mm256_setzero_pd();
  __m256i counts = _mm256_setzero_si256();
  const __m256i one = _mm256_set1_epi64x(1);
  for (; i + 4 <= size; i += 4) {
    __m256d values = _mm256_cvtps_pd(_mm_loadu_ps(&dose[i]));
    // a NaN compares unequal to itself, so this picks out the samples with a dosage
    __m256d ok = _mm256_cmp_pd(values, values, _CMP_EQ_OQ);
    acc = _mm256_add_pd(acc, _mm256_and_pd(values, ok));
    counts = _mm256_add_epi64(counts, _mm256_and_si256(_mm256_castpd_si256(ok), one));
  }
  double lanes[4];
  _mm256_storeu_pd(lanes, acc);
  std::uint64_t cl[4];
  _mm256_storeu_si256((__m256i*) &cl[0], counts);
  called += cl[0] + cl[1] + cl[2] + cl[3];
  return (lanes[0] + lanes[1]) + (lanes[2] + lanes[3]);
}

// sum every stride'th dosage, skipping missing samples. Compiled for AVX2 in isolation, so
// this is only ever entered after a runtime check for AVX2 support.
//
// The samples wanted are stride apart, so a plain load cannot reach them and a gather does
// the addressing instead. Eight indices go in and eight dosages come out, which replaces
// eight loads, eight NaN branches and eight dependent adds.
BGEN_TARGET_AVX2
static double strided_sum_avx2(const float * dose, std::uint32_t count,
                               std::uint32_t start, std::uint32_t stride,
                               std::uint32_t & i, std::uint64_t & called) {
  const int step = (int) (stride * 8);
  __m256i offs = _mm256_setr_epi32(0, (int) stride, (int) stride * 2, (int) stride * 3,
                                   (int) stride * 4, (int) stride * 5, (int) stride * 6,
                                   (int) stride * 7);
  const __m256i steps = _mm256_set1_epi32(step);
  __m256d acc = _mm256_setzero_pd();
  __m256i counts = _mm256_setzero_si256();
  const __m256i one = _mm256_set1_epi64x(1);

  for (; i + 8 <= count; i += 8) {
    __m256 values = _mm256_i32gather_ps(&dose[start], offs, 4);
    offs = _mm256_add_epi32(offs, steps);
    // widen in two halves, since the accumulator is double to match the scalar loop
    __m256d lo = _mm256_cvtps_pd(_mm256_castps256_ps128(values));
    __m256d hi = _mm256_cvtps_pd(_mm256_extractf128_ps(values, 1));
    // a NaN compares unequal to itself, so this picks out the samples with a dosage
    __m256d ok_lo = _mm256_cmp_pd(lo, lo, _CMP_EQ_OQ);
    __m256d ok_hi = _mm256_cmp_pd(hi, hi, _CMP_EQ_OQ);
    acc = _mm256_add_pd(acc, _mm256_and_pd(lo, ok_lo));
    acc = _mm256_add_pd(acc, _mm256_and_pd(hi, ok_hi));
    counts = _mm256_add_epi64(counts, _mm256_and_si256(_mm256_castpd_si256(ok_lo), one));
    counts = _mm256_add_epi64(counts, _mm256_and_si256(_mm256_castpd_si256(ok_hi), one));
  }
  double lanes[4];
  _mm256_storeu_pd(lanes, acc);
  std::uint64_t cl[4];
  _mm256_storeu_si256((__m256i*) &cl[0], counts);
  called += cl[0] + cl[1] + cl[2] + cl[3];
  return (lanes[0] + lanes[1]) + (lanes[2] + lanes[3]);
}

#elif defined(__aarch64__)

// sum ploidy states with NEON. NEON is part of the aarch64 baseline, so there is no
// runtime check to make before entering any of these.
static std::uint64_t ploidy_sum_neon(std::uint8_t * x, std::uint32_t & size,
                                     std::uint32_t & i) {
  uint64x2_t acc = vdupq_n_u64(0);
  for (; i + 16 <= size; i += 16) {
    // widen pairwise from 8 to 16 to 32 bits, then accumulate into 64-bit lanes, so
    // nothing can overflow whatever the ploidy values or the cohort size
    acc = vpadalq_u32(acc, vpaddlq_u16(vpaddlq_u8(vld1q_u8(&x[i]))));
  }
  return vgetq_lane_u64(acc, 0) + vgetq_lane_u64(acc, 1);
}

// get min and max of the ploidy values with NEON
static void range_neon(std::uint8_t * x, std::uint32_t & size, size_t & i,
                       std::uint8_t & min_val, std::uint8_t & max_val) {
  if (i + 16 > size) {
    return;
  }
  uint8x16_t mins = vdupq_n_u8(255);
  uint8x16_t maxs = vdupq_n_u8(0);
  for (; i + 16 <= size; i += 16) {
    uint8x16_t values = vld1q_u8(&x[i]);
    mins = vminq_u8(mins, values);
    maxs = vmaxq_u8(maxs, values);
  }
  // reduce across the lanes in one instruction each
  min_val = std::min(min_val, vminvq_u8(mins));
  max_val = std::max(max_val, vmaxvq_u8(maxs));
}

// scan 16 bytes at a time for the missingness flag.
//
// NEON has no equivalent of movemask, so the flags cannot be turned into a bitmask to
// walk. Instead the maximum across the masked lanes says whether any of the 16 samples is
// missing, which skips a clear block in one test. Missing samples are rare, so that is
// the case worth making fast, and a block that does flag something falls back to checking
// its bytes.
static void missing_scan_neon(const char * data, std::uint32_t & size, std::uint32_t & x,
                              std::vector<std::uint32_t> & missing) {
  const std::uint8_t * bytes = reinterpret_cast<const std::uint8_t *>(data);
  const uint8x16_t flag = vdupq_n_u8(0x80);
  for (; x + 16 <= size; x += 16) {
    uint8x16_t values = vld1q_u8(&bytes[x]);
    if (vmaxvq_u8(vandq_u8(values, flag)) == 0) {
      continue;
    }
    for (std::uint32_t k = 0; k < 16; k++) {
      if (bytes[x + k] & 0x80) {
        missing.push_back(x + k);
      }
    }
  }
}

// sum dosages into 64-bit lanes, skipping missing samples.
//
// The lanes are grouped and summed in the same order as the AVX2 path, so both give the
// same total for the same input.
static double dosage_sum_neon(float * dose, std::uint32_t & size, std::uint32_t & i,
                              std::uint64_t & called) {
  float64x2_t acc_lo = vdupq_n_f64(0.0);
  float64x2_t acc_hi = vdupq_n_f64(0.0);
  uint64x2_t counts_lo = vdupq_n_u64(0);
  uint64x2_t counts_hi = vdupq_n_u64(0);
  const uint64x2_t one = vdupq_n_u64(1);
  for (; i + 4 <= size; i += 4) {
    float32x4_t values = vld1q_f32(&dose[i]);
    // widen to double, since the accumulator has to match the scalar loop
    float64x2_t lo = vcvt_f64_f32(vget_low_f32(values));
    float64x2_t hi = vcvt_high_f64_f32(values);
    // a NaN compares unequal to itself, so this picks out the samples with a dosage
    uint64x2_t ok_lo = vceqq_f64(lo, lo);
    uint64x2_t ok_hi = vceqq_f64(hi, hi);
    acc_lo = vaddq_f64(acc_lo, vreinterpretq_f64_u64(
                                   vandq_u64(vreinterpretq_u64_f64(lo), ok_lo)));
    acc_hi = vaddq_f64(acc_hi, vreinterpretq_f64_u64(
                                   vandq_u64(vreinterpretq_u64_f64(hi), ok_hi)));
    counts_lo = vaddq_u64(counts_lo, vandq_u64(ok_lo, one));
    counts_hi = vaddq_u64(counts_hi, vandq_u64(ok_hi, one));
  }
  called += vgetq_lane_u64(counts_lo, 0) + vgetq_lane_u64(counts_lo, 1)
          + vgetq_lane_u64(counts_hi, 0) + vgetq_lane_u64(counts_hi, 1);
  double lanes[4] = {vgetq_lane_f64(acc_lo, 0), vgetq_lane_f64(acc_lo, 1),
                     vgetq_lane_f64(acc_hi, 0), vgetq_lane_f64(acc_hi, 1)};
  return (lanes[0] + lanes[1]) + (lanes[2] + lanes[3]);
}

#endif

/// @brief binomial coefficient C(n, k)
///
/// This counts the genotypes a sample can have, so it sizes the probability arrays and
/// the byte count a variant is checked against. Both arguments come from the file, built
/// from the allele count and the ploidy, so they can exceed any real variant.
///
/// @param n total items to choose from
/// @param k number of items chosen
/// @return C(n, k)
std::uint32_t n_choose_k(int n, int k) {
  if ((n < 0) || (k < 0) || (k > n)) {
    throw std::invalid_argument("cannot count the ways to choose " + std::to_string(k) +
                                " from " + std::to_string(n));
  }

  // Since C(n, k) = C(n, n-k)
  if ( k > n - k ) {
    k = n - k;
  }

  std::uint64_t res = 1;
  for (std::uint32_t i = 0; i < (std::uint32_t) k; ++i) {
    // throw error if next round would exceed integer limit
    if (res > UINT64_MAX / (std::uint64_t) (n - i)) {
      throw std::invalid_argument("variant needs more than 2^32 probabilities");
    }
    res *= (std::uint64_t) (n - i);
    res /= (i + 1);
  }

  if (res > UINT32_MAX) {
    throw std::invalid_argument("variant needs more than 2^32 probabilities");
  }
  return (std::uint32_t) res;
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
#elif defined(__aarch64__)
  total += ploidy_sum_neon(x, size, i);
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
#elif defined(__aarch64__)
  range_neon(x, size, i, min_val, max_val);
#endif

  // include the remainder not used during vectorised operations
  for ( ; i < size; i++) {
    min_val = std::min(min_val, x[i]);
    max_val = std::max(max_val, x[i]);
  }
  return {min_val, max_val};
}

/// @brief record the samples flagged as missing, in order
///
/// Layout 2 stores one byte per sample, with the top bit set when the sample is missing and
/// the low bits holding its ploidy. A variant with constant ploidy only needs the
/// missingness out of that, so this scans for the flag without touching the ploidy.
///
/// Missing samples are rare, so the useful case is a stretch with none, which this clears
/// 32 bytes at a time. A block that does flag something gets its indices straight out of
/// the movemask result, so the data is only read once either way.
///
/// @param data the ploidy byte per sample, at the start of the variant's block
/// @param size number of samples
/// @param missing indices of missing samples get appended, in increasing order
void fast_missing_scan(const char * data, std::uint32_t size,
                       std::vector<std::uint32_t> & missing) {
  std::uint32_t x = 0;

#if defined(__x86_64__)
  if (__builtin_cpu_supports("avx2")) {
    missing_scan_avx2(data, size, x, missing);
  }
#elif defined(__aarch64__)
  missing_scan_neon(data, size, x, missing);
#endif

  // include the remainder not used during vectorised operations
  for ( ; x < size; x++) {
    if (data[x] & 0x80) {
      missing.push_back(x);
    }
  }
}

/// @brief sum dosages and count the callable samples, skipping missing ones
///
/// find_minor_allele needs the dosage summed over the cohort, which it divides by the
/// alleles seen. Written as a plain loop this cannot vectorise, since accumulating floats
/// into a double is order dependent, and a branch per sample skips the missing ones.
///
/// Reordering the sum is safe here. Every dosage is a multiple of 1/maxval for the bit
/// depth, so the totals come out bit identical to the ordered sum for cohorts far larger
/// than any real one, and the caller only compares the frequency against 0.5. A missing
/// sample is marked with NaN and is dropped by adding zero rather than by branching, since
/// a NaN compares unequal to itself.
///
/// @param dose dosages, with missing samples already set to nan
/// @param size number of samples
/// @param n_called set to the number of samples which were not missing
/// @return sum of the dosages of the callable samples
double fast_dosage_sum(float * dose, std::uint32_t size, std::uint64_t & n_called) {
  double total = 0;
  std::uint64_t called = 0;
  std::uint32_t i = 0;

#if defined(__x86_64__)
  if (__builtin_cpu_supports("avx2")) {
    total += dosage_sum_avx2(dose, size, i, called);
  }
#elif defined(__aarch64__)
  total += dosage_sum_neon(dose, size, i, called);
#endif

  // include the remainder not used during vectorised operations
  for ( ; i < size; i++) {
    if (std::isnan(dose[i])) {
      continue;
    }
    total += dose[i];
    called++;
  }
  n_called = called;
  return total;
}

/// @brief sum every stride'th dosage, counting the callable samples
///
/// find_minor_allele samples the cohort in strides before it commits to a full pass, so
/// that a cohort ordered by genotype cannot fool the early exit. This does one such stride.
///
/// The reordering is safe for the same reason as fast_dosage_sum, and the result is bit
/// identical to summing the same samples in index order.
///
/// The gather addresses samples as signed 32-bit byte offsets from the start of the walk,
/// so it is only used while the span stays inside that range. Beyond it the scalar loop
/// runs instead, which cannot happen for a cohort that fits in memory but is checked
/// rather than assumed.
///
/// There is no aarch64 path here, unlike the other helpers in this file. NEON has no
/// gather, so the four lanes have to be loaded one at a time and inserted, which measured
/// no faster than the scalar loop it would replace.
///
/// @param dose dosages, with missing samples already set to nan
/// @param n_samples number of samples in the cohort
/// @param start index of the first sample to visit
/// @param stride gap between visited samples
/// @param n_checked incremented by the number of visited samples which were not missing
/// @return sum of the dosages of the visited callable samples
double fast_strided_sum(const float * dose, std::uint32_t n_samples, std::uint32_t start,
                        std::uint32_t stride, std::uint64_t & n_checked) {
  double total = 0;
  if (start >= n_samples) {
    return total;
  }
  // how many samples this walk visits, counting the one it starts on
  std::uint32_t count = (n_samples - 1 - start) / stride + 1;
  std::uint32_t i = 0;
  std::uint64_t called = 0;

#if defined(__x86_64__)
  // the largest byte offset the gather forms is stride * 4 * (count - 1) plus the seven
  // lanes of the last block, so cap on the span in bytes staying positive as an int32
  std::uint64_t span = (std::uint64_t) stride * 4 * count;
  if (__builtin_cpu_supports("avx2") && span < 0x7FFFFFFF) {
    total += strided_sum_avx2(dose, count, start, stride, i, called);
  }
#endif

  // include the remainder not used during vectorised operations
  for ( ; i < count; i++) {
    float v = dose[start + (std::uint64_t) i * stride];
    if (std::isnan(v)) {
      continue;
    }
    total += v;
    called++;
  }
  n_checked += called;
  return total;
}

} // namespace bgen
