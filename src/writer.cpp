
#include <array>
#include <cstring>
#include <cmath>
#include <fstream>
#include <stdexcept>
#include <vector>

#include <iostream>

#include "zstd/lib/zstd.h"
#include <zlib.h>

#if defined(__x86_64__)
  #include <immintrin.h>
#endif

#if defined(__aarch64__)
  #include <arm_neon.h>
#endif

#include "writer.h"
#include "genotypes.h"
#include "utils.h"

namespace bgen {

// append raw bytes to a buffer destined for the file handle
static void append_bytes(std::vector<char> &buf, const void *data, std::size_t size) {
  const char * bytes = reinterpret_cast<const char *>(data);
  buf.insert(buf.end(), bytes, bytes + size);
}

/// @brief the current write position, or a clear error if there is not one
///
/// These positions become the header's variant offset and the .bgi index, so a wrong one
/// silently misplaces every read. tellp() returns -1 on an already failed stream without
/// tripping the exception mask, and throws an unhelpful message on a pipe
static std::uint64_t current_position(std::ofstream &handle) {
  std::streamoff pos = -1;
  if (!handle.fail()) {
    try {
      pos = handle.tellp();
    } catch (const std::exception &) {
      // the mask turned the failed seek into a throw. Fall through to report it in
      // terms of what the caller can act on
      pos = -1;
    }
  }
  if (pos < 0) {
    throw std::ios_base::failure("cannot find the current position in the bgen. A bgen "
                                 "has to be written to a file that can seek, rather than "
                                 "to a pipe, and no more can be written after a failed "
                                 "write");
  }
  return (std::uint64_t) pos;
}

// write a 32-bit value at a given file offset
static void write_at_offset(std::ofstream &handle, std::uint32_t &val, std::uint32_t offset=0) {
  std::uint64_t orig_pos = current_position(handle);
  handle.seekp(offset);
  handle.write(reinterpret_cast<char *>(&val), 4);
  handle.seekp(orig_pos);
}

static void write_variants_offset(std::ofstream &handle, std::uint32_t &offset) {
  write_at_offset(handle, offset, 0);
}

static void write_nvariants(std::ofstream &handle, std::uint32_t &offset, std::uint32_t &n_variants) {
  write_at_offset(handle, n_variants, offset);
}

// finish the file off by writing the number of variants, and where the variant data
// starts. This can throw on a full disk, so call it rather than relying on the destructor
void CppBgenWriter::close() {
  if (closed) { return; }
  closed = true;
  write_variants_offset(handle, variant_data_offset);
  write_nvariants(handle, nvars_offset, n_variants);
  handle.close();
}

// a destructor must not throw, so any error while finishing the file can only
// be swallowed here. Call close() directly to have those errors reported.
CppBgenWriter::~CppBgenWriter() {
  try {
    close();
  } catch (...) {
  }
}

void CppBgenWriter::write_header(std::string &free_data,
                              std::vector<std::string> &samples) {
  std::uint32_t header_len = 20 + free_data.size();
  variant_data_offset = header_len;
  write_variants_offset(handle, variant_data_offset);
  handle.seekp(4);

  // figure out the length of the header block
  handle.write(reinterpret_cast<char *>(&header_len), 4);

  // write zero variants for now, is fixed while closing
  handle.write(reinterpret_cast<char *>(&n_variants), 4);
  handle.write(reinterpret_cast<char *>(&n_samples), 4);
  handle << "bgen";
  handle << free_data;

  // check and write flags
  if (compression > 2) {
    throw std::invalid_argument("compression flag must be 0, 1, or 2");
  }
  if ((layout < 1) || layout > 2) {
    throw std::invalid_argument("layout flag must be 1, or 2");
  }
  std::uint32_t sample_id_flag = samples.size() > 0;
  std::uint32_t flags = 0;
  flags |= compression;
  flags |= layout << 2;
  flags |= sample_id_flag << 31;
  handle.write(reinterpret_cast<char *>(&flags), 4);
  }

/// @brief check a string fits the 16-bit length field the bgen stores it behind
///
/// The length would otherwise wrap while the whole string is still written, leaving a
/// file which looks consistent, but whose fields all read from the wrong offsets
static void check_length_fits(const std::string &value, const char *field) {
  if (value.size() > UINT16_MAX) {
    throw std::invalid_argument(std::string(field) + " is too long for a bgen - "
                                "it is " + std::to_string(value.size()) +
                                " characters, but the maximum is " +
                                std::to_string(UINT16_MAX));
  }
}

void CppBgenWriter::add_samples(std::vector<std::string> &samples) {
  if (samples.size() == 0) { return; }

  if (samples.size() != n_samples) {
    throw std::invalid_argument("samples vector length doesn't match the sample count in file");
  }

  // count the number of characters across all sample IDs
  std::uint64_t nchars = 0;
  for (auto &x : samples) {
    check_length_fits(x, "sample ID");
    nchars += x.size();
  }

  // write the length of the sample ID block, and number of sample IDs
  std::uint64_t block_len = 8 + 2 * (std::uint64_t) samples.size() + nchars;
  if (block_len > UINT32_MAX) {
    // the block length is a four byte field, so a block this big would wrap and
    // leave the variants starting at the wrong offset
    throw std::invalid_argument("sample IDs need " + std::to_string(block_len) +
                                " bytes, which overflows the bgen sample block "
                                "length field");
  }
  std::uint32_t samples_len = (std::uint32_t) block_len;
  handle.write(reinterpret_cast<char *>(&samples_len), 4);
  std::uint32_t size = samples.size();
  handle.write(reinterpret_cast<char *>(&size), 4);

  // write each sample ID to the file, preceeded by each ID length
  std::uint16_t id_size;
  for (auto &x : samples) {
    id_size = x.size();
    handle.write(reinterpret_cast<char *>(&id_size), 2);
    handle << x;
  }
  // the offset is a four byte field, so check it in 64 bits first, since the check above
  // bounds the sample block alone, and the header pushes the two together over
  std::uint64_t data_offset = current_position(handle) - 4;
  if (data_offset > UINT32_MAX) {
    throw std::invalid_argument("the header and sample IDs take up " +
                                std::to_string(data_offset) + " bytes, which overflows "
                                "the bgen variant data offset field");
  }
  variant_data_offset = (std::uint32_t) data_offset;
  write_variants_offset(handle, variant_data_offset);
}

std::uint64_t CppBgenWriter::write_variant_header(std::string &varid,
                                                  std::string &rsid,
                                                  std::string &chrom,
                                                  std::uint32_t &pos,
                                                  std::vector<std::string> &alleles,
                                                  std::uint32_t _n_samples) {
  std::uint64_t var_offset = current_position(handle);
  if (_n_samples != n_samples) {
    throw std::invalid_argument("number of samples doesn't match sample count in file");
  }
  if ((layout == 1) && alleles.size() != 2) {
    throw std::invalid_argument("layout 1 requires exactly two alleles.");
  }
  // check these before anything is written, so a rejected variant leaves the
  // file as it was rather than partly extended by a variant we then refuse
  check_length_fits(varid, "variant ID");
  check_length_fits(rsid, "rsID");
  check_length_fits(chrom, "chromosome");
  if (alleles.size() > UINT16_MAX) {
    // the allele count is a two byte field, so more than this would wrap and
    // leave the reader looking for the wrong number of alleles
    throw std::invalid_argument("variant has " + std::to_string(alleles.size()) +
                                " alleles, but the maximum is " +
                                std::to_string(UINT16_MAX));
  }
  n_variants += 1;
  if (layout == 1) {
    handle.write(reinterpret_cast<char *>(&_n_samples), 4);
    // handle << _n_samples;
  }
  std::uint16_t tmp;
  tmp = varid.size();
  handle.write(reinterpret_cast<char *>(&tmp), 2);
  handle << varid;
  tmp = rsid.size();
  handle.write(reinterpret_cast<char *>(&tmp), 2);
  handle << rsid;
  tmp = chrom.size();
  handle.write(reinterpret_cast<char *>(&tmp), 2);
  handle << chrom;
  handle.write(reinterpret_cast<char *>(&pos), 4);
  
  if (layout != 1) {
    std::uint16_t n_alleles = alleles.size();
    handle.write(reinterpret_cast<char *>(&n_alleles), 2);
  }

  std::uint32_t allele_size;
  for (auto &x : alleles) {
    allele_size = x.size();
    handle.write(reinterpret_cast<char *>(&allele_size), 4);
    handle << x;
  }
  handle.flush();
  return var_offset;
}

std::uint64_t CppBgenWriter::write_variant_direct(std::vector<std::uint8_t> & data) {
  std::uint64_t var_offset = current_position(handle);
  n_variants += 1;
  // write() rather than a std::ostreambuf_iterator, which reports a failed write only in
  // its own state, leaving the stream, and so the exception mask, none the wiser
  handle.write(reinterpret_cast<char *>(data.data()), data.size());
  return var_offset;
}

// compress a char array with zlib
static void zlib_compress(char * input, int input_len, std::vector<char> &output) {
  z_stream strm;
  strm.zalloc = Z_NULL;
  strm.zfree = Z_NULL;
  strm.opaque = Z_NULL;

  strm.avail_in = input_len;      // size of input
  strm.next_in = (Bytef *) input; // input char array
  strm.avail_out = output.size();        // size of output
  strm.next_out = (Bytef *) &output[0]; // output char array

  int ret;
  deflateInit(&strm, 6);
  ret = deflate(&strm, Z_FINISH);
  deflateEnd(&strm);

  if (ret != Z_STREAM_END) {
    throw(std::invalid_argument("zlib compression encountered an error"));
  }

  output.resize(strm.total_out);
}

// compress a char array with zstd
static void zstd_compress(char *input, int input_len, std::vector<char> &output) {
  std::size_t total_out = ZSTD_compress(&output[0], output.size(), input, input_len, 3);

  if (ZSTD_isError(total_out)) {
    throw(std::invalid_argument("zstd compression encountered an error"));
  }

  output.resize(total_out);
}

/// Compress genotype data for a variant.
///
/// Compression is handled internally by either zlib_decompress, or zstd_decompress,
/// depending on compression scheme.
static std::vector<char> compress(std::vector<std::uint8_t> &uncompressed, std::uint32_t compression) {
  std::size_t bound;
  if (compression == 1) { // zlib
    bound = compressBound(uncompressed.size());
  } else { // zstd
    bound = ZSTD_compressBound(uncompressed.size());
  }
  std::vector<char> compressed(bound);
  if (compression == 1) { // zlib
    zlib_compress(reinterpret_cast<char *>(&uncompressed[0]), (int)uncompressed.size(), compressed);
  } else if (compression == 2) { // zstd
    zstd_compress(reinterpret_cast<char *>(&uncompressed[0]), (int)uncompressed.size(), compressed);
  }
  return compressed;
}

/// tolerance for genotype probabilities which sit just outside the legal range, since
/// they often arrive via the float32 the reader hands back, whose epsilon of 1.2e-07
/// accumulates across a sample. Values within this are clamped, further out is an error
const double PROB_TOLERANCE = 1e-6;

/// @brief report an unusable genotype probability, kept out of line so the checks below
/// stay small enough to inline into the encoding loops
static void raise_probability_error(double g) {
  if (std::isnan(g)) {
    throw std::invalid_argument("samples with any missing genotype must encode all "
                                "as missing (i.e. float(nan))");
  }
  throw std::invalid_argument("genotype probability must be between 0 and 1, not " +
                              std::to_string(g));
}

static void raise_cumulative_error(double cumulative) {
  throw std::invalid_argument("genotype probabilities for a sample sum to more "
                              "than 1: " + std::to_string(cumulative));
}

/// @brief check a genotype probability is in the range that layout 1 and 2 can encode
///
/// The negated range test also fails a nan, which here means a partly missing sample,
/// since fully missing ones are caught before this and encoded as zeroes
static inline void check_probability(double g) {
  if (!((g >= -PROB_TOLERANCE) && (g <= 1.0 + PROB_TOLERANCE))) {
    raise_probability_error(g);
  }
}

/// @brief check the probabilities stored for a sample do not sum above one
///
/// Layout 2 infers the last probability of a group as the remainder, so a group summing
/// above one cannot be encoded. Summing below one is left alone, as a partial distribution
static inline void check_cumulative(double cumulative) {
  if (cumulative > 1.0 + PROB_TOLERANCE) {
    raise_cumulative_error(cumulative);
  }
}

/// @brief round a non-negative double to the nearest integer, halves away from zero
///
/// std::round is an out of line libm call, and this sits in the innermost encode loop.
/// The truncation and subtraction are exact, so this avoids the std::floor(x + 0.5) trap.
/// Only defined for finite non-negative x below 2^53, which is what the callers pass
static inline std::uint64_t round_nonneg(double x) {
  std::uint64_t truncated = (std::uint64_t) x;
  return truncated + (x - (double) truncated >= 0.5);
}

/// @brief decide whether a sample's genotypes are all missing, rejecting a partial row
///
/// inline, since the encode paths call this per sample, for as few as three values. The
/// count stays branchless, as an early exit measured slower once half are missing
static inline bool missing_genotypes(double *genotypes, std::uint32_t size) {
  std::uint32_t nan_count = 0;
  for (std::uint32_t i=0; i<size; i++) {
    nan_count += std::isnan(genotypes[i]);
  }
  if ((nan_count > 0) && (nan_count < size)) {
    throw std::invalid_argument("samples with any missing genotype must encode all as missing (i.e. float(nan))");
  }
  return nan_count == size;
}

static std::vector<std::uint8_t> encode_layout1(
                    double *genotypes,
                    std::uint32_t geno_len) {
  // genotypes are encoded as 16-bit uints, so resize to n_genotypes * 2
  std::vector<std::uint8_t> encoded(geno_len * 2 + 8);

  std::uint32_t i = 0;
  std::int32_t scaled32;
  std::uint16_t scaled;
  bool missing;
  double g;
  double cumulative;
  for (std::uint32_t j=0; j < geno_len; j+=3) {
    missing = missing_genotypes(&genotypes[j], 3);
    cumulative = 0.0;
    for (std::uint32_t k=0; k<3; k++) {
      g = genotypes[j + k];
      if (missing) {
        g = 0;
      }
      // layout 1 could hold values up to 2.0, but is held to the same range as layout 2,
      // so both accept the same genotypes
      check_probability(g);
      cumulative += g;
      check_cumulative(cumulative);
      // check_probability allows a small negative, which round_nonneg is not defined
      // for, and which scales to zero anyway
      double product = g * 32768;
      scaled32 = (product > 0.0) ? (std::int32_t) round_nonneg(product) : 0;
      // check the value is in bounds
      if ((scaled32 < 0) || (scaled32 > 65535)) {
        throw std::invalid_argument("scaled genotype is out of bounds");
      }
      scaled = scaled32;
      std::memcpy(&encoded[i], &scaled, 2);
      i += 2;
    }
  }
  encoded.resize(geno_len * 2);
  return encoded;
}

/// @brief scale a running total of probabilities into the encoded integer range
///
/// Layout 2 infers a sample's last probability as max - sum(stored), so the stored values
/// must sum to at most max. Scaling each alone breaks that: at 8 bits two probabilities of
/// 0.5 both round to 128, which sums past 255. The bounds go on the double, before the
/// cast, which is only defined for values the integer type can represent
///
/// @param cumulative running total of the probabilities so far, including this one
/// @param factor scaling factor, the maximum encoded value for the bit depth. Passed by
///        value since it is only read, which also lets a caller hand over a constant
/// @param previous value this returned for the preceding probability, so that
///        clamping cannot make the current probability encode as negative
/// @return the running total scaled to the encoded integer range
static std::uint64_t scale_cumulative(double cumulative,
                                      double factor,
                                      std::uint64_t previous)
{
  double scaled = cumulative * factor;
  // the comparison order also catches nan, which must not reach the cast below
  if (!(scaled >= 0.0)) {
    return previous;
  }
  if (scaled > factor) {
    return (std::uint64_t) factor;
  }
  std::uint64_t value = round_nonneg(scaled);
  // previous is either zero or an earlier return of this function, so it is never
  // above factor, and returning it cannot exceed the range for the bit depth
  if (value < previous) {
    return previous;
  }
  return value;
}

/// @brief figure out the 64-bit pattern to insert an encoded genotype probability
/// @param value encoded genotype probability, already scaled to the bit depth
/// @param encoded vector of proabilties, set up as 8-bit. We pull a 64-bit slice
///        of this at the pointer position, in order to swap in the bits for the
///        encoded genotype at the correct offset
/// @param bit_remainder bit offset to place the encoded genotype at
/// @return data with probability inserted
static std::uint64_t emplace_probability(std::uint64_t value,
                                  std::uint8_t *encoded,
                                  std::uint32_t &bit_remainder)
{
  std::uint64_t window;

  // read the full 64 bits the caller writes back, so the bits already in the buffer
  // survive, and a 32 bit probability shifted across 39 bits keeps its top
  std::memcpy(&window, encoded, sizeof(window));
  window |= (value << bit_remainder);
  return window;
}

/// @brief precompute number of probabilities per sample for every ploidy
///
/// Only entries up to max_ploidy are filled, since higher ones cannot occur, and are not
/// even representable (12 alleles at a ploidy of 63 exceeds a 32-bit count)
///
/// @param n_alleles allele count for the variant
/// @param max_ploidy highest ploidy any sample has
/// @param phased whether the genotypes are phased
/// @return counts indexed by ploidy, valid up to max_ploidy
static std::array<std::uint32_t, 64> probs_per_ploidy(int n_alleles, int max_ploidy,
                                                     bool phased) {
  std::array<std::uint32_t, 64> counts{};
  for (int ploid = 0; ploid <= max_ploidy; ploid++) {
    int p = ploid;
    counts[ploid] = get_max_probs(p, n_alleles, phased);
  }
  return counts;
}

/// @brief encode a range of samples of the biallelic, unphased, ploidy 2, 8 bit shape
///
/// Split out so the vector paths below can hand back any batch they decline, which keeps
/// every rounding, clamping and error decision in one place. Encodes samples [first, last)
///
/// @param out first output byte of the variant, two per sample
/// @param flags per sample ploidy bytes, whose top bit marks a missing sample
/// @param first first sample to encode
/// @param last one past the last sample to encode
/// @param genotypes three probabilities per sample
static void encode_biallelic_8bit_range(std::uint8_t *out,
                     std::uint8_t *flags,
                     std::uint32_t first,
                     std::uint32_t last,
                     double *genotypes)
{
  const double factor = 255.0;
  for (std::uint32_t n=first; n < last; n++) {
    double *probs = &genotypes[3 * n];
    std::uint8_t *bytes = out + 2 * n;
    if (missing_genotypes(probs, 3)) {
      flags[n] |= 0x80;
      // as in encode_unphased, both bytes of a missing sample encode as zero and the
      // buffer already holds zeroes, so there is nothing to write
      continue;
    }
    // the probabilities for a sample sum to 1.0, so scale their running total
    check_probability(probs[0]);
    double cumulative = probs[0];
    check_cumulative(cumulative);
    std::uint64_t running = scale_cumulative(cumulative, factor, 0);
    bytes[0] = (std::uint8_t) running;

    check_probability(probs[1]);
    cumulative += probs[1];
    check_cumulative(cumulative);
    std::uint64_t previous = running;
    running = scale_cumulative(cumulative, factor, previous);
    bytes[1] = (std::uint8_t) (running - previous);

    // the third probability is inferred by the reader rather than stored, but it still
    // has to be checked, or a row summing above one would be written as something else
    check_probability(probs[2]);
    cumulative += probs[2];
    check_cumulative(cumulative);
  }
}

#if defined(__x86_64__)

/// @brief deinterleave 4 samples' worth of stride 3 doubles into one vector per probability
///
/// The genotypes arrive as p0,p1,p2 per sample, but the arithmetic wants all four samples'
/// p0 in one register. Three loads cover 12 doubles, and each output lane is then selected
/// from whichever load holds it: p0 comes from A[0], A[3], B[2], C[1], and so on.
BGEN_TARGET_AVX2
static inline void deinterleave3_avx2(const double *src, __m256d &p0, __m256d &p1,
                                      __m256d &p2) {
  #define BGEN_SEL(a, b, c, d) (((d) << 6) | ((c) << 4) | ((b) << 2) | (a))
  __m256d A = _mm256_loadu_pd(src);
  __m256d B = _mm256_loadu_pd(src + 4);
  __m256d C = _mm256_loadu_pd(src + 8);
  p0 = _mm256_blend_pd(
      _mm256_blend_pd(_mm256_permute4x64_pd(A, BGEN_SEL(0, 3, 0, 0)),
                      _mm256_permute4x64_pd(B, BGEN_SEL(0, 0, 2, 0)), 0x4),
      _mm256_permute4x64_pd(C, BGEN_SEL(0, 0, 0, 1)), 0x8);
  p1 = _mm256_blend_pd(
      _mm256_blend_pd(_mm256_permute4x64_pd(A, BGEN_SEL(1, 0, 0, 0)),
                      _mm256_permute4x64_pd(B, BGEN_SEL(0, 0, 3, 0)), 0x6),
      _mm256_permute4x64_pd(C, BGEN_SEL(0, 0, 0, 2)), 0x8);
  p2 = _mm256_blend_pd(
      _mm256_blend_pd(_mm256_permute4x64_pd(A, BGEN_SEL(2, 0, 0, 0)),
                      _mm256_permute4x64_pd(B, BGEN_SEL(0, 1, 0, 0)), 0x2),
      _mm256_permute4x64_pd(C, BGEN_SEL(0, 0, 0, 3)), 0xC);
  #undef BGEN_SEL
}

/// @brief round non-negative doubles half away from zero, matching round_nonneg
///
/// _mm256_round_pd's nearest mode is half to even, which differs from round_nonneg at every
/// exact .5, and probabilities that are multiples of 1/510 land there exactly. So floor and
/// then add one wherever the fraction reaches a half, which is what the scalar version does.
BGEN_TARGET_AVX2
static inline __m256d round_nonneg_avx2(__m256d x) {
  __m256d floored = _mm256_round_pd(x, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);
  __m256d frac = _mm256_sub_pd(x, floored);
  __m256d carry = _mm256_and_pd(_mm256_cmp_pd(frac, _mm256_set1_pd(0.5), _CMP_GE_OQ),
                                _mm256_set1_pd(1.0));
  return _mm256_add_pd(floored, carry);
}

/// @brief encode the biallelic 8 bit shape four samples at a time
///
/// Each sample is independent, so the scaling vectorises, but the checks throw, so a batch
/// with anything unencodable goes to the scalar range function, which reaches the same
/// check in the same order. Missing samples are handled here, being too common to defer
BGEN_TARGET_AVX2
static void encode_biallelic_8bit_avx2(std::uint8_t *out, std::uint8_t *flags,
                                       std::uint32_t n_samples, double *genotypes) {
  const __m256d lo_ok = _mm256_set1_pd(-PROB_TOLERANCE);
  const __m256d hi_ok = _mm256_set1_pd(1.0 + PROB_TOLERANCE);
  const __m256d factor = _mm256_set1_pd(255.0);
  const __m256d zero = _mm256_setzero_pd();

  std::uint32_t n = 0;
  for (; n + 4 <= n_samples; n += 4) {
    __m256d p0, p1, p2;
    deinterleave3_avx2(&genotypes[3 * n], p0, p1, p2);

    // the same range check_probability applies. Both comparisons are false for nan, so a
    // missing sample fails this too and is dealt with below
    __m256d ok = _mm256_and_pd(_mm256_cmp_pd(p0, lo_ok, _CMP_GE_OQ),
                               _mm256_cmp_pd(p0, hi_ok, _CMP_LE_OQ));
    ok = _mm256_and_pd(ok, _mm256_and_pd(_mm256_cmp_pd(p1, lo_ok, _CMP_GE_OQ),
                                         _mm256_cmp_pd(p1, hi_ok, _CMP_LE_OQ)));
    ok = _mm256_and_pd(ok, _mm256_and_pd(_mm256_cmp_pd(p2, lo_ok, _CMP_GE_OQ),
                                         _mm256_cmp_pd(p2, hi_ok, _CMP_LE_OQ)));
    __m256d c1 = _mm256_add_pd(p0, p1);
    __m256d total = _mm256_add_pd(c1, p2);
    // and the check_cumulative bound on both running totals. The first probability alone
    // cannot exceed it, having already passed the range check
    ok = _mm256_and_pd(ok, _mm256_cmp_pd(c1, hi_ok, _CMP_LE_OQ));
    ok = _mm256_and_pd(ok, _mm256_cmp_pd(total, hi_ok, _CMP_LE_OQ));

    int ok_mask = _mm256_movemask_pd(ok);
    int miss_mask = 0;
    if (ok_mask != 0xF) {
      // a lane can fail simply because the sample is missing, so identify those: nan is
      // the only value not equal to itself
      __m256d all_nan = _mm256_and_pd(
          _mm256_cmp_pd(p0, p0, _CMP_UNORD_Q),
          _mm256_and_pd(_mm256_cmp_pd(p1, p1, _CMP_UNORD_Q),
                        _mm256_cmp_pd(p2, p2, _CMP_UNORD_Q)));
      miss_mask = _mm256_movemask_pd(all_nan);
      if ((ok_mask | miss_mask) != 0xF) {
        // something is partially missing or out of range, so let the scalar loop throw
        encode_biallelic_8bit_range(out, flags, n, n + 4, genotypes);
        continue;
      }
      for (int lane = 0; lane < 4; lane++) {
        if (miss_mask & (1 << lane)) {
          flags[n + lane] |= 0x80;
        }
      }
      // substitute zero for a missing lane, as the scalar path does, rather than relying
      // on how each architecture happens to convert the nan
      p0 = _mm256_blendv_pd(p0, zero, all_nan);
      c1 = _mm256_blendv_pd(c1, zero, all_nan);
    }

    // bound the doubles before any cast, exactly as scale_cumulative does
    __m256d s0 = _mm256_min_pd(_mm256_max_pd(_mm256_mul_pd(p0, factor), zero), factor);
    __m256d s1 = _mm256_min_pd(_mm256_max_pd(_mm256_mul_pd(c1, factor), zero), factor);
    __m256d v0 = round_nonneg_avx2(s0);
    // the monotonic clamp, so the stored difference below cannot go negative
    __m256d v1 = _mm256_max_pd(round_nonneg_avx2(s1), v0);

    __m128i i0 = _mm256_cvttpd_epi32(v0);
    __m128i i1 = _mm256_cvttpd_epi32(_mm256_sub_pd(v1, v0));
    // narrow to bytes, then interleave so each sample's pair ends up adjacent. The
    // narrowing saturates, so a missing lane writes zeroes over zeroes
    __m128i packed = _mm_packus_epi16(_mm_packus_epi32(i0, i1), _mm_setzero_si128());
    __m128i pairs = _mm_unpacklo_epi8(packed, _mm_srli_si128(packed, 4));
    std::memcpy(out + 2 * n, &pairs, 8);
  }
  // however many samples are left over after the last full batch
  encode_biallelic_8bit_range(out, flags, n, n_samples, genotypes);
}
#endif

#if defined(__aarch64__)

/// @brief round non-negative doubles half away from zero, matching round_nonneg
///
/// vrndnq_f64 rounds halves to even, which differs from round_nonneg at every exact .5, so
/// floor and add one wherever the fraction reaches a half, as the scalar version does.
static inline float64x2_t round_nonneg_neon(float64x2_t x) {
  float64x2_t floored = vrndmq_f64(x);
  float64x2_t frac = vsubq_f64(x, floored);
  // a true comparison lane is all ones, so masking 1.0 with it gives the increment
  uint64x2_t carry = vandq_u64(vcgeq_f64(frac, vdupq_n_f64(0.5)),
                               vreinterpretq_u64_f64(vdupq_n_f64(1.0)));
  return vaddq_f64(floored, vreinterpretq_f64_u64(carry));
}

/// @brief encode the biallelic 8 bit shape two samples at a time
///
/// The same structure as the avx2 version, and for the same reasons, but a NEON double
/// vector holds two lanes rather than four. In exchange the stride 3 layout costs nothing
/// to unpack, since vld3q_f64 deinterleaves as it loads.
static void encode_biallelic_8bit_neon(std::uint8_t *out, std::uint8_t *flags,
                                       std::uint32_t n_samples, double *genotypes) {
  const float64x2_t lo_ok = vdupq_n_f64(-PROB_TOLERANCE);
  const float64x2_t hi_ok = vdupq_n_f64(1.0 + PROB_TOLERANCE);
  const float64x2_t factor = vdupq_n_f64(255.0);
  const float64x2_t zero = vdupq_n_f64(0.0);

  std::uint32_t n = 0;
  for (; n + 2 <= n_samples; n += 2) {
    float64x2x3_t loaded = vld3q_f64(&genotypes[3 * n]);
    float64x2_t p0 = loaded.val[0], p1 = loaded.val[1], p2 = loaded.val[2];

    uint64x2_t ok = vandq_u64(vcgeq_f64(p0, lo_ok), vcleq_f64(p0, hi_ok));
    ok = vandq_u64(ok, vandq_u64(vcgeq_f64(p1, lo_ok), vcleq_f64(p1, hi_ok)));
    ok = vandq_u64(ok, vandq_u64(vcgeq_f64(p2, lo_ok), vcleq_f64(p2, hi_ok)));
    float64x2_t c1 = vaddq_f64(p0, p1);
    float64x2_t total = vaddq_f64(c1, p2);
    ok = vandq_u64(ok, vcleq_f64(c1, hi_ok));
    ok = vandq_u64(ok, vcleq_f64(total, hi_ok));

    // vminvq across the mask reinterpreted as 32-bit lanes is zero unless every lane passed
    if (vminvq_u32(vreinterpretq_u32_u64(ok)) == 0) {
      // vceqq_f64(x, x) is false only for nan, so a lane whose three probabilities are all
      // nan has every bit of this clear
      uint64x2_t present = vorrq_u64(vorrq_u64(vceqq_f64(p0, p0), vceqq_f64(p1, p1)),
                                     vceqq_f64(p2, p2));
      uint64x2_t all_nan = vceqzq_u64(present);
      // both lanes have to be either encodable or wholly missing to stay here
      if (vminvq_u32(vreinterpretq_u32_u64(vorrq_u64(ok, all_nan))) == 0) {
        encode_biallelic_8bit_range(out, flags, n, n + 2, genotypes);
        continue;
      }
      if (vgetq_lane_u64(all_nan, 0) != 0) {
        flags[n] |= 0x80;
      }
      if (vgetq_lane_u64(all_nan, 1) != 0) {
        flags[n + 1] |= 0x80;
      }
      // as in the avx2 path, substitute zero so a nan lane stores the zeroes already in
      // the buffer, rather than depending on how the conversion handles a nan
      p0 = vbslq_f64(all_nan, zero, p0);
      c1 = vbslq_f64(all_nan, zero, c1);
    }

    float64x2_t s0 = vminq_f64(vmaxq_f64(vmulq_f64(p0, factor), zero), factor);
    float64x2_t s1 = vminq_f64(vmaxq_f64(vmulq_f64(c1, factor), zero), factor);
    float64x2_t v0 = round_nonneg_neon(s0);
    float64x2_t v1 = vmaxq_f64(round_nonneg_neon(s1), v0);

    // vcvtq_u64_f64 saturates rather than wrapping, so it plays the same part as the x86
    // packus: a difference that came out negative lands on zero
    uint64x2_t u0 = vcvtq_u64_f64(v0);
    uint64x2_t u1 = vcvtq_u64_f64(vsubq_f64(v1, v0));
    out[2 * n] = (std::uint8_t) vgetq_lane_u64(u0, 0);
    out[2 * n + 1] = (std::uint8_t) vgetq_lane_u64(u1, 0);
    out[2 * n + 2] = (std::uint8_t) vgetq_lane_u64(u0, 1);
    out[2 * n + 3] = (std::uint8_t) vgetq_lane_u64(u1, 1);
  }
  encode_biallelic_8bit_range(out, flags, n, n_samples, genotypes);
}
#endif

/// @brief encode the common case of biallelic, unphased, ploidy 2, at a bit depth of 8
///
/// Almost all real data has this shape, which drops the generic loop's per sample
/// overhead: three probabilities in, two bytes out, a factor of 255, and byte aligned
/// output. scale_cumulative still decides the rounding, so only the loop is specialised
///
/// @param encoded buffer to write into, already zero filled
/// @param genotype_offset where this variant's probability bytes start
/// @param ploidy_offset where the per sample ploidy bytes start
/// @param n_samples number of samples
/// @param genotypes three probabilities per sample
/// @return offset one past the last byte written
static std::uint32_t encode_biallelic_8bit(std::vector<std::uint8_t> &encoded,
                     std::uint32_t genotype_offset,
                     std::uint32_t ploidy_offset,
                     std::uint32_t n_samples,
                     double *genotypes)
{
  std::uint8_t *out = &encoded[genotype_offset];
  std::uint8_t *flags = &encoded[ploidy_offset];
#if defined(__x86_64__)
  // avx2 is checked at runtime, since the build targets baseline x86-64
  if (__builtin_cpu_supports("avx2")) {
    encode_biallelic_8bit_avx2(out, flags, n_samples, genotypes);
  } else {
    encode_biallelic_8bit_range(out, flags, 0, n_samples, genotypes);
  }
#elif defined(__aarch64__)
  // every aarch64 has NEON, so there is nothing to check
  encode_biallelic_8bit_neon(out, flags, n_samples, genotypes);
#else
  encode_biallelic_8bit_range(out, flags, 0, n_samples, genotypes);
#endif
  return genotype_offset + 2 * n_samples;
}

static std::uint32_t encode_unphased(std::vector<std::uint8_t> &encoded,
                     std::uint32_t genotype_offset,
                     std::uint32_t ploidy_offset,
                     std::uint32_t n_samples,
                     std::uint16_t n_alleles,
                     bool constant_ploidy,
                     std::uint32_t max_ploidy,
                     double *genotypes,
                     std::uint8_t &bit_depth)
{
  int _ploid = (int)max_ploidy;
  int _n_alleles = (int)n_alleles;
  bool phased = false;
  std::uint32_t max_probs = get_max_probs(_ploid, _n_alleles, phased);
  std::uint32_t n_probs = max_probs;  // for storing probs per person

  // only needed when the ploidy varies, since a constant ploidy reuses max_probs
  std::array<std::uint32_t, 64> ploidy_probs;
  if (!constant_ploidy) {
    ploidy_probs = probs_per_ploidy(_n_alleles, (int) max_ploidy, phased);
  }

  double factor = std::pow(2, bit_depth) - 1;
  bool missing;
  std::uint32_t bit_idx=0;
  std::uint32_t byte_idx;
  std::uint32_t bit_remainder;
  std::uint64_t window;
  double g;
  double cumulative;
  std::uint64_t running, previous, value;
  // counted alongside i rather than recovered as i / max_probs, which divides by a
  // runtime value once or twice per sample
  std::uint32_t sample_idx = 0;
  for (std::uint32_t i=0; i<(n_samples*max_probs); i+= max_probs, sample_idx++) {
    if (!constant_ploidy) {
      _ploid = (int)(encoded[ploidy_offset + sample_idx] &= 63);
      n_probs = ploidy_probs[_ploid];
    } else {
      n_probs = max_probs;
    }
    missing = missing_genotypes(&genotypes[i], n_probs);
    if (missing) {
      encoded[ploidy_offset + sample_idx] |= 0x80;
      // every probability of a missing sample encodes as zero, and the buffer is already
      // zero here, so stepping over the sample keeps the test out of the inner loop
      bit_idx += (n_probs - 1) * bit_depth;
      continue;
    }
    // the probabilities for a sample sum to 1.0, so scale their running total
    cumulative = 0.0;
    running = 0;
    for (std::uint32_t j = 0; j < (n_probs - 1); j++) {
      g = genotypes[i + j];
      check_probability(g);
      cumulative += g;
      check_cumulative(cumulative);
      previous = running;
      running = scale_cumulative(cumulative, factor, previous);
      value = running - previous;
      byte_idx = genotype_offset + (bit_idx / 8);
      if (bit_depth == 8) {
        // fast path for 8-bit genotype data
        encoded[byte_idx] = (std::uint8_t) value;
      } else {
        bit_remainder = bit_idx % 8;
        window = emplace_probability(value, &encoded[byte_idx], bit_remainder);
        std::memcpy(&encoded[byte_idx], &window, 8);
      }
      bit_idx += bit_depth;
    }
    // the reader infers the last probability rather than reading it, so check it as part
    // of the total, or a row summing over one silently changes its final probability
    g = genotypes[i + n_probs - 1];
    check_probability(g);
    cumulative += g;
    check_cumulative(cumulative);
  }
  return genotype_offset + (bit_idx / 8) + (std::uint32_t)((bit_idx % 8) > 0);
}

static std::uint32_t encode_phased(std::vector<std::uint8_t> &encoded,
                            std::uint32_t genotype_offset,
                            std::uint32_t ploidy_offset,
                            std::uint32_t n_samples,
                            std::uint16_t n_alleles,
                            bool constant_ploidy,
                            std::uint32_t max_ploidy,
                            double *genotypes,
                            std::uint8_t &bit_depth)
{
  int _ploid = (int)max_ploidy;
  int _n_alleles = (int)n_alleles;
  bool phased = true;
  std::uint32_t max_probs = get_max_probs(_ploid, _n_alleles, phased);
  std::uint32_t n_probs = max_probs; // for storing probs per person

  // as in encode_unphased, only the varying ploidy path needs the table
  std::array<std::uint32_t, 64> ploidy_probs;
  if (!constant_ploidy) {
    ploidy_probs = probs_per_ploidy(_n_alleles, (int) max_ploidy, phased);
  }

  double factor = std::pow(2, bit_depth) - 1;
  bool missing;
  std::uint32_t bit_idx = 0;
  std::uint32_t byte_idx, bit_remainder;
  std::uint64_t window;
  double g;
  double cumulative;
  std::uint64_t running, previous, value;
  std::uint32_t i = 0;
  std::uint32_t sample_idx=0;
  while (i < (n_samples * max_probs * max_ploidy)) {
    if (!constant_ploidy) {
      _ploid = (int)(encoded[ploidy_offset + sample_idx] &= 63);
      n_probs = ploidy_probs[_ploid];
    } else {
      _ploid = max_ploidy;
      n_probs = max_probs;
    }
    missing = missing_genotypes(&genotypes[i], n_probs);
    if (missing) {
      encoded[ploidy_offset + sample_idx] |= 0x80;
      // as in encode_unphased, every stored probability of a missing sample is zero,
      // and the buffer is already zero here, so step over the whole sample
      bit_idx += (std::uint32_t)_ploid * (n_probs - 1) * bit_depth;
      i += max_probs * max_ploidy;
      sample_idx += 1;
      continue;
    }
    // phased data is received in n_alleles * n_ploidy values, but is stored in
    // n_alleles * (n_ploidy - 1) values, where n_ploidy can differ per person.
    for (std::uint32_t k = 0; k < (std::uint32_t)_ploid; k++) {
      // repeat for each haplotype. The allele probabilities sum to 1.0 within a
      // haplotype, so the running total restarts for each one
      cumulative = 0.0;
      running = 0;
      for (std::uint32_t j = 0; j < (n_probs - 1); j++) {
        // repeat for each allele
        g = genotypes[i];
        check_probability(g);
        cumulative += g;
        check_cumulative(cumulative);
        previous = running;
        running = scale_cumulative(cumulative, factor, previous);
        value = running - previous;
        byte_idx = genotype_offset + (bit_idx / 8);
        bit_remainder = bit_idx % 8;
        window = emplace_probability(value, &encoded[byte_idx], bit_remainder);
        std::memcpy(&encoded[byte_idx], &window, 8);
        bit_idx += bit_depth;
        i += 1;
      }
      // as above, the final probability of the haplotype is inferred, not stored
      g = genotypes[i];
      check_probability(g);
      cumulative += g;
      check_cumulative(cumulative);
      i += 1;
    }
    i += (max_probs * max_ploidy) - (n_probs * _ploid);
    sample_idx += 1;
  }
  return genotype_offset + (bit_idx / 8) + (std::uint32_t)((bit_idx % 8) > 0);
}

static std::vector<std::uint8_t> encode_layout2(
                    std::uint32_t n_samples,
                    std::uint16_t n_alleles,
                    double *genotypes,
                    std::uint32_t geno_len,
                    uint8_t *ploidy,
                    std::uint8_t min_ploidy,
                    std::uint8_t max_ploidy,
                    bool phased,
                    std::uint8_t &bit_depth
                    ) 
{
  int _max_ploid = (int)max_ploidy;
  int _n_alleles = (int)n_alleles;
  std::uint32_t max_probs = get_max_probs(_max_ploid, _n_alleles, phased);
  if (phased) {
    max_probs *= max_ploidy;
  }
  if ((std::uint64_t) geno_len != (std::uint64_t) max_probs * n_samples) {
    throw std::invalid_argument("genotypes does not match n_samples * per_person_probs");
  }

  std::uint32_t probs_len = (n_samples * bit_depth) * (max_probs - 1);
  bool remainder = (probs_len % 8) > 0;
  probs_len = (probs_len / 8) + remainder;

  std::uint32_t encoded_size = 10 + n_samples + probs_len;
  std::vector<std::uint8_t> encoded(encoded_size + 8);  // extend slightly to help with variable bit depths
  std::uint32_t i=0;
  std::memcpy(&encoded[i], &n_samples, 4);
  i += 4;
  std::memcpy(&encoded[i], &n_alleles, 2);
  i += 2;
  encoded[i] = min_ploidy;
  i += 1;
  encoded[i] = max_ploidy;
  i += 1;

  // set the individuals ploidy values. We'll fill samples with missing data 
  // when we run through the genotypes
  bool constant_ploidy = min_ploidy == max_ploidy;
  const std::uint32_t ploidy_offset = i;
  if (constant_ploidy) {
    std::memset(&encoded[i], max_ploidy, n_samples);
    i += n_samples;
  } else {
    for (size_t j=0; j<n_samples; j++) {
      encoded[i] = ploidy[j];
      i += 1;
    }
  }

  encoded[i] = phased;
  i += 1;
  encoded[i] = bit_depth;
  i += 1;

  if (!phased) {
    if ((n_alleles == 2) && constant_ploidy && (max_ploidy == 2) &&
        (bit_depth == 8)) {
      // the shape nearly all real data has, which avoids the generic loop's per sample
      // overhead. It writes the same bytes encode_unphased would
      encoded_size = encode_biallelic_8bit(encoded, i, ploidy_offset, n_samples,
                                           genotypes);
    } else {
      encoded_size = encode_unphased(encoded, i, ploidy_offset, n_samples, n_alleles,
                                constant_ploidy, max_ploidy, genotypes, bit_depth);
    }
  } else {
    encoded_size = encode_phased(encoded, i, ploidy_offset, n_samples, n_alleles,
                                   constant_ploidy, max_ploidy, genotypes, bit_depth);
  }

  encoded.resize(encoded_size);
  return encoded;
}

// convenience function for constant ploidy
void CppBgenWriter::encode_genotype_data(std::uint16_t n_alleles,
                                         double *genotypes,
                                         std::uint32_t geno_len,
                                         std::uint8_t ploidy,
                                         bool phased,
                                         std::uint8_t bit_depth)
{
  std::uint8_t *ploidy_vector = {};
  encode_genotype_data(n_alleles, genotypes, geno_len, ploidy_vector, ploidy, ploidy, phased, bit_depth);
}

// encode (and compress) a genotype block, ready for write_genotype_data(). This never
// touches the file handle, so an error it raises leaves the bgen unchanged
void CppBgenWriter::encode_genotype_data(std::uint16_t n_alleles,
                                         double *genotypes,
                                         std::uint32_t geno_len,
                                         uint8_t *ploidy,
                                         std::uint8_t min_ploidy,
                                         std::uint8_t max_ploidy,
                                         bool phased,
                                         std::uint8_t bit_depth)
{
  if ((layout != 1) && (layout != 2)) {
    throw std::invalid_argument("layout must be 1 or 2");
  }
  if ((layout == 1) && (compression == 2)) {
    throw std::invalid_argument("you cannot use zstd compression with layout 1");
  }

  std::vector<std::uint8_t> encoded;
  if (layout == 1) {
    encoded = encode_layout1(genotypes, geno_len);
  } else {
    encoded = encode_layout2(n_samples, n_alleles, genotypes, geno_len, ploidy, 
                   min_ploidy, max_ploidy, phased, bit_depth);
  }

  std::vector<char> compressed;
  if (compression != 0) {
    compressed = compress(encoded, compression);
  }
  std::uint32_t compressed_len = compressed.size();

  // assemble the block exactly as it needs to appear on disk, so that writing
  // it later is a single copy which cannot fail partway on bad input
  pending.clear();
  std::uint32_t size;
  if (layout == 1) {
    if (compression != 0) {
      append_bytes(pending, &compressed_len, 4);
    }
  } else {
    if (compression == 0) {
      size = encoded.size();
      append_bytes(pending, &size, 4);
    } else {
      size = compressed_len + 4;
      append_bytes(pending, &size, 4);
      size = encoded.size();
      append_bytes(pending, &size, 4);
    }
  }

  if (compression == 0) {
    append_bytes(pending, encoded.data(), encoded.size());
  } else {
    append_bytes(pending, compressed.data(), compressed.size());
  }
}

// write the block prepared by encode_genotype_data()
std::uint64_t CppBgenWriter::write_genotype_data() {
  handle.write(pending.data(), pending.size());
  pending.clear();
  return current_position(handle);
}

}  // namespace bgen
