
#include <array>
#include <cstring>
#include <cmath>
#include <fstream>
#include <stdexcept>
#include <vector>

#include <iostream>

#include "zstd/lib/zstd.h"
#include <zlib.h>

#include "writer.h"
#include "genotypes.h"

namespace bgen {

// append raw bytes to a buffer destined for the file handle
static void append_bytes(std::vector<char> &buf, const void *data, std::size_t size) {
  const char * bytes = reinterpret_cast<const char *>(data);
  buf.insert(buf.end(), bytes, bytes + size);
}

/// @brief the current write position, or a clear error if there is not one
///
/// The positions tellp() reports are written into the bgen header as the variant data
/// offset, and handed back to python to build the .bgi index, so a wrong one silently
/// produces a file whose index points at the wrong bytes.
///
/// tellp() has two ways to not give an answer, and only one of them is noisy. On a
/// stream that has already failed it returns -1 immediately, without touching the stream
/// state, so the exception mask has nothing to report and (std::uint64_t) -1 becomes a
/// recorded offset of 18446744073709551615. On a healthy output that cannot seek, such as
/// a pipe, the seek itself fails and the mask does throw, but only with "basic_ios::clear:
/// iostream error", which does not say what went wrong.
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

// finish the file off by writing the number of variants, and where the variant
// data starts
//
// This can throw (e.g. if the disk is full), so call it explicitly to have
// write errors reported, rather than relying on the destructor.
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
/// varids, rsids, chromosomes and sample IDs are all written as a two byte
/// length followed by the characters, so anything longer than 65535 has no valid
/// encoding. Without this the length silently wraps while the full string is
/// still written, which leaves a file whose block lengths are self-consistent but
/// whose fields are all read from the wrong offsets - so the bgen looks fine
/// until something tries to read it back.
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
  // The bgen stores this as a four byte field, so an offset past that has no valid
  // encoding. Casting a larger one down wraps it, which leaves a file whose header
  // sends readers to the wrong place for the variants, so check it in 64 bits first.
  // The block length check above bounds the sample block alone, but the offset also
  // carries the header, so the two together can still overflow.
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
  // write() rather than a std::ostreambuf_iterator, because the iterator reports a
  // failed write only in its own state. It leaves the stream looking good, so the
  // exception mask that covers every other write here would not see the bytes go
  // missing, and this function would still return a plausible offset for the index.
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

/// tolerance for genotype probabilities which sit just outside the legal range.
///
/// Probabilities frequently arrive after a round trip through float32, since
/// that is what the reader hands back, and that can shift a value by a few
/// multiples of the float32 epsilon of 1.2e-07. The error also accumulates when
/// summing the probabilities for a sample. Values within this much of the limit
/// are treated as legal (and are clamped when encoded), anything further out is
/// reported as a mistake.
const double PROB_TOLERANCE = 1e-6;

/// @brief report an unusable genotype probability
///
/// Kept out of line so that the checks below stay small enough to inline into
/// the encoding loops, since they run once per probability written.
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
/// The comparison is written as a single negated range test so that nan fails it
/// too, rather than slipping through as it would with a plain
/// 'g < -PROB_TOLERANCE' test. A nan here means only part of a sample is
/// missing, since fully missing samples are spotted before this and encoded as
/// zeroes.
static inline void check_probability(double g) {
  if (!((g >= -PROB_TOLERANCE) && (g <= 1.0 + PROB_TOLERANCE))) {
    raise_probability_error(g);
  }
}

/// @brief check the probabilities stored for a sample do not sum above one
///
/// Layout 2 leaves the last probability of each group out of the file, and the
/// reader infers it as one minus the sum of the stored values, so a group
/// summing above one has no valid encoding. Note this only bounds the stored
/// values from above: probabilities which sum to less than one are still
/// written, since the caller may legitimately be storing a partial
/// distribution.
static inline void check_cumulative(double cumulative) {
  if (cumulative > 1.0 + PROB_TOLERANCE) {
    raise_cumulative_error(cumulative);
  }
}

/// @brief round a non-negative double to the nearest integer, halves away from zero
///
/// std::round is an out of line libm call on the platforms this builds for, and it sits
/// in the innermost encode loop, once per stored probability. Truncation plus a
/// comparison inlines to a few instructions instead. Both the truncation and the
/// subtraction below are exact, so this avoids the usual std::floor(x + 0.5) trap where
/// x + 0.5 rounds up in the addition and pushes the result to the wrong integer.
///
/// Only defined for finite non-negative x below 2^53, which is what the callers pass.
static inline std::uint64_t round_nonneg(double x) {
  std::uint64_t truncated = (std::uint64_t) x;
  return truncated + (x - (double) truncated >= 0.5);
}

/// @brief decide whether a sample's genotypes are all missing, rejecting a partial row
///
/// Marked inline because the hot encode paths call it once per sample, for as few as three
/// values, where the call itself costs more than the work. The nan count is deliberately
/// kept branchless: an early exit on the first non-nan measured faster when few samples
/// are missing but far slower around half, where the branch stops being predictable.
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
      // layout 1 scales by 32768 and stores 16-bit values, so the field itself
      // would hold values up to 2.0. Hold it to the same range as layout 2, so
      // that the two layouts accept the same genotypes
      check_probability(g);
      cumulative += g;
      check_cumulative(cumulative);
      // check_probability allows a small negative, down to -PROB_TOLERANCE, which
      // round_nonneg is not defined for. Such a value scales to zero anyway, which is
      // what std::round followed by the cast to an integer used to produce
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
/// Layout 2 stores all but the last probability for a sample (or for a
/// haplotype, if phased), and the reader infers the last one as
/// max - sum(stored), so the stored values have to sum to at most max. Scaling
/// each probability on its own breaks that: at a bit depth of 8, two
/// probabilities of 0.5 each round to 128, which sums to 256 when the maximum
/// is 255, and the reader then infers a negative final probability. Scaling the
/// running total instead means the stored values sum to round(total * max), so
/// the constraint holds for any probabilities that sum to 1.0.
///
/// The bounds are applied to the double, before the cast to an integer, since a cast
/// is only defined for values the integer type can represent. That also keeps nan and
/// inf away from it, though check_probability and check_cumulative mean neither can
/// arrive here in the first place.
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

  // read the full 64 bits that the caller writes back, so that the bits already
  // in the buffer survive. At a bit depth of 32 the shifted probability can span
  // 39 bits, so a narrower read would drop the top of it, and would blank
  // whatever the window covers beyond the bytes it read.
  std::memcpy(&window, encoded, sizeof(window));
  window |= (value << bit_remainder);
  return window;
}

/// @brief precompute number of probabilities per sample for every ploidy
///
/// Only entries up to max_ploidy are filled. Higher ploidies cannot occur, and for a
/// variant with many alleles they are not even representable (12 alleles at a ploidy
/// of 63 exceeds a 32-bit count), so computing them would refuse a valid variant.
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
  // counted alongside i rather than recovered as i / max_probs, which is a division by
  // a value only known at runtime, once or twice for every sample. encode_phased
  // already tracks the index this way
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
      // Every probability of a missing sample encodes as zero: the loop below would
      // substitute zero for each value, leaving the running total at zero, so each
      // stored value is zero too. The buffer is zero initialised and nothing has
      // written this sample's bits yet, so stepping over them writes the same bytes
      // the loop would have. That also keeps the missing test out of the inner loop
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
    // the last probability is inferred by the reader as one minus the stored
    // values, rather than being stored itself, so check it as part of the total.
    // Otherwise a sample whose stored values are individually fine but whose
    // total runs over one would be written with a different final probability
    // than the caller passed in.
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
    encoded_size = encode_unphased(encoded, i, ploidy_offset, n_samples, n_alleles,
                                constant_ploidy, max_ploidy, genotypes, bit_depth);
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

// encode (and compress) a genotype block, ready for write_genotype_data()
//
// This never touches the file handle, so any error it raises (e.g. mismatched
// genotype lengths, or inconsistent missingness) leaves the bgen unchanged.
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
