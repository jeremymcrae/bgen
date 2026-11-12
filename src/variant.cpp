
#include <algorithm>
#include <stdexcept>
#include <utility>

#include "variant.h"

namespace bgen {

/// read a fixed width value from the bgen, throwing out_of_range if the read failed,
/// which surfaces in python as the IndexError that ends iteration
template <typename T>
static void read_checked(std::istream & handle, T & value) {
  if (!read_value(handle, value)) {
    throw std::out_of_range("reached end of file");
  }
}

/// read a string which is prefixed by its length, and check the reads succeeded
template <typename LenType>
static void read_checked_string(std::istream & handle, std::string & value) {
  if (!read_prefixed_string<LenType>(handle, value)) {
    throw std::out_of_range("reached end of file");
  }
}

/// initialise a single variant with chrom, pos, rsID identifiers
///
/// This starts a Genotypes object, but this doesn't parse the genotypes until
/// required, just starts it so we can get the offset of the next variant, so as
/// to parse the bgen variants at speed.
///
///  @param handle std::istream for bgen file, shared with the CppBgenReader so
///     the file stays open for as long as this Variant might read from it
///  @param varoffset start byte for variant in bgen file
///  @param layout bgen layout version (1 or 2)
///  @param compression compression scheme (0=no compression, 1=zlib, 2=zstd)
///  @param expected_n number of samples for variant
Variant::Variant(std::shared_ptr<std::istream> handle, std::uint64_t varoffset, int layout, int compression, std::uint32_t expected_n, bool is_stdin) {
  offset = varoffset;
  if (!is_stdin) {
    handle->clear();
    handle->seekg(offset);
  }
  if (handle->eof()) {
    // check for end-of-file before reading, which is how iteration over a stdin bgen
    // terminates, since we cannot seek back to check whether another variant follows
    throw std::out_of_range("reached end of file");
  }
  std::uint32_t n_samples = expected_n;
  if (layout == 1) {
    read_checked(*handle, n_samples);
  }
  
  if (n_samples != expected_n) {
    throw std::invalid_argument("number of samples doesn't match");
  }
  
  read_checked_string<std::uint16_t>(*handle, varid);
  read_checked_string<std::uint16_t>(*handle, rsid);
  read_checked_string<std::uint16_t>(*handle, chrom);
  
  read_checked(*handle, pos);
  std::uint16_t n_alleles = 2;
  if (layout != 1) {
    read_checked(*handle, n_alleles);
  }
  
  // read each allele directly into the vector
  alleles.resize(n_alleles);
  for (auto & allele : alleles) {
    read_checked_string<std::uint32_t>(*handle, allele);
  }
  
  // 64-bit, since an uncompressed layout 1 block of n_samples * 6 can exceed 32 bits
  std::uint64_t length;
  if ((layout == 1) && (compression == 0)) {
    length = (std::uint64_t) n_samples * 6;
  } else {
    std::uint32_t stored;
    read_checked(*handle, stored);
    length = stored;
  }
  std::uint64_t geno_offset = 0;
  if (!is_stdin) {
    geno_offset = (std::uint64_t) handle->tellg();
  }
  geno.initialize(std::move(handle), layout, compression, n_alleles, n_samples, geno_offset, length, is_stdin);
  next_variant_offset = geno_offset + length;
}

std::uint32_t Variant::probs_per_sample() {
  geno.load_data_and_parse_header();
  return geno.max_probs;
}

bool Variant::phased() {
  geno.load_data_and_parse_header();
  return geno.phased;
}

/// lowest ploidy the variant gives for its samples
int Variant::min_ploidy() {
  geno.load_data_and_parse_header();
  return geno.min_ploidy;
}

/// highest ploidy the variant gives for its samples
int Variant::max_ploidy() {
  geno.load_data_and_parse_header();
  return geno.max_ploidy;
}

std::uint8_t * Variant::ploidy() {
  geno.load_data_and_parse_header();
  // the array is only filled on demand, since a constant ploidy variant does not otherwise
  // need it, and filling it costs more than decoding the genotypes
  geno.materialise_ploidy();
  return geno.ploidy.get();
}

/// get genotype probabilities for the variant as a 1-dimensional vector
///
/// This makes it easy to pass the data via cython into a numpy array, which can
/// be reshaped to a 2-D array.
void Variant::probs_1d(float * probs) {
  geno.probabilities(probs);
}

/// whether the last decode saw probabilities summing above the bit depth's maximum. The
/// decoders infer the final probability from the remainder, so an over-large sum clamps
bool Variant::probs_above_max() {
  return geno.probs_above_max;
}

/// get dosage of the alt allele (only works for biallelic variants)
void Variant::alt_dosage(float * dose) {
  geno.get_allele_dosage(dose, false);
}

/// get dosage of the minor allele (only works for biallelic variants)
void Variant::minor_allele_dosage(float * dose) {
  geno.get_allele_dosage(dose, true);
}

/// the least common of a biallelic variant's alleles
///
/// Which allele is the minor one depends on the genotypes, so this reads them
///
/// @return the minor allele
std::string Variant::get_minor_allele() {
  int idx2 = geno.get_minor_idx();
  if ((idx2 < 0) || ((std::size_t) idx2 >= alleles.size())) {
    // unreachable, since get_minor_idx only returns 0 or 1 for a biallelic variant,
    // but guards a future change from reading past the end of the alleles
    throw std::invalid_argument("bgen variant has no allele for the minor index");
  }
  return alleles[idx2];
}

std::vector<std::uint8_t> Variant::copy_data() {
  std::istream * handle = geno.stream();
  // as in Genotypes::decompress, the badbit marks a closed bgen, while other
  // error states are recoverable and just need clearing before the seek
  if (handle->bad()) {
    throw std::invalid_argument("cannot read from closed bgen file");
  }
  handle->clear();
  // 64-bit, since long alleles can take a variant past 4 GB
  std::uint64_t length = next_variant_offset - offset;
  handle->seekg(offset);
  // read in stages, so a corrupt length only allocates what the file can supply
  std::vector<std::uint8_t> data;
  std::uint64_t filled = 0;
  while (filled < length) {
    std::uint64_t chunk = std::min(length - filled, BLOCK_READ_CHUNK);
    data.resize(filled + chunk);
    if (!handle->read(reinterpret_cast<char *>(data.data() + filled), chunk)) {
      throw std::invalid_argument("could not read variant data - is the bgen truncated?");
    }
    filled += chunk;
  }
  return data;
}

} // namespace bgen
