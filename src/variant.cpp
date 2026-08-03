
#include <stdexcept>
#include <cmath>
#include <iostream>

#include "variant.h"

namespace bgen {

/// read a fixed width value from the bgen, throwing if the read failed
///
/// Throws out_of_range, which surfaces in python as IndexError, so that reaching
/// the end of the file ends iteration.
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
///  @param _handle std::istream for bgen file, shared with the CppBgenReader so
///     the file stays open for as long as this Variant might read from it
///  @param varoffset start byte for variant in bgen file
///  @param layout bgen layout version (1 or 2)
///  @param compression compression scheme (0=no compression, 1=zlib, 2=zstd)
///  @param expected_n number of samples for variant
Variant::Variant(std::shared_ptr<std::istream> _handle, std::uint64_t & varoffset, int layout, int compression, int expected_n, bool is_stdin) : handle(_handle) {
  offset = varoffset;
  if (!is_stdin) {
    handle->clear();
    handle->seekg(offset);
  }
  if (handle->eof()) {
    // check for end-of-file before reading, so we don't try to read after EOF.
    // This is how iteration over a stdin bgen terminates, since we cannot seek
    // back on stdin to check whether another variant follows.
    throw std::out_of_range("reached end of file");
  }
  if (layout == 1) {
    read_checked(*handle, n_samples);
  } else {
    n_samples = expected_n;
  }
  
  if ((int) n_samples != expected_n) {
    throw std::invalid_argument("number of samples doesn't match");
  }
  
  read_checked_string<std::uint16_t>(*handle, varid);
  read_checked_string<std::uint16_t>(*handle, rsid);
  read_checked_string<std::uint16_t>(*handle, chrom);
  
  read_checked(*handle, pos);
  if (layout == 1) {
    n_alleles = 2;
  } else {
    read_checked(*handle, n_alleles);
  }
  
  alleles.reserve(n_alleles);
  for (int x=0; x < n_alleles; x++) {
    std::string allele;
    read_checked_string<std::uint32_t>(*handle, allele);
    alleles.push_back(allele);
  }
  
  std::uint32_t length;
  if ((layout == 1) && (compression == 0)) {
    length = n_samples * 6;
  } else {
    read_checked(*handle, length);
  }
  std::uint64_t geno_offset = 0;
  if (!is_stdin) {
    geno_offset = (std::uint64_t) handle->tellg();
  }
  geno.initialize(handle, layout, compression, n_alleles, n_samples, geno_offset, length, is_stdin);
  next_variant_offset = geno_offset + length;
}

int Variant::probs_per_sample() {
  if (geno.max_probs == 0) {
    geno.load_data_and_parse_header();
  }
  return geno.max_probs;
}

bool Variant::phased() {
  if (geno.max_probs == 0) {
    geno.load_data_and_parse_header();
  }
  return geno.phased;
}

std::uint8_t * Variant::ploidy() {
  if (geno.max_probs == 0) {
    geno.load_data_and_parse_header();
  }
  return geno.ploidy.get();
}

/// get genotype probabilities for the variant as a 1-dimensional vector
///
/// This makes it easy to pass the data via cython into a numpy array, which can
/// be reshaped to a 2-D array.
void Variant::probs_1d(float * probs) {
  geno.probabilities(probs);
}

/// get dosage of the alt allele (only works for biallelic variants)
void Variant::alt_dosage(float * dose) {
  geno.get_allele_dosage(dose, true, false);
  minor_allele = alleles[geno.minor_idx];
}

/// get dosage of the minor allele (only works for biallelic variants)
void Variant::minor_allele_dosage(float * dose) {
  geno.get_allele_dosage(dose, false, true);
  minor_allele = alleles[geno.minor_idx];
}

/// the least common of a biallelic variant's alleles
///
/// Which allele is the minor one depends on the genotypes, so this reads them,
/// rather than reporting whatever a previous dosage call happened to leave behind.
///
/// @return the minor allele
std::string Variant::get_minor_allele() {
  int idx2 = geno.get_minor_idx();
  if ((idx2 < 0) || ((std::size_t) idx2 >= alleles.size())) {
    // get_minor_idx only returns 0 or 1, and it throws for anything other than a
    // biallelic variant, so this is unreachable. It is here so that a future
    // change to either cannot turn into a read past the end of the alleles.
    throw std::invalid_argument("bgen variant has no allele for the minor index");
  }
  minor_allele = alleles[idx2];
  return minor_allele;
}

std::vector<std::uint8_t> Variant::copy_data() {
  // as in Genotypes::decompress, the badbit marks a closed bgen, while other
  // error states are recoverable and just need clearing before the seek
  if (handle->bad()) {
    throw std::invalid_argument("cannot read from closed bgen file");
  }
  handle->clear();
  std::uint32_t length = next_variant_offset - offset;
  std::vector<std::uint8_t> data(length);
  handle->seekg(offset);
  if (!handle->read(reinterpret_cast<char *>(data.data()), length)) {
    throw std::invalid_argument("could not read variant data - is the bgen truncated?");
  }
  return data;
}

} // namespace bgen
