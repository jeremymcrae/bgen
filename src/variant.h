#ifndef BGEN_VARIANT_H_
#define BGEN_VARIANT_H_

#include <cstdint>
#include <istream>
#include <memory>
#include <string>
#include <vector>

#include "genotypes.h"
#include "utils.h"

namespace bgen {

class Variant {
  Genotypes geno;
public:
  Variant(std::shared_ptr<std::istream> handle, std::uint64_t varoffset, int layout, int compression, std::uint32_t expected_n, bool is_stdin=false);
  Variant() {}
  std::uint32_t probs_per_sample();
  void alt_dosage(float * dosage);
  void minor_allele_dosage(float * dosage);
  std::string get_minor_allele();
  void probs_1d(float * probs);
  bool phased();
  int min_ploidy();
  int max_ploidy();
  /// check for probabilities summing above bit-depth max, indicates malformed bgen
  bool probs_above_max();
  std::uint8_t * ploidy();
  std::vector<std::uint8_t> copy_data();
  
  std::uint64_t offset = 0;
  std::string varid;
  std::string rsid;
  std::string chrom;
  std::uint32_t pos = 0;
  std::vector<std::string> alleles;
  std::uint64_t next_variant_offset = 0;
};

} // namespace bgen

#endif  // BGEN_VARIANT_H_
