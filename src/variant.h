#ifndef BGEN_VARIANT_H_
#define BGEN_VARIANT_H_

#include <cstdint>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <iterator>

#include "genotypes.h"

namespace bgen {

class Variant {
  std::uint64_t offset;
  Genotypes geno;
public:
  Variant(std::ifstream & handle, int layout, int compression);
  std::uint64_t next_variant_offset();
  std::string alt();
  std::vector<std::vector<double>> genotypes();
  std::vector<double> alt_dosage();
  
  std::uint32_t n_samples = 0;
  std::string varid;
  std::string rsid;
  std::string chrom;
  std::uint32_t pos;
  std::uint16_t n_alleles;
  std::vector<std::string> alleles;
};

} // namespace bgen

#endif  // BGEN_VARIANT_H_
