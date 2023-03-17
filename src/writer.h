#ifndef BGEN_BGEN_H_
#define BGEN_BGEN_H_

#include <fstream>
#include <stdexcept>
#include <vector>

namespace bgen {

class BgenWriter {
  std::ofstream handle;
  std::uint32_t n_samples;
  std::uint32_t compression;
  std::uint32_t layout;
  std::uint32_t n_variants=0;
  std::uint32_t nvars_offset=8;
  std::uint32_t variant_data_offset=0;
public:
  BgenWriter(std::string &path,
             std::uint32_t n_samples,
             std::string &free_data,
             uint32_t compression,
             uint32_t layout,
             std::vector<std::string> &samples) : n_samples(n_samples),
                                                  compression(compression),
                                                  layout(layout)
  {
    handle.open(path, std::ios::out | std::ios::binary);
    write_header(free_data, samples);
    add_samples(samples);
  };
  BgenWriter() {};
  ~BgenWriter();
  void write_header(std::string &free_data,
                    std::vector<std::string> &samples);
  void add_samples(std::vector<std::string> &samples);
  void write_variant_header(std::string &varid,
                            std::string &rsid,
                            std::string &chrom,
                            std::uint32_t &pos,
                            std::vector<std::string> &alleles,
                            std::uint32_t _n_samples);
  void add_genotype_data(std::uint16_t n_alleles,
                         double *genotypes,
                         std::uint32_t geno_len,
                         std::uint8_t ploidy = 2,
                         bool phased = 0,
                         std::uint8_t bit_depth = 8);
  void add_genotype_data(std::uint16_t n_alleles,
                         double *genotypes,
                         std::uint32_t geno_len,
                         std::vector<uint8_t> &ploidy,
                         std::uint8_t min_ploidy = 2,
                         std::uint8_t max_ploidy = 2,
                         bool phased = 0,
                         std::uint8_t bit_depth = 8);
};

} // namespace bgen

#endif  // BGEN_BGEN_H_
