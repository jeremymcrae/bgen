#ifndef BGEN_WRITER_H_
#define BGEN_WRITER_H_

#include <cstddef>
#include <cstdint>
#include <fstream>
#include <string>
#include <vector>

namespace bgen {

class CppBgenWriter {
  std::ofstream handle;
  std::uint32_t n_samples;
  std::uint32_t compression;
  std::uint32_t layout;
  std::uint32_t n_variants=0;
  bool closed=false;
  // genotype block (length prefixes included) held between encoding and writing
  std::vector<char> pending;
  void write_header(const std::string &free_data,
                    const std::vector<std::string> &samples);
  void add_samples(const std::vector<std::string> &samples);
public:
  CppBgenWriter(const std::string &path,
             std::uint32_t _n_samples,
             const std::string &free_data,
             std::uint32_t _compression,
             std::uint32_t _layout,
             const std::vector<std::string> &samples) : n_samples(_n_samples),
                                                        compression(_compression),
                                                        layout(_layout)
  {
    handle.open(path, std::ios::out | std::ios::binary);
    if (!handle.is_open()) {
      throw std::ios_base::failure("cannot open '" + path + "' for writing");
    }
    handle.exceptions(std::ios::badbit | std::ios::failbit);
    write_header(free_data, samples);
    add_samples(samples);
  }
  ~CppBgenWriter();
  std::uint64_t write_variant_header(const std::string &varid,
                                     const std::string &rsid,
                                     const std::string &chrom,
                                     std::uint32_t pos,
                                     const std::vector<std::string> &alleles);
  std::uint64_t write_variant_direct(const char * data, std::size_t size);
  void encode_genotype_data(std::uint16_t n_alleles,
                            double *genotypes,
                            std::uint32_t geno_len,
                            const std::uint8_t *ploidy,
                            std::uint8_t min_ploidy,
                            std::uint8_t max_ploidy,
                            bool phased,
                            std::uint8_t bit_depth);
  std::uint64_t write_genotype_data();
  void close();
};

} // namespace bgen

#endif  // BGEN_WRITER_H_
