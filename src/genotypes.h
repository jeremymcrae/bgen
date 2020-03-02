#ifndef BGEN_GENOTYPES_H_
#define BGEN_GENOTYPES_H_

#include <cstdint>
#include <fstream>
#include <vector>
#include <array>
#include <string>
#include <sstream>

namespace bgen {

class Genotypes {
  std::uint64_t offset;
  std::ifstream* handle;
  int layout;
  int compression;
  int n_alleles;
  int n_samples;
public:
  Genotypes(std::ifstream* handle, int lay, int compr, int n_alleles, int n_samples) :
    handle(handle), layout(lay), compression(compr), n_alleles(n_alleles), n_samples(n_samples) {
      std::uint32_t length;
      handle->read(reinterpret_cast<char*>(&length), sizeof(length));
      offset = handle->tellg();
      next_var_offset = offset + length;
  };
  Genotypes() {};
  // std::vector<char> decompress(std::vector<char> bytes, int length);
  std::vector<char> decompress(char * bytes, int compressed_len, int decompressed_len);
  void parse_layout1(std::vector<char>);
  void parse_layout2(std::vector<char>);
  std::vector<std::vector<double>> genotypes();
  std::vector<int> ploidy;
  std::vector<std::vector<double>> parsed;
  std::uint64_t next_var_offset;
};

} // namespace bgen

#endif  // BGEN_GENOTYPES_H_
