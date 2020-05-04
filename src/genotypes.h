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
  float * probs;
  std::vector<int> missing;
public:
  Genotypes(std::ifstream* handle, int lay, int compr, int n_alleles, int n_samples) :
    handle(handle), layout(lay), compression(compr), n_alleles(n_alleles), n_samples(n_samples) {
      std::uint32_t length;
      handle->read(reinterpret_cast<char*>(&length), sizeof(length));
      offset = handle->tellg();
      next_var_offset = offset + length;
  };
  Genotypes() {};
  ~Genotypes() { clear_probs(); };
  void decompress(char * bytes, int compressed_len, char * decompressed, int decompressed_len);
  void parse_ploidy(char * uncompressed, int & idx);
  float * parse_layout1(char *);
  float * parse_layout2(char *);
  float * probabilities();
  void clear_probs();
  bool phased;
  int max_probs = 0;
  bool constant_ploidy;
  int min_ploidy;
  int max_ploidy;
  std::uint8_t * ploidy;
  std::uint64_t next_var_offset;
};

} // namespace bgen

#endif  // BGEN_GENOTYPES_H_
