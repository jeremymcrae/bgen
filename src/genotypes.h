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
public:
  Genotypes() {}
  ~Genotypes() { clear_probs(); }
  void initialize(std::istream* _handle, 
           int lay,
           int compr,
           int _n_alleles, 
           std::uint32_t _n_samples,
           std::uint64_t _offset,
           std::uint32_t _length,
           bool _is_stdin=false) {
      handle = _handle;
      layout = lay;
      compression = compr;
      n_alleles = _n_alleles;
      n_samples = _n_samples;
      file_offset = _offset;
      length = _length;
      is_stdin = _is_stdin;
      if (is_stdin) {
        load_data_and_parse_header();
      }
    }
  void load_data_and_parse_header();
  void probabilities(float * probs);
  void get_allele_dosage(float * dose, bool use_alt=true, bool use_minor=false);
  bool phased=false;
  std::uint32_t max_probs=0;
  int min_ploidy=0;
  int max_ploidy=0;
  int minor_idx=0;
  std::uint8_t * ploidy={};
private:
  void decompress();
  void parse_ploidy();
  void probabilities_layout1(char * uncompressed, std::uint32_t idx, float * probs, std::uint32_t & nrows);
  void probabilities_layout2(char * uncompressed, std::uint32_t idx, float * probs, std::uint32_t & nrows);
  void fast_haplotype_probs(char * uncompressed, std::uint32_t idx, float * probs, std::uint32_t & nrows);
  void ref_dosage_fast(char * uncompressed, std::uint32_t idx, float * dose, std::uint32_t nrows);
  void ref_dosage_slow_unphased(char * uncompressed, std::uint32_t idx, float * dose, std::uint32_t nrows);
  void ref_dosage_slow_phased(char * uncompressed, std::uint32_t idx, float * dose, std::uint32_t nrows);
  void swap_allele_dosage_simple(float * dose);
  void swap_allele_dosage_complex(float * dose);
  int find_minor_allele(float * dose);
  void clear_probs();
  std::istream* handle;
  int layout;
  int compression;
  int n_alleles;
  std::uint32_t n_samples;
  std::uint64_t file_offset;
  std::uint32_t length;
  bool is_stdin;
  std::uint32_t bit_depth=0;
  std::uint32_t idx=0;
  char * uncompressed={};
  bool is_decompressed = false;
  bool constant_ploidy=true;
  bool has_ploidy = false;
  std::vector<std::uint32_t> missing;
};

std::uint32_t get_max_probs(int &max_ploidy, int &n_alleles, bool &phased);

} // namespace bgen

#endif  // BGEN_GENOTYPES_H_
