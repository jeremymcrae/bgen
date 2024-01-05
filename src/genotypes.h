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
  Genotypes(std::ifstream* handle, int lay, int compr, int n_alleles, std::uint32_t n_samples, std::uint64_t offset, std::uint32_t length) :
     handle(handle), layout(lay), compression(compr), n_alleles(n_alleles), n_samples(n_samples), file_offset(offset), length(length) {};
  Genotypes() {};
  ~Genotypes() { clear_probs(); };
  void load_data_and_parse_header();
  float * probabilities();
  float * get_allele_dosage(bool use_alt=true, bool use_minor=false);
  bool phased=false;
  std::uint32_t max_probs=0;
  int min_ploidy=0;
  int max_ploidy=0;
  int minor_idx=0;
  std::uint8_t * ploidy={};
private:
  void decompress();
  void parse_ploidy();
  float * probabilities_layout1();
  float * probabilities_layout2();
  void fast_haplotype_probs(char * uncompressed, float * probs, std::uint32_t idx, std::uint32_t & nrows);
  void ref_dosage_fast(char * uncompressed, std::uint32_t idx, float * dose);
  void ref_dosage_slow(char * uncompressed, std::uint32_t idx, float * dose);
  void swap_allele_dosage(float * dose);
  int find_minor_allele(float * dose);
  void clear_probs();
  std::ifstream* handle;
  int layout;
  int compression;
  int n_alleles;
  std::uint32_t n_samples;
  std::uint64_t file_offset;
  std::uint32_t length;
  std::uint32_t bit_depth=0;
  std::uint32_t idx=0;
  char * uncompressed={};
  float * probs={};
  float * alt_dose={};
  float * minor_dose={};
  bool is_decompressed = false;
  bool constant_ploidy=true;
  bool has_ploidy = false;
  bool probs_parsed = false;
  bool minor_dosage_parsed = false;
  bool alt_dosage_parsed = false;
  std::vector<std::uint32_t> missing;
};

std::uint32_t get_max_probs(int &max_ploidy, int &n_alleles, bool &phased);

} // namespace bgen

#endif  // BGEN_GENOTYPES_H_
