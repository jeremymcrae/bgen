#ifndef BGEN_GENOTYPES_H_
#define BGEN_GENOTYPES_H_

#include <cstdint>
#include <fstream>
#include <memory>
#include <vector>
#include <array>
#include <string>
#include <sstream>

namespace bgen {

/// bytes of padding to allocate past the end of the decompressed genotype data
///
/// probabilities_layout2 reads 8 bytes at a time so that a probability can be
/// extracted in a single read, even at bit depths which straddle byte
/// boundaries. Without this padding, the read for the final probability would
/// run past the end of the buffer.
const std::uint32_t PROBS_READ_PAD = 8;

/// genotype data for one variant, which owns its decompressed buffers
///
/// The buffers are held in unique_ptrs, so a Genotypes can be moved but not
/// copied. That matters because a Variant holds one of these by value, and
/// parse_all_variants moves Variants into a vector. Owning the buffers through
/// raw pointers instead would let a copy duplicate them, and both copies would
/// then free the same memory.
class Genotypes {
public:
  Genotypes() {}
  void initialize(std::shared_ptr<std::istream> _handle,
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
  int get_minor_idx();
  bool phased=false;
  std::uint32_t max_probs=0;
  int min_ploidy=0;
  int max_ploidy=0;
  int minor_idx=0;
  // whether a sample's stored probabilities were found to sum above the maximum the bit
  // depth can hold, which means the bgen is malformed. The final probability of a group is
  // not stored but inferred as the remainder, so an over-large sum makes it negative. The
  // negative values are clamped to zero.
  bool probs_above_max = false;
  std::unique_ptr<std::uint8_t[]> ploidy;
  // fills the ploidy array when every sample shares a ploidy, in which case parse_ploidy
  // leaves it empty. Public because variant.cpp hands the array out to callers
  void materialise_ploidy();
private:
  void decompress();
  void parse_ploidy();
  std::uint64_t probability_bytes();
  void check_block_size();
  void probabilities_layout1(char * uncompressed, std::uint32_t idx, float * probs, std::uint32_t & nrows);
  void probabilities_layout2(char * uncompressed, std::uint32_t idx, float * probs, std::uint32_t & nrows);
  void fast_haplotype_probs(char * uncompressed, std::uint32_t idx, float * probs, std::uint32_t & nrows);
  void ref_dosage_fast(char * uncompressed, std::uint32_t idx, float * dose, std::uint32_t nrows);
  void ref_dosage_slow_unphased(char * uncompressed, std::uint32_t idx, float * dose, std::uint32_t nrows);
  void ref_dosage_slow_phased(char * uncompressed, std::uint32_t idx, float * dose, std::uint32_t nrows);
  void swap_allele_dosage_simple(float * dose);
  void swap_allele_dosage_complex(float * dose);
  int find_minor_allele(float * dose);
  std::shared_ptr<std::istream> handle;
  int layout = 0;
  int compression = 0;
  int n_alleles = 0;
  std::uint32_t n_samples = 0;
  std::uint64_t file_offset = 0;
  std::uint32_t length = 0;
  bool is_stdin = false;
  std::uint32_t bit_depth=0;
  std::uint32_t idx=0;
  std::unique_ptr<char[]> uncompressed;
  // size of the decompressed genotype block, so that the reads below can be
  // bounded by the data which is actually present. The block length comes from
  // the bgen itself, so it cannot be assumed to match what the other header
  // fields say the variant needs.
  std::uint32_t uncompressed_len = 0;
  // whether the block has been checked as big enough for the variant. The check
  // has to be remembered, since the header is only parsed once, so a caller
  // which catches the error and asks for the probabilities again would otherwise
  // skip straight past it to the reads which the check exists to prevent.
  bool block_checked = false;
  bool is_decompressed = false;
  bool constant_ploidy=true;
  bool has_ploidy = false;
  // whether minor_idx has been worked out for this variant. Finding it needs the
  // dosages of the whole cohort, so the answer is kept once known, both to save
  // repeating that and because the layout 1 dosage path appends to the missing
  // list as it goes.
  bool minor_known = false;
  std::vector<std::uint32_t> missing;
};

std::uint32_t get_max_probs(int &max_ploidy, int &n_alleles, bool &phased);

} // namespace bgen

#endif  // BGEN_GENOTYPES_H_
