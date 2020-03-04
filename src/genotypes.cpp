
#include <algorithm>
#include <iterator>
#include <stdexcept>
#include <bitset>
#include <cmath>
#include <cassert>

// #include "zstd.h"
#include <zlib.h>

#include "genotypes.h"
#include "utils.h"

namespace bgen {

std::vector<char> zlib_uncompress(char * input, int compressed_len, int decompressed_len) {
  /* uncompress a char array with zlib
  */
  char decompressed[decompressed_len];
  
  z_stream infstream;
  infstream.zalloc = Z_NULL;
  infstream.zfree = Z_NULL;
  infstream.opaque = Z_NULL;
  
  infstream.avail_in = compressed_len; // size of input
  infstream.next_in = (Bytef *) input; // input char array
  infstream.avail_out = (uInt) sizeof(decompressed); // size of output
  infstream.next_out = (Bytef *) decompressed; // output char array
  
  inflateInit(&infstream);
  inflate(&infstream, Z_NO_FLUSH);
  inflateEnd(&infstream);
  
  return std::vector<char> {decompressed, decompressed+infstream.total_out};
}

std::vector<char> Genotypes::decompress(char * bytes, int compressed_len, int decompressed_len) {
  /* decompress the probabilty data
  */
  std::vector<char> decompressed;
  switch (compression) {
    case 0: { // no compression
      decompressed = std::vector<char> {bytes, bytes + compressed_len};
      break;
    }
    case 1: { // zlib
      decompressed = zlib_uncompress(bytes, compressed_len, decompressed_len);
      break;
    }
    case 2: { //zstd
      throw std::invalid_argument("zstd decompression not implemented yet");
      break;
    }
  }
  return decompressed;
}

void Genotypes::parse_layout1(std::vector<char> probs) {
  /* parse probabilities for layout1
  */
  ploidy = {};
  parsed = {};
  int end;
  std::vector<float> sample;
  float prob;
  for (int start=0; start < n_samples; start++) {
    end = start + 6;
    for (int x=start; x<end; x+2) {
      std::uint16_t value = *reinterpret_cast<const std::uint16_t*>(&probs[x]);
      prob = (float) value / 32768;
      sample.push_back(prob);
    }
    if ((sample[0] == 0.0) & (sample[1] == 0.0) & (sample[2] == 0.0)) {
      sample = {std::nan("1"), std::nan("1"), std::nan("1")};
    }
    
    parsed.push_back(sample);
    ploidy.push_back(2);
  }
}

void Genotypes::parse_layout2(std::vector<char> uncompressed) {
  /* parse probabilities for layout2
  */
  int idx = 0;
  std::uint32_t nn_samples = *reinterpret_cast<const std::uint32_t*>(&uncompressed[idx]);
  idx += sizeof(std::uint32_t);
  std::uint16_t allele_check = *reinterpret_cast<const std::uint16_t*>(&uncompressed[idx]);
  idx += sizeof(std::uint16_t);
  if (allele_check != n_alleles) {
    throw std::invalid_argument("number of alleles doesn't match!");
  }
  
  int min_ploidy = (int) *reinterpret_cast<const std::uint8_t*>(&uncompressed[idx]);
  idx += sizeof(std::uint8_t);
  int max_ploidy = (int) *reinterpret_cast<const std::uint8_t*>(&uncompressed[idx]);
  idx += sizeof(std::uint8_t);
  
  // get ploidy and missing states. this uses 15 milliseconds for 500k samples
  ploidy = std::vector<std::uint8_t>(n_samples);
  std::vector<bool> missing(n_samples);
  std::uint8_t flags;
  std::uint8_t mask = 63;
  for (int x=0; x < n_samples; x++) {
    flags = uncompressed[idx];
    ploidy[x] = mask & flags;
    missing[x] = flags >> 7;
    idx += 1;
  }
  
  int phased = (int) *reinterpret_cast<const std::uint8_t*>(&uncompressed[idx]);
  idx += sizeof(std::uint8_t);
  int bit_depth = (int) *reinterpret_cast<const std::uint8_t*>(&uncompressed[idx]);
  if (bit_depth < 1 | bit_depth > 32) {
    throw std::invalid_argument("probabilities bit depth out of bounds");
  }
  
  if (!(bit_depth == 8 | bit_depth == 16 | bit_depth == 32)) {
    throw std::invalid_argument("probabilities bit depth isn't a standard type (8, 16 or 32)");
  }
  
  idx += sizeof(std::uint8_t);
  float divisor = (float) (std::pow(2, (int) bit_depth)) - 1;
  
  // figure out the maximum number of probabilities across the individuals. This
  // allows preallocating the probs vector rather than clearing at each sample.
  int max_sample_len;
  if (phased) {
    max_sample_len = max_ploidy * (n_alleles - 1);
  } else {
    max_sample_len = n_choose_k(max_ploidy + n_alleles - 1, n_alleles - 1) - 1;
  }
  
  // get genotype/allele probabilities
  int bit_len = (int) bit_depth / 8;
  std::uint32_t n_probs;
  float prob;
  float remainder;
  std::vector<float> sample(max_sample_len + 1);
  parsed = std::vector<std::vector<float>>(n_samples, sample);
  for (int start=0; start < n_samples; start++) {
    // calculate the number of probabilities per sample (depends on whether the
    // data is phased, the sample ploidy and the number of alleles)
    if (phased) {
      n_probs = ploidy[start] * (n_alleles - 1);
    } else {
      n_probs = n_choose_k(ploidy[start] + n_alleles - 1, n_alleles - 1) - 1;
    }
    remainder = 1.0;
    for (int x=0; x<n_probs; x++) {
      if (bit_depth == 8) {
        prob = *reinterpret_cast<const std::uint8_t*>(&uncompressed[idx]) / divisor;
      } else if (bit_depth == 16) {
        prob = *reinterpret_cast<const std::uint16_t*>(&uncompressed[idx]) / divisor;
      } else if (bit_depth == 32) {
        prob = *reinterpret_cast<const std::uint32_t*>(&uncompressed[idx]) / divisor;
      }
      idx += bit_len;
      remainder -= prob;
      parsed[start][x] = prob;
    }
    parsed[start][n_probs] = remainder;
    
    // if the sample has missing data, just set values to NA
    if (missing[start]) {
      for (int x=0; x<(n_probs + 1); x++) {
        parsed[start][x] = std::nan("1");
      }
    }
  }
}

std::vector<std::vector<float>> & Genotypes::probabilities() {
  /* parse genotype data for a single variant
  */
  handle->seekg(offset);  // about 1 microsecond
  
  bool decompressed_field = false;
  std:uint32_t decompressed_len;
  if (compression != 0) {
    if (layout == 1) {
      decompressed_len = n_samples * 6;
    } else if (layout == 2) {
      decompressed_field = true;
      handle->read(reinterpret_cast<char*>(&decompressed_len), sizeof(std::uint32_t));
    }
  }
  
  std::uint32_t compressed_len = next_var_offset - offset - decompressed_field * 4;
  char geno[compressed_len];
  handle->read(&geno[0], compressed_len); // about 70 microseconds
  auto uncompressed = decompress(geno, (int) compressed_len, (int) decompressed_len);  // about 2-3 milliseconds
  
  switch (layout) {
    case 1: {
      parse_layout1(uncompressed);
      break;
    }
    case 2: {
      parse_layout2(uncompressed);
      break;
    }
  }
  return parsed;
}

} //namespace bgen
