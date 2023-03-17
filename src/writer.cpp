
#include <cstring>
#include <cmath>
#include <fstream>
#include <stdexcept>
#include <vector>

#include <bitset>
#include <iostream>

#include "zstd/lib/zstd.h"
#include <zlib.h>

#include "writer.h"
#include "genotypes.h"

namespace bgen {

// write a 32-bit value at a given file offset
void write_at_offset(std::ofstream &handle, std::uint32_t &val, std::uint32_t offset=0) {
  std::uint64_t orig_pos = handle.tellp();
  handle.seekp(offset);
  handle.write(reinterpret_cast<char *>(&val), 4);
  handle.seekp(orig_pos);
}

void write_variants_offset(std::ofstream &handle, std::uint32_t &offset) {
  write_at_offset(handle, offset, 0);
}

void write_nvariants(std::ofstream &handle, std::uint32_t &offset, std::uint32_t &n_variants) {
  write_at_offset(handle, n_variants, offset);
}

// when the object is removed, finally write the number of variants, and where 
// the variant data starts
BgenWriter::~BgenWriter() {
  write_variants_offset(handle, variant_data_offset);
  write_nvariants(handle, nvars_offset, n_variants);
  handle.close();
}

void BgenWriter::write_header(std::string &free_data,
                              std::vector<std::string> &samples) {
  std::uint32_t header_len = 20 + free_data.size();
  variant_data_offset = header_len;
  write_variants_offset(handle, variant_data_offset);
  handle.seekp(4);

  // figure out the length of the header block
  handle.write(reinterpret_cast<char *>(&header_len), 4);

  // write zero variants for now, is fixed while closing
  handle.write(reinterpret_cast<char *>(&n_variants), 4);
  handle.write(reinterpret_cast<char *>(&n_samples), 4);
  handle << "bgen";
  handle << free_data;

  // check and write flags
  if (compression > 2) {
    throw std::invalid_argument("compression flag must be 0, 1, or 2");
  }
  if ((layout < 1) || layout > 2) {
    throw std::invalid_argument("layout flag must be 1, or 2");
  }
  std::uint32_t sample_id_flag = samples.size() > 0;
  std::uint32_t flags = 0;
  flags |= compression;
  flags |= layout << 2;
  flags |= sample_id_flag << 31;
  handle.write(reinterpret_cast<char *>(&flags), 4);
  }

void BgenWriter::add_samples(std::vector<std::string> &samples) {
  if (samples.size() == 0) { return; }

  if (samples.size() != n_samples) {
    throw std::invalid_argument("samples vector length doesn't match the sample count in file");
  }

  // count the number of chanaracters across all sample IDs
  std::uint32_t nchars = 0;
  for (auto &x : samples) { nchars += x.size(); }

  // write the length of the sample ID block, and number of sample IDs
  std::uint32_t samples_len = 8 + 2 * samples.size() + nchars;
  handle.write(reinterpret_cast<char *>(&samples_len), 4);
  std::uint32_t size = samples.size();
  handle.write(reinterpret_cast<char *>(&size), 4);

  // write each sample ID to the file, preceeded by each ID length
  std::uint16_t id_size;
  for (auto &x : samples) {
    id_size = x.size();
    handle.write(reinterpret_cast<char *>(&id_size), 2);
    handle << x;
  }
  variant_data_offset = handle.tellp() - 4;
  write_variants_offset(handle, variant_data_offset);
}

void BgenWriter::write_variant_header(std::string &varid,
                             std::string &rsid,
                             std::string &chrom,
                             std::uint32_t &pos,
                             std::vector<std::string> &alleles,
                             std::uint32_t _n_samples) {
  n_variants += 1;
  if (_n_samples != n_samples) {
    throw std::invalid_argument("number of samples doesn't match sample count in file");
  }
  if (layout == 1) {
    handle.write(reinterpret_cast<char *>(&_n_samples), 4);
    // handle << _n_samples;
  }
  std::uint16_t tmp;
  tmp = varid.size();
  handle.write(reinterpret_cast<char *>(&tmp), 2);
  handle << varid;
  tmp = rsid.size();
  handle.write(reinterpret_cast<char *>(&tmp), 2);
  handle << rsid;
  tmp = chrom.size();
  handle.write(reinterpret_cast<char *>(&tmp), 2);
  handle << chrom;
  handle.write(reinterpret_cast<char *>(&pos), 4);
  
  if (layout != 1) {
    std::uint16_t n_alleles = alleles.size();
    handle.write(reinterpret_cast<char *>(&n_alleles), 2);
  }

  if ((layout == 1) && alleles.size() != 2) {
    throw std::invalid_argument("layout 1 requires exactly two alleles.");
  }

  std::uint32_t allele_size;
  for (auto &x : alleles) {
    allele_size = x.size();
    handle.write(reinterpret_cast<char *>(&allele_size), 4);
    handle << x;
  }
  handle.flush();
}

// uncompress a char array with zlib
void zlib_compress(char * input, int input_len, std::vector<char> &output) {
  z_stream strm;
  strm.zalloc = Z_NULL;
  strm.zfree = Z_NULL;
  strm.opaque = Z_NULL;

  strm.avail_in = input_len;      // size of input
  strm.next_in = (Bytef *) input; // input char array
  strm.avail_out = output.size();        // size of output
  strm.next_out = (Bytef *) &output[0]; // output char array

  deflateInit(&strm, 6);
  deflate(&strm, Z_FINISH);
  deflateEnd(&strm);

  output.resize(strm.total_out);
}

// uncompress a char array with zstd
void zstd_compress(char * input, int input_len, std::vector<char> &output) {
  std::size_t total_out = ZSTD_compress(&output[0], output.size(), input, input_len, 3);
  output.resize(total_out);
}
 
/// Read genotype data for a variant from disk and decompress.
///
/// The decompressed data is stored in the 'uncompressed' member. Decompression
/// is handled internally by either zlib_decompress, or zstd_decompress,
/// depending on compression scheme.
std::vector<char> compress(std::vector<std::uint8_t> &uncompressed, std::uint32_t compression) {
  std::vector<char> compressed(uncompressed.size() * 5 + 20);
  if (compression == 1) { // zlib
    zlib_compress(reinterpret_cast<char *>(&uncompressed[0]), (int)uncompressed.size(), compressed);
  } else if (compression == 2) { // zstd
    zstd_compress(reinterpret_cast<char *>(&uncompressed[0]), (int)uncompressed.size(), compressed);
  }
  return compressed;
}

bool missing_genotypes(double * genotypes, std::uint16_t len) {
  std::uint16_t nan_count = 0;
    for (std::uint32_t i=0; i<len; i++) {
      nan_count += std::isnan(genotypes[i]);
    }
    if ((nan_count > 0) && (nan_count < len)) {
      throw std::invalid_argument("samples with any missing genotype must encode all as missing (i.e. float(nan))");
    }
    return nan_count == len;
}

std::vector<std::uint8_t> encode_layout1(
                    double *genotypes,
                    std::uint32_t geno_len) {
  // genotypes are encoded as 16-bit uints, so resize to n_genotypes * 2
  std::vector<std::uint8_t> encoded(geno_len * 2 + 8);

  std::uint32_t i = 0;
  std::int32_t scaled32;
  std::uint16_t scaled;
  bool missing;
  double g;
  for (std::uint32_t j=0; j < geno_len; j+=3) {
    missing = missing_genotypes(&genotypes[j], 3);
    for (std::uint32_t k=0; k<3; k++) {
      g = genotypes[j + k];
      if (missing) {
        g = 0;
      }
      scaled32 = (std::int32_t)std::round(g * 32768);
      // check the value is in bounds
      if ((scaled32 < 0) || (scaled32 > 65535)) {
        throw std::invalid_argument("scaled genotype is out of bounds");
      }
      scaled = scaled32;
      std::memcpy(&encoded[i], &scaled, 2);
      i += 2;
    }
  }
  encoded.resize(geno_len * 2);
  return encoded;
}

std::vector<std::uint8_t> encode_layout2(
                    std::uint32_t n_samples,
                    std::uint16_t n_alleles,
                    double *genotypes,
                    std::uint32_t geno_len,
                    std::vector<uint8_t> &ploidy,
                    std::uint8_t min_ploidy,
                    std::uint8_t max_ploidy,
                    bool phased,
                    std::uint8_t &bit_depth
                    ) 
{
  int _max_ploid = (int)max_ploidy;
  int _n_alleles = (int)n_alleles;
  std::uint32_t max_probs = get_max_probs(_max_ploid, _n_alleles, phased);
  if ((geno_len / max_probs) != n_samples) {
    throw std::invalid_argument("genotypes and ploidy lengths don't match");
  }

  std::uint32_t probs_len = (n_samples * bit_depth) * (max_probs - 1);
  bool remainder = (probs_len % 8) > 0;
  probs_len = (probs_len / 8) + remainder;

  std::uint32_t encoded_size = 10 + n_samples + probs_len;
  std::vector<std::uint8_t> encoded(encoded_size + 8);  // extend slightly to help with variable bit depths
  std::uint32_t i=0;
  std::memcpy(&encoded[i], &n_samples, 4);
  i += 4;
  std::memcpy(&encoded[i], &n_alleles, 2);
  i += 2;
  encoded[i] = min_ploidy;
  i += 1;
  encoded[i] = max_ploidy;
  i += 1;

  // set the individuals ploidy values. We'll fill samples with missing data 
  // when we run through the genotypes
  const std::uint32_t ploidy_offset = i;
  if (min_ploidy == max_ploidy) {
    std::memset(&encoded[i], max_ploidy, n_samples);
    i += n_samples;
  } else {
    for (auto &x : ploidy) {
      encoded[i] = x;
      i += 1;
    }
  }

  encoded[i] = phased;
  i += 1;
  encoded[i] = bit_depth;
  i += 1;

  double factor = std::pow(2, bit_depth) - 1;
  bool missing;
  const std::uint32_t genotype_offset=i;
  std::uint32_t bit_idx=0;
  std::uint32_t byte_idx;
  std::uint32_t bit_remainder;
  std::uint64_t converted;
  std::uint64_t window;
  double sample_max;
  double g;
  for (std::uint32_t i=0; i<(n_samples*max_probs); i+= max_probs) {
    missing = missing_genotypes(&genotypes[i], max_probs);
    if (missing) {
      encoded[ploidy_offset + (i / 3)] |= 0x80;
    }
    sample_max = 0;
    for (std::uint32_t j = 0; j < (max_probs - 1); j++) {
      g = genotypes[i + j];
      if (missing) {
        g = 0;
      }
      sample_max = std::max(sample_max, g);
    }
    for (std::uint32_t j = 0; j < (max_probs - 1); j++) {
      g = genotypes[i + j];
      if (missing) {
        g = 0;
      }
      byte_idx = genotype_offset + (bit_idx / 8);
      bit_remainder = bit_idx % 8;
      window = *reinterpret_cast<const std::uint32_t* >(&encoded[byte_idx]);
      if (g == sample_max) {
        converted = std::ceil(g * factor);
      } else {
        converted = std::floor(g * factor);
      }
      window |= (converted << bit_remainder);

      std::memcpy(&encoded[byte_idx], &window, 8);
      bit_idx += bit_depth;
    }
  }
  encoded.resize(encoded_size);
  return encoded;
}

void BgenWriter::add_genotype_data(std::uint16_t n_alleles,
                                   double *genotypes,
                                   std::uint32_t geno_len,
                                   std::uint8_t ploidy,
                                   bool phased,
                                   std::uint8_t bit_depth)
{
  std::vector<std::uint8_t> ploidy_vector;
  add_genotype_data(n_alleles, genotypes, geno_len, ploidy_vector, ploidy, ploidy, phased, bit_depth);
}

void BgenWriter::add_genotype_data(std::uint16_t n_alleles,
                                   double *genotypes,
                                   std::uint32_t geno_len,
                                   std::vector<uint8_t> &ploidy,
                                   std::uint8_t min_ploidy,
                                   std::uint8_t max_ploidy,
                                   bool phased,
                                   std::uint8_t bit_depth)
{
  if ((layout == 1) && (compression == 2)) {
    throw std::invalid_argument("you cannot use zstd compression with layout 1");
  }

  std::vector<std::uint8_t> encoded;
  if (layout == 1) {
    encoded = encode_layout1(genotypes, geno_len);
  } else {
    encoded = encode_layout2(n_samples, n_alleles, genotypes, geno_len, ploidy, 
                   min_ploidy, max_ploidy, phased, bit_depth);
  }

  std::vector<char> compressed;
  if (compression != 0) {
    compressed = compress(encoded, compression);
  }
  std::uint32_t compressed_len = compressed.size();

  std::uint32_t size;
  if (layout == 1) {
    if (compression == 0) {
      for (auto &x : encoded) {
        handle << x;
      }
    } else {
      handle.write(reinterpret_cast<char *>(&compressed_len), 4);
      for (auto &x : compressed) {
        handle << x;
      }
    }
  } else if (layout == 2) {
    // layout 2
    if (compression == 0) {
      size = encoded.size();
      handle.write(reinterpret_cast<char *>(&size), 4);
      for (auto &x : encoded) {
        handle << x;
      }
    } else {
      size = compressed_len + 4;
      handle.write(reinterpret_cast<char *>(&size), 4);
      size = encoded.size();
      handle.write(reinterpret_cast<char *>(&size), 4);
      for (auto &x : compressed) {
        handle << x;
      }
    }
  } else {
    throw std::invalid_argument("layout must be 1 or 2");
  }
}

}  // namespace bgen