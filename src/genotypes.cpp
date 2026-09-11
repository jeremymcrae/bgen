
#include <algorithm>
#include <iterator>
#include <stdexcept>
#include <bitset>
#include <cmath>
#include <cassert>
#include <cstring>

#include <iostream>

#if defined(__x86_64__)
  #include <immintrin.h>
#endif

#if defined(__aarch64__)
  #include <arm_neon.h>
#endif

#include "zstd/lib/zstd.h"
#include "zlib.h"

#include "genotypes.h"
#include "utils.h"

namespace bgen {

// create lookup table of all probs for 8 bit integers. This has 511 entries, in
// order to allow for looking up values for the minor allele dosage, which can
// be up to 2.0 if a sample contains both copies are of the minor allele
// (2.0=255 + 255).
const float lut8[511] = {0.0000000, 0.0039216, 0.0078431, 0.0117647, 0.0156863, 0.0196078,
  0.0235294, 0.0274510, 0.0313725, 0.0352941, 0.0392157, 0.0431373, 0.0470588,
  0.0509804, 0.0549020, 0.0588235, 0.0627451, 0.0666667, 0.0705882, 0.0745098,
  0.0784314, 0.0823529, 0.0862745, 0.0901961, 0.0941176, 0.0980392, 0.1019608,
  0.1058824, 0.1098039, 0.1137255, 0.1176471, 0.1215686, 0.1254902, 0.1294118,
  0.1333333, 0.1372549, 0.1411765, 0.1450980, 0.1490196, 0.1529412, 0.1568627,
  0.1607843, 0.1647059, 0.1686275, 0.1725490, 0.1764706, 0.1803922, 0.1843137,
  0.1882353, 0.1921569, 0.1960784, 0.2000000, 0.2039216, 0.2078431, 0.2117647,
  0.2156863, 0.2196078, 0.2235294, 0.2274510, 0.2313725, 0.2352941, 0.2392157,
  0.2431373, 0.2470588, 0.2509804, 0.2549020, 0.2588235, 0.2627451, 0.2666667,
  0.2705882, 0.2745098, 0.2784314, 0.2823529, 0.2862745, 0.2901961, 0.2941176,
  0.2980392, 0.3019608, 0.3058824, 0.3098039, 0.3137255, 0.3176471, 0.3215686,
  0.3254902, 0.3294118, 0.3333333, 0.3372549, 0.3411765, 0.3450980, 0.3490196,
  0.3529412, 0.3568627, 0.3607843, 0.3647059, 0.3686275, 0.3725490, 0.3764706,
  0.3803922, 0.3843137, 0.3882353, 0.3921569, 0.3960784, 0.4000000, 0.4039216,
  0.4078431, 0.4117647, 0.4156863, 0.4196078, 0.4235294, 0.4274510, 0.4313725,
  0.4352941, 0.4392157, 0.4431373, 0.4470588, 0.4509804, 0.4549020, 0.4588235,
  0.4627451, 0.4666667, 0.4705882, 0.4745098, 0.4784314, 0.4823529, 0.4862745,
  0.4901961, 0.4941176, 0.4980392, 0.5019608, 0.5058824, 0.5098039, 0.5137255,
  0.5176471, 0.5215686, 0.5254902, 0.5294118, 0.5333333, 0.5372549, 0.5411765,
  0.5450980, 0.5490196, 0.5529412, 0.5568627, 0.5607843, 0.5647059, 0.5686275,
  0.5725490, 0.5764706, 0.5803922, 0.5843137, 0.5882353, 0.5921569, 0.5960784,
  0.6000000, 0.6039216, 0.6078431, 0.6117647, 0.6156863, 0.6196078, 0.6235294,
  0.6274510, 0.6313725, 0.6352941, 0.6392157, 0.6431373, 0.6470588, 0.6509804,
  0.6549020, 0.6588235, 0.6627451, 0.6666667, 0.6705882, 0.6745098, 0.6784314,
  0.6823529, 0.6862745, 0.6901961, 0.6941176, 0.6980392, 0.7019608, 0.7058824,
  0.7098039, 0.7137255, 0.7176471, 0.7215686, 0.7254902, 0.7294118, 0.7333333,
  0.7372549, 0.7411765, 0.7450980, 0.7490196, 0.7529412, 0.7568627, 0.7607843,
  0.7647059, 0.7686275, 0.7725490, 0.7764706, 0.7803922, 0.7843137, 0.7882353,
  0.7921569, 0.7960784, 0.8000000, 0.8039216, 0.8078431, 0.8117647, 0.8156863,
  0.8196078, 0.8235294, 0.8274510, 0.8313725, 0.8352941, 0.8392157, 0.8431373,
  0.8470588, 0.8509804, 0.8549020, 0.8588235, 0.8627451, 0.8666667, 0.8705882,
  0.8745098, 0.8784314, 0.8823529, 0.8862745, 0.8901961, 0.8941176, 0.8980392,
  0.9019608, 0.9058824, 0.9098039, 0.9137255, 0.9176471, 0.9215686, 0.9254902,
  0.9294118, 0.9333333, 0.9372549, 0.9411765, 0.9450980, 0.9490196, 0.9529412,
  0.9568627, 0.9607843, 0.9647059, 0.9686275, 0.9725490, 0.9764706, 0.9803922,
  0.9843137, 0.9882353, 0.9921569, 0.9960784, 1.0000000, 1.0039216, 1.0078431,
  1.0117647, 1.0156863, 1.0196078, 1.0235294, 1.0274510, 1.0313725, 1.0352941,
  1.0392157, 1.0431373, 1.0470588, 1.0509804, 1.0549020, 1.0588235, 1.0627451,
  1.0666667, 1.0705882, 1.0745098, 1.0784314, 1.0823529, 1.0862745, 1.0901961,
  1.0941176, 1.0980392, 1.1019608, 1.1058824, 1.1098039, 1.1137255, 1.1176471,
  1.1215686, 1.1254902, 1.1294118, 1.1333333, 1.1372549, 1.1411765, 1.1450980,
  1.1490196, 1.1529412, 1.1568627, 1.1607843, 1.1647059, 1.1686275, 1.1725490,
  1.1764706, 1.1803922, 1.1843137, 1.1882353, 1.1921569, 1.1960784, 1.2000000,
  1.2039216, 1.2078431, 1.2117647, 1.2156863, 1.2196078, 1.2235294, 1.2274510,
  1.2313725, 1.2352941, 1.2392157, 1.2431373, 1.2470588, 1.2509804, 1.2549020,
  1.2588235, 1.2627451, 1.2666667, 1.2705882, 1.2745098, 1.2784314, 1.2823529,
  1.2862745, 1.2901961, 1.2941176, 1.2980392, 1.3019608, 1.3058824, 1.3098039,
  1.3137255, 1.3176471, 1.3215686, 1.3254902, 1.3294118, 1.3333333, 1.3372549,
  1.3411765, 1.3450980, 1.3490196, 1.3529412, 1.3568627, 1.3607843, 1.3647059,
  1.3686275, 1.3725490, 1.3764706, 1.3803922, 1.3843137, 1.3882353, 1.3921569,
  1.3960784, 1.4000000, 1.4039216, 1.4078431, 1.4117647, 1.4156863, 1.4196078,
  1.4235294, 1.4274510, 1.4313725, 1.4352941, 1.4392157, 1.4431373, 1.4470588,
  1.4509804, 1.4549020, 1.4588235, 1.4627451, 1.4666667, 1.4705882, 1.4745098,
  1.4784314, 1.4823529, 1.4862745, 1.4901961, 1.4941176, 1.4980392, 1.5019608,
  1.5058824, 1.5098039, 1.5137255, 1.5176471, 1.5215686, 1.5254902, 1.5294118,
  1.5333333, 1.5372549, 1.5411765, 1.5450980, 1.5490196, 1.5529412, 1.5568627,
  1.5607843, 1.5647059, 1.5686275, 1.5725490, 1.5764706, 1.5803922, 1.5843137,
  1.5882353, 1.5921569, 1.5960784, 1.6000000, 1.6039216, 1.6078431, 1.6117647,
  1.6156863, 1.6196078, 1.6235294, 1.6274510, 1.6313725, 1.6352941, 1.6392157,
  1.6431373, 1.6470588, 1.6509804, 1.6549020, 1.6588235, 1.6627451, 1.6666667,
  1.6705882, 1.6745098, 1.6784314, 1.6823529, 1.6862745, 1.6901961, 1.6941176,
  1.6980392, 1.7019608, 1.7058824, 1.7098039, 1.7137255, 1.7176471, 1.7215686,
  1.7254902, 1.7294118, 1.7333333, 1.7372549, 1.7411765, 1.7450980, 1.7490196,
  1.7529412, 1.7568627, 1.7607843, 1.7647059, 1.7686275, 1.7725490, 1.7764706,
  1.7803922, 1.7843137, 1.7882353, 1.7921569, 1.7960784, 1.8000000, 1.8039216,
  1.8078431, 1.8117647, 1.8156863, 1.8196078, 1.8235294, 1.8274510, 1.8313725,
  1.8352941, 1.8392157, 1.8431373, 1.8470588, 1.8509804, 1.8549020, 1.8588235,
  1.8627451, 1.8666667, 1.8705882, 1.8745098, 1.8784314, 1.8823529, 1.8862745,
  1.8901961, 1.8941176, 1.8980392, 1.9019608, 1.9058824, 1.9098039, 1.9137255,
  1.9176471, 1.9215686, 1.9254902, 1.9294118, 1.9333333, 1.9372549, 1.9411765,
  1.9450980, 1.9490196, 1.9529412, 1.9568627, 1.9607843, 1.9647059, 1.9686275,
  1.9725490, 1.9764706, 1.9803922, 1.9843137, 1.9882353, 1.9921569, 1.9960784,
  2.0000000};

// uncompress a char array with zlib
static void zlib_uncompress(char * input, int compressed_len, char * decompressed, int decompressed_len) {
  z_stream infstream;
  infstream.zalloc = Z_NULL;
  infstream.zfree = Z_NULL;
  infstream.opaque = Z_NULL;
  
  infstream.avail_in = compressed_len; // size of input
  infstream.next_in = (Bytef *) input; // input char array
  infstream.avail_out = decompressed_len; // size of output
  infstream.next_out = (Bytef *) decompressed; // output char array
  
  inflateInit(&infstream);
  inflate(&infstream, Z_NO_FLUSH);
  inflateEnd(&infstream);
  
  if (decompressed_len != (int) infstream.total_out) {
    throw std::invalid_argument("zlib decompression gave data of wrong length");
  }
}

// uncompress a char array with zstd
//
// ZSTD_decompress allocates and frees a context on every call, at a fixed cost which
// dominates a bgen of many small variants (~30% faster reusing one across 20000 variants
// of 500 samples, against ~7% across 300 variants of 20000 samples). The context is
// thread_local, since a Genotypes exists per variant, and a shared one is not thread safe
static ZSTD_DCtx * borrow_zstd_dctx() {
  // the unique_ptr frees the context when the thread exits
  struct Deleter {
    void operator()(ZSTD_DCtx * ctx) const { ZSTD_freeDCtx(ctx); }
  };
  static thread_local std::unique_ptr<ZSTD_DCtx, Deleter> ctx;
  if (!ctx) {
    ctx.reset(ZSTD_createDCtx());
    if (!ctx) {
      throw std::runtime_error("cannot allocate a zstd decompression context");
    }
  }
  return ctx.get();
}

static void zstd_uncompress(char * input, int compressed_len, char * decompressed,  int decompressed_len) {
  std::size_t total_out = ZSTD_decompressDCtx(borrow_zstd_dctx(), decompressed,
                                              decompressed_len, input, compressed_len);
  if (ZSTD_isError(total_out)) {
    throw std::invalid_argument(std::string("zstd decompression failed: ") +
                                ZSTD_getErrorName(total_out));
  }
  if (decompressed_len != (int) total_out) {
    throw std::invalid_argument("zstd decompression gave data of wrong length");
  }
}

/// Read a variant's genotype data block from the bgen, without decompression. The
/// block is padded like the decompressed buffer, so an uncompressed one needs no copy
void Genotypes::load_block() {
  if (block) {
    // the block is always allocated with padding, so it is never a null pointer
    // after it has been read, even for a zero length block
    return;
  }
  
  // a closed bgen is marked with the badbit, so refuse to read from it
  if (handle->bad()) {
    throw std::invalid_argument("cannot read from closed bgen file");
  }
  
  if (!is_stdin) {
    // any other error state (e.g. the failbit and eofbit left behind by a read
    // which ran to the end of the file) is recoverable, since we are about to
    // seek to a known good offset. Only safe when we can seek, so stdin keeps
    // whatever state it is in and fails the read below instead.
    handle->clear();
    handle->seekg(file_offset);  // about 1 microsecond
  }
  
  // the block length is read from the bgen, so guard the buffer size against
  // wrapping, as for the decompressed length below
  if (length > (UINT32_MAX - PROBS_READ_PAD)) {
    throw std::invalid_argument("bgen genotype data claims an implausible "
                                "length: " + std::to_string(length));
  }
  // hold the buffer in a unique_ptr, so a read which throws on a malformed bgen
  // cannot leak it
  std::unique_ptr<char[]> buffer(new char[length + PROBS_READ_PAD]);
  std::memset(buffer.get() + length, 0, PROBS_READ_PAD);
  if (! handle->read(buffer.get(), length)) {
    throw std::invalid_argument("couldn't read the compressed data");
  }
  block = std::move(buffer);
}

/// Decompress a variant's genotype data.
///
/// The decompressed data is stored in the 'uncompressed' member. Decompression
/// is handled internally by either zlib_decompress, or zstd_decompress,
/// depending on compression scheme.
void Genotypes::decompress() {
  if (is_decompressed) {
    // don't decompress if already available
    return;
  }
  
  load_block();
  
  bool decompressed_field = false;
  std::uint32_t decompressed_len = length;
  if (compression != 0) {
    if (layout == 1) {
      decompressed_len = n_samples * 6;
    } else if (layout == 2) {
      decompressed_field = true;
      if (length < sizeof(std::uint32_t)) {
        // the field is part of the block, so a block this short would underflow
        // the compressed length below
        throw std::invalid_argument("bgen genotype data is too short to hold a "
                                    "decompressed length");
      }
      std::memcpy(&decompressed_len, block.get(), sizeof(std::uint32_t));
    }
  }
  
  if (compression == 0) {
    // the block is the genotype data, and was padded when it was read, so hand it
    // over as it stands rather than copying it into a second buffer
    uncompressed = std::move(block);
    uncompressed_len = length;
    is_decompressed = true;
    return;
  }
  
  std::uint32_t compressed_len = length - decompressed_field * 4;
  char * compressed = block.get() + decompressed_field * 4;
  // the decompressed length is read from the bgen, so guard against wrapping, since a
  // length near the 32-bit limit would allocate a few bytes then memset past the end
  if (decompressed_len > (UINT32_MAX - PROBS_READ_PAD)) {
    throw std::invalid_argument("bgen genotype data claims an implausible "
                                "decompressed length: " +
                                std::to_string(decompressed_len));
  }
  // zero the padding, so the trailing bits probabilities_layout2 reads are deterministic
  std::unique_ptr<char[]> buffer(new char[decompressed_len + PROBS_READ_PAD]);
  std::memset(buffer.get() + decompressed_len, 0, PROBS_READ_PAD);
  
  if (compression == 1) { // zlib
    zlib_uncompress(compressed, (int) compressed_len, buffer.get(), (int) decompressed_len);  // about 2 milliseconds
  } else if (compression == 2) { // zstd
    zstd_uncompress(compressed, (int) compressed_len, buffer.get(), (int) decompressed_len);
  }
  // only take ownership once the data has decompressed cleanly, so a failed
  // parse leaves no half filled buffer behind for a later call to read from.
  // The block is kept until then, so a retry decompresses the same bytes.
  uncompressed = std::move(buffer);
  uncompressed_len = decompressed_len;
  is_decompressed = true;
  block.reset();
}

/// figure out the maximum number of probabilities across the individuals
///
/// @param max_ploidy maximum ploidy value across all samples
/// @param n_alleles number of alleles for the variant
/// @param phased whether the genotypes are phased
/// @return integer for maximum number of possible probabilites per sample
std::uint32_t get_max_probs(int & max_ploidy, int & n_alleles, bool & phased) {
  std::uint32_t max_probs;
  if (phased) {
    max_probs = n_alleles;
  } else {
    max_probs = n_choose_k(max_ploidy + n_alleles - 1, n_alleles - 1);
  }
  return max_probs;
}

/// @brief number of bytes the packed probability data needs for this variant
///
/// Every sample stores all but one probability per group, at bit_depth bits each, and
/// phased data stores a group per haplotype. The counts come from the file, so a 32-bit
/// total could overflow
std::uint64_t Genotypes::probability_bytes() {
  if (layout == 1) {
    // layout 1 is a fixed three 16-bit probabilities per sample, with nothing
    // omitted, and decompress() already sized the buffer as n_samples * 6
    return (std::uint64_t) n_samples * 6;
  }
  
  std::uint64_t values = 0;
  if (constant_ploidy) {
    std::uint64_t groups = phased ? (std::uint64_t) max_ploidy : 1;
    std::uint64_t per_group = phased ? (std::uint64_t) n_alleles : max_probs;
    if (per_group < 1) {
      throw std::invalid_argument("bgen variant stores no probabilities per sample");
    }
    values = (std::uint64_t) n_samples * groups * (per_group - 1);
  } else if (phased || (n_alleles == 2)) {
    // both come to one stored value per (haplotype x allele), so the total is just the
    // sum of the ploidy values, which vectorises. parse_ploidy filled the array, since
    // this is only reached when the ploidy varies
    std::uint64_t total_ploidy = fast_ploidy_sum(ploidy.get(), n_samples);
    values = phased ? total_ploidy * ((std::uint64_t) n_alleles - 1) : total_ploidy;
  } else {
    // unphased, multi-allelic and varying ploidy, so each sample stores a different
    // number of genotypes. The ploidy is masked to 6 bits, so 64 counts cover them all
    //
    // Only fill up to max_ploidy, because higher counts cannot occur and may
    // cause problems e.g. 12 alleles at a ploidy of 63 overflows a 32-bit count
    std::uint32_t counts[64];
    for (int ploid=0; ploid <= max_ploidy; ploid++) {
      int p = ploid;
      std::uint32_t n_probs = get_max_probs(p, n_alleles, phased);
      if (n_probs < 1) {
        throw std::invalid_argument("bgen variant stores no probabilities per sample");
      }
      counts[ploid] = n_probs - 1;
    }
    for (std::uint32_t n=0; n < n_samples; n++) {
      values += counts[ploidy[n]];
    }
  }
  std::uint64_t bits = values * bit_depth;
  return (bits / 8) + ((bits % 8) > 0);
}

/// @brief check the genotype block is big enough for the variant it describes
///
/// The block length and the fields sizing the data it needs are stored separately, and
/// the parsing indexes on those fields, so a malformed bgen would read past the buffer
void Genotypes::check_block_size() {
  // the accessors in variant.cpp only parse the header when max_probs is zero, so clear
  // it on failure, or a rejected variant hands back data on the next call
  std::uint64_t required;
  try {
    required = (std::uint64_t) idx + probability_bytes();
  } catch (...) {
    max_probs = 0;
    throw;
  }
  if (required > (std::uint64_t) uncompressed_len) {
    max_probs = 0;
    throw std::invalid_argument("bgen variant is truncated - its genotype data "
                                "is " + std::to_string(uncompressed_len) +
                                " bytes, but the variant needs " +
                                std::to_string(required));
  }
  block_checked = true;
}

/// parse the initial data that defines the ploidy and phased status.
///
/// Layout 1 doesn't store any information before the genotype probabilities,
/// so layout 1 just receives default values.
void Genotypes::load_data_and_parse_header() {
  decompress();
  if (has_ploidy) {
    if (!block_checked) {
      // the header was parsed on an earlier call which then failed the size
      // check, so fail the same way rather than fall through to the reads
      check_block_size();
    }
    return;
  }
  idx = 0;
  if (layout == 1) {
    phased = false;
    min_ploidy = 2;
    max_ploidy = 2;
    bit_depth = 16;
  } else if (layout == 2) {
    // the preamble is 10 bytes plus a ploidy byte per sample, and the reads below are
    // sized from the bgen header, so check it is present before reading any of it
    const std::uint64_t preamble = 10 + (std::uint64_t) n_samples;
    if (preamble > (std::uint64_t) uncompressed_len) {
      throw std::invalid_argument("bgen variant is truncated - its genotype data "
                                  "is " + std::to_string(uncompressed_len) +
                                  " bytes, but " + std::to_string(preamble) +
                                  " are needed for the ploidy of " +
                                  std::to_string(n_samples) + " samples");
    }
    if (n_alleles < 1) {
      // n_choose_k below loops on n_alleles - 1, so zero alleles would run it
      // about four billion times before returning a meaningless count
      throw std::invalid_argument("bgen variant has no alleles");
    }
    std::uint32_t nn_samples = *reinterpret_cast<const std::uint32_t*>(&uncompressed[idx]);
    idx += sizeof(std::uint32_t);
    std::uint16_t allele_check = *reinterpret_cast<const std::uint16_t*>(&uncompressed[idx]);
    idx += sizeof(std::uint16_t);
    if (nn_samples != (std::uint32_t) n_samples) {
      throw std::invalid_argument("number of samples doesn't match!");
    }
    if (allele_check != n_alleles) {
      throw std::invalid_argument("number of alleles doesn't match!");
    }
    
    min_ploidy = (int) *reinterpret_cast<const std::uint8_t*>(&uncompressed[idx]);
    idx += sizeof(std::uint8_t);
    max_ploidy = (int) *reinterpret_cast<const std::uint8_t*>(&uncompressed[idx]);
    idx += sizeof(std::uint8_t);

    if (max_ploidy < min_ploidy) {
      throw std::invalid_argument("maximum ploidy is less than minimum ploidy");
    } else if (max_ploidy > 63) {
      throw std::invalid_argument("maximum ploidy is greater than 63");
    }
  }
  
  constant_ploidy = (min_ploidy == max_ploidy);
  parse_ploidy();
  
  if (layout == 2) {
    phased = (bool) *reinterpret_cast<const std::uint8_t*>(&uncompressed[idx]);
    idx += sizeof(std::uint8_t);
    bit_depth = (int) *reinterpret_cast<const std::uint8_t*>(&uncompressed[idx]);
    if ((bit_depth < 1) | (bit_depth > 32)) {
      throw std::invalid_argument("probabilities bit depth out of bounds");
    }
    idx += sizeof(std::uint8_t);
  }
  max_probs = get_max_probs(max_ploidy, n_alleles, phased);
  check_block_size();
}

/// get ploidy state for all samples (and missingness for layout2).
///
/// ploidy state is stored in ploidy member. Missingness is stored as indices of
/// samples This takes 100 microseconds for 500k samples. Layout 1 variants lack
/// ploidy info, and are assigned the same ploidy for all samples.
///
/// A constant ploidy is left for materialise_ploidy() to fill, since max_ploidy already
/// says it. Only the missingness scan has to happen now, while the block is here
void Genotypes::parse_ploidy() {
  if (has_ploidy) {
    return;
  }
  
  has_ploidy = true;
  if (layout == 1) {
    // no per-sample ploidy is stored, so there is nothing to scan and nothing to fill
    // beyond the constant, which materialise_ploidy handles
    return;
  }
  
  // we want to avoid parsing the ploidy states if  every sample has the same
  // ploidy. If we have a constant ploidy, set all entries to the same value
  std::uint8_t mask = 63;
  if (constant_ploidy) {
    // only the missingness matters, since the ploidy is known, and the scan is vectorised
    fast_missing_scan(&uncompressed[idx], n_samples, missing);
  } else {
    ploidy.reset(new std::uint8_t[n_samples]);
    for (std::uint32_t x=0; x < n_samples; x++) {
      ploidy[x] = mask & uncompressed[idx + x];
      if (uncompressed[idx + x] & 0x80) {
        missing.push_back(x);
      }
    }
    // the ploidy sizes the probability reads, but the arrays they fill are sized from
    // max_ploidy, so a sample outside the declared bounds writes past the end of its row
    Range range = fast_range(ploidy.get(), n_samples);
    if ((range._max > (std::uint8_t) max_ploidy) ||
        (range._min < (std::uint8_t) min_ploidy)) {
      // this throws partway through the header parse, so clear the state a later call
      // would take as parsed, or it carries on with a half parsed header
      has_ploidy = false;
      missing.clear();
      throw std::invalid_argument("bgen variant has sample ploidies from " +
                                  std::to_string((int) range._min) + " to " +
                                  std::to_string((int) range._max) +
                                  ", outside the " + std::to_string(min_ploidy) +
                                  " to " + std::to_string(max_ploidy) +
                                  " range the variant declares");
    }
  }
  idx += n_samples;
}

/// @brief fill the per-sample ploidy array, if it is not already filled
///
/// parse_ploidy leaves it empty when every sample shares a ploidy, since filling half a
/// million entries costs more than the decode. Prefer max_ploidy to forcing this
void Genotypes::materialise_ploidy() {
  if (ploidy) {
    return;
  }
  ploidy.reset(new std::uint8_t[n_samples]);
  std::memset(ploidy.get(), max_ploidy, n_samples);
}

/// parse probabilities for layout1.
///
/// Layout 1 genotype probabilities are simple to parse, they are just 16-bit ints
/// one for each of the three possible genotypes. Missingness is encoded for a
/// sample when all probabilties are zero.
///
/// @return 1D float array of genotype probabilties (each from 0.0-1.0).
void Genotypes::probabilities_layout1(char * uncompressed, std::uint32_t idx, float * probs, std::uint32_t & nrows) {
  float factor = 1.0 / 32768;
  for (std::uint32_t offset=0; offset<nrows * max_probs; offset+=max_probs) {
    probs[offset] = *reinterpret_cast<const std::uint16_t*>(&uncompressed[idx]) * factor;
    probs[offset + 1] = *reinterpret_cast<const std::uint16_t*>(&uncompressed[idx + 2]) * factor;
    probs[offset + 2] = *reinterpret_cast<const std::uint16_t*>(&uncompressed[idx + 4]) * factor;
    idx += 6;
    
    if ((probs[offset] == 0.0) & (probs[offset + 1] == 0.0) & (probs[offset + 2] == 0.0)) {
      probs[offset] = std::nan("1");
      probs[offset + 1] = std::nan("1");
      probs[offset + 2] = std::nan("1");
    }
  }
}

#if defined(__x86_64__)
// the AVX2 half of fast_haplotype_probs below. Compiled for AVX2 in isolation,
// so it is only ever entered after a runtime check for AVX2 support.
BGEN_TARGET_AVX2
static void haplotype_probs_avx2(char * uncompressed, std::uint32_t & idx,
                                 float * probs, std::uint32_t & nrows,
                                 std::uint32_t & n) {
  const std::uint32_t c = 255;
  const float f = 1.0 / 255.0f;
  __m256i k = _mm256_set_epi32(c, c, c, c, c, c, c, c);
  __m256 j = _mm256_set_ps(f, f, f, f, f, f, f, f);

  __m128i initial;
  __m256i widened, partner;
  __m256 store_lo, store_hi;

  // run through most of the samples, but make sure we stay well away from the
  // end of the float array, so the mm_loadu doesn't attempt to load beyond
  // the end of the array
  for ( ; n + 32 < nrows * 2; n+=32) {
    // load 16 values into a m128 register
    initial = _mm_loadu_si128((const __m128i*) &uncompressed[idx]);

    // expand the first 8 values to 32-bits, and get the corresponding delta
    widened = _mm256_cvtepu8_epi32(initial);
    partner = _mm256_sub_epi32(k, widened);

    // interleave the registers to get paired values beside each other, then
    // convert to float and multiply by inverted constant to get probability
    store_lo = _mm256_mul_ps(_mm256_cvtepi32_ps(_mm256_unpacklo_epi32(widened, partner)), j);
    store_hi = _mm256_mul_ps(_mm256_cvtepi32_ps(_mm256_unpackhi_epi32(widened, partner)), j);

    // _mm256_unpacklo_epi32 interleaves values from two registers, but it
    // does it within the top part of each m128 portion of the register
    // we start with
    // widened = [a1, b1, c1, d1, e1, f1, g1, h1]
    // partner = [a2, b2, c2, d2, e2, f2, g2, h2]  # (1 - widened)

    // we want to store the sequential partners e.g.
    // [a1, a2, b1, b2, c1, c2, d1, d1]

    // but _mm256_unpacklo_epi32 gives registers with mixed order
    // [a1, a2, b1, b2, e1, e2, f1, f2]
    // [c1, c2, d1, d2, g1, g2, h1, h2]
    // this means we need to store the first and last four values in separate
    // steps to ensure the correct ordering
    _mm_storeu_ps(&probs[n], _mm256_castps256_ps128(store_lo));
    _mm_storeu_ps(&probs[n + 4], _mm256_castps256_ps128(store_hi));
    _mm_storeu_ps(&probs[n + 8], _mm256_extractf128_ps(store_lo, 1));
    _mm_storeu_ps(&probs[n + 12], _mm256_extractf128_ps(store_hi, 1));

    // now repeat for the other 8 values of the original register
    widened = _mm256_cvtepu8_epi32(_mm_bsrli_si128(initial, 8));
    partner = _mm256_sub_epi32(k, widened);

    store_lo = _mm256_mul_ps(_mm256_cvtepi32_ps(_mm256_unpacklo_epi32(widened, partner)), j);
    store_hi = _mm256_mul_ps(_mm256_cvtepi32_ps(_mm256_unpackhi_epi32(widened, partner)), j);

    _mm_storeu_ps(&probs[n + 16], _mm256_castps256_ps128(store_lo));
    _mm_storeu_ps(&probs[n + 20], _mm256_castps256_ps128(store_hi));
    _mm_storeu_ps(&probs[n + 24], _mm256_extractf128_ps(store_lo, 1));
    _mm_storeu_ps(&probs[n + 28], _mm256_extractf128_ps(store_hi, 1));

    idx += 16;
  }
}

/// @brief AVX2 half of the unphased biallelic diploid 8 bit decode
///
/// All three probabilities come out of lut8 by gather rather than arithmetic, since its
/// seven digit literals differ from value * (1.0f/255.0f) in 365 of 511 entries, so
/// computing them would change every value the library returns. Leaves n as the samples
/// handled, and sets above_max if any sample's stored values summed above the maximum
BGEN_TARGET_AVX2
static void unphased_probs_avx2(char * uncompressed, std::uint32_t & idx,
                                float * probs, std::uint32_t & nrows,
                                std::uint32_t & n, bool & above_max) {
  const __m256i k255 = _mm256_set1_epi32(255);
  const __m256i zero = _mm256_setzero_si256();
  __m256i signs = _mm256_setzero_si256();
  // gather the two stored bytes of 8 samples from the 16 bytes they occupy
  const __m128i even = _mm_setr_epi8(0, 2, 4, 6, 8, 10, 12, 14,
                                     -1, -1, -1, -1, -1, -1, -1, -1);
  const __m128i odd = _mm_setr_epi8(1, 3, 5, 7, 9, 11, 13, 15,
                                    -1, -1, -1, -1, -1, -1, -1, -1);

  const std::uint8_t * in = reinterpret_cast<const std::uint8_t *>(uncompressed + idx);
  float * out = probs;
  for ( ; n + 8 <= nrows; n += 8) {
    // the samples are stored as interleaved pairs, so split them into the first and
    // second probability of each sample before widening to the 32 bit gather indices
    __m128i raw = _mm_loadu_si128((const __m128i *) in);
    __m256i first = _mm256_cvtepu8_epi32(_mm_shuffle_epi8(raw, even));
    __m256i second = _mm256_cvtepu8_epi32(_mm_shuffle_epi8(raw, odd));

    // the third probability is the remainder, 255 - first - second. A bgen may store
    // values summing above 255, which would index before the table, so clamp at zero
    // exactly as remainder_lut8 does
    __m256i third = _mm256_sub_epi32(_mm256_sub_epi32(k255, first), second);
    // keep the sign bits before clamping, since the clamp is what hides a malformed sample
    signs = _mm256_or_si256(signs, third);
    third = _mm256_max_epi32(third, zero);

    __m256 p0 = _mm256_i32gather_ps(lut8, first, 4);
    __m256 p1 = _mm256_i32gather_ps(lut8, second, 4);
    __m256 p2 = _mm256_i32gather_ps(lut8, third, 4);

    // AVX2 has no scatter, and the output stride of three floats per sample does not
    // divide the register width, so spill the lanes and interleave them
    float b0[8], b1[8], b2[8];
    _mm256_storeu_ps(b0, p0);
    _mm256_storeu_ps(b1, p1);
    _mm256_storeu_ps(b2, p2);
    for (int j = 0; j < 8; j++) {
      out[j * 3] = b0[j];
      out[j * 3 + 1] = b1[j];
      out[j * 3 + 2] = b2[j];
    }
    in += 16;
    out += 24;
  }
  // a negative remainder left its sign bit in the accumulator
  if (_mm256_movemask_ps(_mm256_castsi256_ps(signs)) != 0) {
    above_max = true;
  }
  idx += n * 2;
}
#endif

/// @brief look up the final probability of an unphased biallelic diploid sample
///
/// The remainder can index before the table, since a bgen may store values summing above
/// the maximum. Clamping keeps the read inside the table, and sets above_max, on the
/// integers already in registers, which cannot misfire on rounding as a float test could
static inline float remainder_lut8(std::uint8_t first, std::uint8_t second,
                                   bool & above_max) {
  int idx = 255 - (int) first - (int) second;
  above_max |= (idx < 0);
  return lut8[(idx > 0) ? idx : 0];
}

/// @brief the dosage count of an unphased biallelic diploid sample
///
/// Malformed files can push first * 2 + second past the 510 the spec allows, so clamp as
/// the SIMD paths do. Returning the count keeps the divide by 255 identical to those
static inline int dosage_count8(std::uint8_t first, std::uint8_t second,
                                bool & above_max) {
  int count = (int) first * 2 + (int) second;
  above_max |= (count > 510);
  return (count < 511) ? count : 510;
}

/// fast path for phased data with ploidy=2, and 8 bits per probability
void Genotypes::fast_haplotype_probs(char * uncompressed, std::uint32_t idx, float * probs,  std::uint32_t & nrows) {
  std::uint32_t n = 0;
#if defined(__x86_64__)
  if (__builtin_cpu_supports("avx2")) {
    haplotype_probs_avx2(uncompressed, idx, probs, nrows, n);
  }
#endif
  // finish off the unvectorized samples
  std::uint8_t first;
  for (; n < nrows * 2; n += 2) {
    first = *reinterpret_cast<const std::uint8_t*>(&uncompressed[idx]);
    probs[n] = lut8[first];
    probs[n + 1] = lut8[255 - first];
    idx += 1;
  }
}

/// parse probabilities for layout2
///
/// This is a fairly complex function, due to how the probabilties are encoded.
/// We need to extract probabilites per sample, but the final probability is
/// inferred as the remainder after all other genotype probabilties are parsed.
/// The number of bits used for encoding probabilties can differ between variants,
/// so there's a fast method for parsing 8-bit encoded probabilties, since that
/// is a common use case.
///
/// @return 1D float array of genotype probabilties (each from 0.0-1.0).
void Genotypes::probabilities_layout2(char * uncompressed, std::uint32_t idx, float * probs, std::uint32_t & nrows) {
  // get genotype/allele probabilities
  std::uint32_t n_probs;
  std::uint32_t max_less_1 = max_probs - 1;
  float prob = 0;
  float remainder;
  
  // define variables for parsing depths not aligned with 8 bit char array
  float factor = 1.0 / ((float) (std::pow(2, (int) bit_depth)) - 1);
  std::uint64_t probs_mask = std::uint64_t(0xFFFFFFFFFFFFFFFF) >> (64 - bit_depth);
  std::uint64_t bit_idx = 0;  // index position in bits
  
  if (constant_ploidy & (max_probs == 3) & (bit_depth == 8)) {
    // A fast path for one scenario: all samples have ploidy=2, with 8 bits per
    // probability. This optimises memory accesses, and avoids looping over a
    // sample. This is ~2.5X faster than the standard route to compute
    // probabilities, and diploid samples with 8 bits/prob is likely the most
    // common use case, so the speed-up justifies this special case.
    std::uint32_t n = 0;
#if defined(__x86_64__)
    if (__builtin_cpu_supports("avx2")) {
      unphased_probs_avx2(uncompressed, idx, probs, nrows, n, probs_above_max);
    }
#endif
    // finish off the samples the vector loop did not cover
    const std::uint8_t * in = reinterpret_cast<const std::uint8_t *>(&uncompressed[idx]);
    for (std::uint32_t offset = n * 3; offset < nrows * 3; offset += 3) {
      std::uint8_t first = in[0];
      std::uint8_t second = in[1];
      probs[offset] = lut8[first];
      probs[offset + 1] = lut8[second];
      probs[offset + 2] = remainder_lut8(first, second, probs_above_max);
      in += 2;
    }
  } else if (constant_ploidy & (max_probs == 2) & (bit_depth == 8)) {
    // fast path for phased data with ploidy=2, and 8 bits per probability
    fast_haplotype_probs(uncompressed, idx, probs, nrows);
  } else {
    for (std::uint32_t offset=0; offset < (nrows * max_probs); offset += max_probs) {
      // calculate the number of probabilities per sample (depends on whether the
      // data is phased, the sample ploidy and the number of alleles)
      if (constant_ploidy) {
        n_probs = max_less_1;
      } else if (phased) {
        n_probs = n_alleles - 1;
      } else if ((ploidy[offset / max_probs] == 2) && (n_alleles == 2)) {
        n_probs = 2;
      } else {
        n_probs = n_choose_k(ploidy[offset / max_probs] + n_alleles - 1, n_alleles - 1) - 1;
      }
      remainder = 1.0;
      for (std::uint32_t x=0; x<n_probs; x++) {
        prob = ((*reinterpret_cast<const std::uint64_t* >(&uncompressed[idx + bit_idx / 8]) >> bit_idx % 8) & probs_mask) * factor ;
        bit_idx += bit_depth;
        remainder -= prob;
        probs[offset + x] = prob;
      }
      // A negative remainder means the sample's stored probabilities sum above the maximum,
      // so the bgen is malformed. Clamp values to zero, and flag the clamping.
      if (remainder < -1e-6f) {
        remainder = 0.0f;
        probs_above_max = true;
      }
      probs[offset + n_probs] = remainder;
      for (std::uint32_t x=(n_probs + 1); x<max_probs; x++) {
        probs[offset + x] = std::nan("1");
      }
    }
  }
  
  std::uint32_t offset;
  std::uint32_t k = (phased) ? max_ploidy : 1;  // phased data needs scaled offset
  // a phased sample's offset is the ploidy of every sample before it, and missing is in
  // ascending order, so carry a running total rather than re-summing, which is quadratic
  std::uint32_t summed_to = 0;   // samples already added to haplotypes
  std::uint32_t haplotypes = 0;  // ploidy summed over those samples
  // for samples with missing data, just set values to NA
  for (auto n: missing) {
    offset = max_probs * n;
    if (phased) {
      if (constant_ploidy) {
          offset *= k;
      } else {
        k = ploidy[n];
        for (; summed_to < n; summed_to++) {
          haplotypes += ploidy[summed_to];
        }
        // this matches summing ploidy[i] * max_probs one sample at a time, since
        // multiplying the total instead gives the same value even if it overflows
        offset = haplotypes * max_probs;
      }
    }
    for (std::uint32_t x=0; x<(max_probs * k); x++) {
      probs[offset + x] = std::nan("1");
    }
  }
}

/// parse genotype data for a single variant
///
/// @return 1D float array of genotype probabilties (each from 0.0-1.0).
void Genotypes::probabilities(float * probs) {
  load_data_and_parse_header();
  
  // cleared per decode, so the flag means this decode saw it, rather than some earlier
  // one, even if a decoder is added which does not check
  probs_above_max = false;
  
  std::uint32_t nrows;
  if (!phased) {
    nrows = n_samples;
  } else {
    // phased probabilities require as many rows per sample as the ploidy
    if (constant_ploidy) {
      nrows = n_samples * max_ploidy;
    } else {
      nrows = fast_ploidy_sum(ploidy.get(), n_samples);
    }
  }
  
  if (layout == 1) {
    probabilities_layout1(uncompressed.get(), idx, probs, nrows);
  } else if (layout == 2) {
    probabilities_layout2(uncompressed.get(), idx, probs, nrows); // about 3 milliseconds
  }
}

/// find which allele corresponds to the minor allele
///
/// Rather than checking every individual to see which is the minor allele, we
/// check subsets, in batches of 100. We obtain alleles for individuals in the
/// batch, then check if a confidence interval for the frequency of the less
/// frequent allele could overlap 0.5. If not, we can be reasonably certain the
/// less frequent allele is the true minor allele, without having to check the
/// full cohort. This can be 600X faster than checking the full cohort in larger
/// populations.
///
/// The frequency is divided by the alleles actually summed, since a stride covers a
/// varying number of samples, each with their own ploidy, and missing ones contribute
/// none. A variant sitting at almost exactly 0.5 never settles, so striding stops at
/// stride_limit and the rest runs sequentially, at about half the cost per sample
///
/// @param dose float array of dosages for the reference (first) allele, with
///     missing samples already set to nan
/// @return index for minor allele (0 or 1)
int Genotypes::find_minor_allele(float * dose) {
  std::uint32_t batchsize = 100;
  std::uint32_t increment = std::max(n_samples / batchsize, (std::uint32_t) 1);
  double total = 0;
  double freq = 0;
  // the ploidy summed over the samples visited, and the count of those samples,
  // which is what bounds the interval below
  std::uint64_t n_alleles_seen = 0;
  std::uint64_t n_checked = 0;
  // hoisted out of the loop, so a constant ploidy variant reads max_ploidy rather than an
  // array which would have to be built for the purpose
  bool varies = !constant_ploidy;
  std::uint8_t fixed_ploidy = (std::uint8_t) max_ploidy;
  // striding stops after this many samples, since only a frequency almost exactly at 0.5
  // gets this far, and that never settles, while costing 16 times the cache lines
  std::uint64_t stride_limit = 10000;
  
  // To make sure we don't hit weird groupings of alleles in individuals, this
  // picks samples uniformly thoughout the population, by using an appropriate
  // step size.
  for (std::uint32_t idx2=0; idx2<increment; idx2++) {
    if (varies) {
      for (std::uint32_t n=idx2; n<n_samples; n += increment) {
        // a missing sample has no called alleles, so counts towards neither the total
        // nor the alleles it is divided by, having been marked by get_allele_dosage
        if (std::isnan(dose[n])) {
          continue;
        }
        total += dose[n];
        n_alleles_seen += ploidy[n];
        n_checked += 1;
      }
    } else {
      // every callable sample carries the same ploidy, so the alleles seen follow from
      // the count of them and the stride can be gathered rather than walked
      std::uint64_t before = n_checked;
      total += fast_strided_sum(dose, n_samples, idx2, increment, n_checked);
      n_alleles_seen += (n_checked - before) * (std::uint64_t) fixed_ploidy;
    }
    if (n_alleles_seen == 0) {
      // nothing callable yet, but a later stride may reach samples with alleles
      continue;
    }
    freq = total / (double) n_alleles_seen;
    if (minor_certain(freq, (int) n_checked, 5.0)) {
      return (freq <= 0.5) ? 0 : 1;
    }
    if (n_checked >= stride_limit) {
      // the strides are not going to resolve this one, so sum every sample in order
      // rather than continuing to hop through the array
      break;
    }
  }
  
  // Not resolved by sampling, so use the whole cohort. Summing from scratch keeps the
  // answer independent of how far the strides got
  if (!varies) {
    // every callable sample carries the same ploidy, so the alleles seen follow from the
    // count of them and the sum can be vectorised
    std::uint64_t n_called = 0;
    total = fast_dosage_sum(dose, n_samples, n_called);
    n_alleles_seen = n_called * (std::uint64_t) fixed_ploidy;
  } else {
    total = 0;
    n_alleles_seen = 0;
    for (std::uint32_t n=0; n<n_samples; n++) {
      if (std::isnan(dose[n])) {
        continue;
      }
      total += dose[n];
      n_alleles_seen += ploidy[n];
    }
  }
  if (n_alleles_seen == 0) {
    // no sample has a called allele, so there is no frequency to compare. Report the
    // first allele, as a zero frequency would
    return 0;
  }
  freq = total / (double) n_alleles_seen;
  if (freq <= 0.5) {
    return 0;
  } else {
    return 1;
  }
}

#if defined(__x86_64__)
// the AVX2 half of ref_dosage_fast below, compiled in isolation, so it is only entered
// after a runtime check for AVX2 support
BGEN_TARGET_AVX2
static void ref_dosage_avx2(char * uncompressed, std::uint32_t & idx,
                            float * dose, std::uint32_t & nrows,
                            std::uint32_t & n, bool & above_max) {
  const float c = 1.0f / 255.0f;
  __m256 k = _mm256_set1_ps(c);
  // weights for the dot product below, 2 against each homozygous count and 1 against
  // each heterozygous. Little endian, so the low byte multiplies the even input byte
  __m256i weights = _mm256_set1_epi16(0x0102);
  // the largest count a valid sample reaches, first * 2 + second, so anything above it
  // means the bgen stores probabilities summing above the maximum, and gets clamped
  __m256i limit = _mm256_set1_epi16(510);
  // running maximum count, checked against the limit once the loop ends
  __m256i peak = _mm256_setzero_si256();

  __m256i initial;
  __m256i total;
  __m256i lo;
  __m256i hi;
  for (; n + 16 <= nrows; n+=16) {
    initial = _mm256_loadu_si256((__m256i *) &uncompressed[idx]);

    // multiply each byte by its weight and add the pairs, which applies the ploidy and
    // sums a sample's two counts. The weighted sum needs 10 bits, so 16-bit results
    total = _mm256_maddubs_epi16(initial, weights);

    // hold the largest count seen rather than comparing every iteration, since only
    // whether any sample passed the limit matters, not which
    peak = _mm256_max_epi16(peak, total);

    // clamp as dosage_count8 does, since a malformed sample would give a ref dosage
    // above the ploidy, which the swap to the alt allele turns negative
    total = _mm256_min_epi16(total, limit);

    // convert each half to floats, via 32 bit ints
    lo = _mm256_cvtepi16_epi32(_mm256_castsi256_si128(total));
    hi = _mm256_cvtepi16_epi32(_mm256_extractf128_si256(total, 1));

    _mm256_storeu_ps(&dose[n], _mm256_mul_ps(_mm256_cvtepi32_ps(lo), k));
    _mm256_storeu_ps(&dose[n + 8], _mm256_mul_ps(_mm256_cvtepi32_ps(hi), k));

    idx += 32;
  }
  // one comparison for the whole variant, rather than one per iteration
  if (!_mm256_testz_si256(_mm256_cmpgt_epi16(peak, limit), _mm256_set1_epi8(-1))) {
    above_max = true;
  }
}
#endif

/// calculate dosage of the reference (first) allele for all samples.
///
/// Writes dosage float values to the dose member.
///
/// This is optimised for 8-bit genotypes and constant_ploidy. This function
/// pulls the homozygous minor allele genotype from the first genotype
/// probability. See fast_dosage_minor_second() for getting dosage when the minor
/// allele is the second allele.
///
/// This uses AVX2 and NEON vectorization to speed up calculations on relevant 
/// x86_64 and aarch64 hardware.
///
/// @param uncompressed char array containing genotype probabilities
/// @param idx uint position where the genotype probabilties begin
void Genotypes::ref_dosage_fast(char *uncompressed, std::uint32_t idx, float *dose, std::uint32_t nrows) {
  std::uint32_t n=0;
#if defined(__x86_64__)
  if (__builtin_cpu_supports("avx2")) {
    ref_dosage_avx2(uncompressed, idx, dose, nrows, n, probs_above_max);
  }
#elif defined(__aarch64__)
  // using this optimised method roughly doubles the speed of computing the ref
  // dosage, but it has a limited impact, since 80-90% of the time is spent
  // decompressing the genotypes array (with zlib compressed data).
  std::uint8_t * buff = reinterpret_cast<std::uint8_t *>(uncompressed);
  const float c = 1.0f / 255.0f;
  float32x4_t k = vdupq_n_f32(c);
  uint8x16x2_t initial;
  uint16x8_t total;
  // as in the AVX2 path, the largest count a valid sample reaches, and a running
  // maximum so the check happens once per variant rather than once per iteration
  uint16x8_t limit = vdupq_n_u16(510);
  uint16x8_t peak = vdupq_n_u16(0);
  uint8x8_t two = vdup_n_u8(2);  // ploidy, applied to the homozygous counts
  for (; n + 16 <= nrows; n += 16) {
    // load data from the array into SIMD registers. This deinterleaves the
    // het and hom counts into separate vector registers
    initial = vld2q_u8(buff + idx);

    // widen the low half of the het counts and add twice each hom count in the same
    // instruction, which applies the ploidy and sums a sample's counts, needing 10 bits
    total = vmlal_u8(vmovl_u8(vget_low_u8(initial.val[1])),
                     vget_low_u8(initial.val[0]), two);
    peak = vmaxq_u16(peak, total);
    total = vminq_u16(total, limit);  // clamp as dosage_count8 does, see the AVX2 path

    // still in the low half, expand 16-bit uints to 32-bit, convert to float, store
    vst1q_f32(dose + n, vmulq_f32(vcvtq_f32_u32(vmovl_u16(vget_low_u16(total))), k));
    vst1q_f32(dose + n + 4, vmulq_f32(vcvtq_f32_u32(vmovl_u16(vget_high_u16(total))), k));

    // repeat for the high half of the vectors
    total = vmlal_u8(vmovl_u8(vget_high_u8(initial.val[1])),
                     vget_high_u8(initial.val[0]), two);
    peak = vmaxq_u16(peak, total);
    total = vminq_u16(total, limit);

    vst1q_f32(dose + n + 8, vmulq_f32(vcvtq_f32_u32(vmovl_u16(vget_low_u16(total))), k));
    vst1q_f32(dose + n + 12, vmulq_f32(vcvtq_f32_u32(vmovl_u16(vget_high_u16(total))), k));

    idx += 32;
  }
  if (vmaxvq_u16(peak) > 510) {
    probs_above_max = true;
  }
#endif
  // Finish off the remaining unvectorized samples. This is 50% slower than SIMD.
  // Also handles if we can't use the code above (not aarch64, x86_64 or lacks avx2)
  for (; n < nrows; n++) {
    dose[n] = dosage_count8(*reinterpret_cast<const std::uint8_t *>(&uncompressed[idx]),
                            *reinterpret_cast<const std::uint8_t *>(&uncompressed[idx + 1]),
                            probs_above_max) * (1.0f / 255.0f);
    idx += 2;
  }
}

/// calculate dosage of the reference (first) allele for all samples for unphased genotpes
///
/// The slow path, loops across samples. Figures out ploidy at each step,
/// and the homozygous genotype for the minor allele. Can use any bit depth.
/// For the most part we just have to use the first and second probabilities for
/// a sample, since they correspond to first hom, and het genotypes. This assumes
/// the first allele is for the ref genotype (but this isn't exposed externally,
/// since the calling function always converts to the minor allele dosage).
///
/// Variants encoded with layout 1 need a little bit more work, since we need to
/// skip over the second homozygous probability, and figure out missingness from
/// the probabilities.
///
/// @param uncompressed char array of genotype probabilities (encoding depends on layout)
/// @param idx uint index position in uncompressed where genotype probabilities start
void Genotypes::ref_dosage_slow_unphased(char * uncompressed, std::uint32_t idx, float * dose, std::uint32_t nrows) {
  std::uint32_t curr_ploidy = max_ploidy;
  std::uint32_t half_ploidy = curr_ploidy / 2;

  std::uint32_t maxval = std::pow(2, (std::uint32_t) (bit_depth)) - 1;
  float factor = (layout == 2) ? 1.0f / (float) maxval : 1.0f / 32768;
  // 64-bit, since a probability fills 32 bits at the maximum bit depth, and a diploid
  // doubles it. Signed makes the float conversion one instruction, and these are masked
  std::int64_t het;
  std::int64_t hom;
  std::uint32_t hom_alt;
  std::uint64_t probs_mask = std::uint64_t(0xFFFFFFFFFFFFFFFF) >> (64 - bit_depth);
  std::uint64_t bit_idx = 0;  // index position in bits
  for (std::uint32_t n=0; n<nrows; n++) {
    if (!constant_ploidy) {
      curr_ploidy = this->ploidy[n];
      half_ploidy = curr_ploidy / 2;
    }
    // check the ploidy before reading, since the reads below are sized from it
    if (curr_ploidy > 2) {
      throw std::invalid_argument("cannot compute dosage with ploidy > 2");
    }
    if (curr_ploidy == 0) {
      // a sample with no alleles stores no probabilities at all, so reading one would
      // shift every later sample along, and no alleles means no dosage either
      dose[n] = 0.0f;
      continue;
    }
    hom = ((*reinterpret_cast<const std::uint64_t* >(&uncompressed[idx + bit_idx / 8]) >> bit_idx % 8) & probs_mask);
    bit_idx += bit_depth;
    
    het = 0;
    if (curr_ploidy == 2) {
      het = ((*reinterpret_cast<const std::uint64_t* >(&uncompressed[idx + bit_idx / 8]) >> bit_idx % 8) & probs_mask);
      bit_idx += bit_depth;
    }
    
    dose[n] = ((hom * curr_ploidy) + het * half_ploidy) * factor;
    if ((layout == 2) && ((std::uint64_t) (hom + het) > probs_mask)) {
      // a sample's stored probabilities must sum to at most probs_mask. Layout 1 stores
      // every probability rather than inferring the last, so it has no such constraint
      dose[n] = curr_ploidy;
      probs_above_max = true;
    }
    if (layout == 1) {
      // layout1 stores hom alt probability, and all zeros indicates missingness
      hom_alt = *reinterpret_cast<const std::uint16_t* >(&uncompressed[idx + bit_idx / 8]);
      bit_idx += bit_depth;
      if ((hom == 0) & (het == 0) & (hom_alt == 0)) {
        missing.push_back(n);
      }
    }
  }
}

/// calculate dosage of the reference (first) allele for all samples for phased genotpes
///
/// The slow path, loops across samples. Figures out ploidy at each step,
/// and the homozygous genotype for the reference allele. Can use any bit depth.
/// For the most part we just have to use the first allele for each haplotype for
/// a sample. This assumes the first allele is for the ref genotype.
///
/// @param uncompressed char array of genotype probabilities (encoding depends on layout)
/// @param idx uint index position in uncompressed where genotype probabilities start
void Genotypes::ref_dosage_slow_phased(char * uncompressed, std::uint32_t idx, float * dose, std::uint32_t nrows) {
  std::uint32_t curr_ploidy = max_ploidy;

  std::uint32_t maxval = std::pow(2, (std::uint32_t) (bit_depth)) - 1;
  float factor = (layout == 2) ? 1.0f / (float) maxval : 1.0f / 32768;
  float prob;
  std::uint64_t probs_mask = std::uint64_t(0xFFFFFFFFFFFFFFFF) >> (64 - bit_depth);
  std::uint64_t bit_idx = 0;  // index position in bits
  for (std::uint32_t n=0; n<nrows; n++) {
    if (!constant_ploidy) {
      curr_ploidy = this->ploidy[n];
    }
    
    prob = 0;
    for (std::uint32_t i=0; i<curr_ploidy; i++) {
      prob += ((*reinterpret_cast<const std::uint64_t* >(&uncompressed[idx + bit_idx / 8]) >> bit_idx % 8) & probs_mask) * factor;
      bit_idx += bit_depth;
    }
    dose[n] = prob;
    
  }
}

#if defined(__x86_64__)
// the AVX half of swap_allele_dosage_simple below, which only needs AVX, and is entered
// after a runtime check for it
BGEN_TARGET_AVX
static void swap_dosage_avx(float * dose, std::uint32_t & n_samples,
                            std::uint32_t & n) {
  __m256 k = {2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f};
  __m256 batch;
  for (; n + 8 < n_samples; n+=8) {
    batch = _mm256_loadu_ps(dose + n);
    _mm256_storeu_ps(dose + n, _mm256_sub_ps(k, batch));
  }
}
#endif

/// swap sample dosages to the opposing allele. Requires a biallelic variant.
///
/// This replaces the values in the dose array with 2.0 - value.
///
/// This uses AVX and NEON vectorization to speed up calculations on relevant
/// x86_64 and aarch64 hardware
void Genotypes::swap_allele_dosage_simple(float * dose) {
  std::uint32_t n=0;
#if defined(__x86_64__)
  if (__builtin_cpu_supports("avx")) {
    swap_dosage_avx(dose, n_samples, n);
  }
#elif defined(__aarch64__)
  float32x4_t k = vdupq_n_f32(2.0f);
  float32x4_t batch;
  for (; n + 4 < n_samples; n += 4) {
    batch = vld1q_f32(dose + n);
    vst1q_f32(dose + n, vsubq_f32(k, batch));
  }
#endif
  for (; n + 4 < n_samples; n+=4) {
    dose[n] = 2.0f - dose[n];
    dose[n+1] = 2.0f - dose[n+1];
    dose[n+2] = 2.0f - dose[n+2];
    dose[n+3] = 2.0f - dose[n+3];
  }
  for (; n<n_samples; n++) {
    dose[n] = 2.0f - dose[n];
  }
}

/// swap sample dosages to the opposing allele. Requires a biallelic variant.
///
/// This replaces the values in the dose array with ploidy - value.
void Genotypes::swap_allele_dosage_complex(float * dose) {
  // reached when the ploidy varies or is not 2, so read the array only in the former case
  if (constant_ploidy) {
    float ploid = (float) max_ploidy;
    for (std::uint32_t n=0; n<n_samples; n++) {
      dose[n] = ploid - dose[n];
    }
    return;
  }
  for (std::uint32_t n=0; n<n_samples; n++) {
    dose[n] = (float) (this->ploidy[n]) - dose[n];
  }
}

/// which of a biallelic variant's two alleles is the minor one
///
/// Finding this needs the whole cohort's dosages, which get_allele_dosage works out as a
/// side effect, so they go to a scratch buffer for callers which only want the allele
///
/// @return index into alleles of the minor allele (0 or 1)
int Genotypes::get_minor_idx() {
  if (!minor_known) {
    std::unique_ptr<float[]> dose(new float[n_samples]);
    // this sets minor_idx, and the dosages themselves are not needed
    get_allele_dosage(dose.get(), true, false);
  }
  return minor_idx;
}

/// calculate allele dosage from the genotype probabilities
///
/// This either computes the alternate allele dosage, or the minor allele dosage
/// depending on the use_alt or use_minor arguments.
///
/// @param use_alt whether to return the alt allele dosages
/// @param use_minor whether to return the minor allele dosages
void Genotypes::get_allele_dosage(float * dose, bool use_alt, bool use_minor) {
  if (use_alt == use_minor) {
    throw std::invalid_argument("one of use_alt or use_minor must be true");
  }
  load_data_and_parse_header();
  
  if (n_alleles != 2) {
    throw std::invalid_argument("can't get allele dosages for non-biallelic var.");
  }
  
  // as in probabilities, cleared per decode rather than latched across calls
  probs_above_max = false;
  
  // the layout 1 dosage path appends missing samples to the list, so clear it first, or
  // a second dosage read holds each of them twice
  if (layout == 1) {
    missing.clear();
  }
  
  // calculate the dosage for the first allele for all samples
  if (constant_ploidy & (max_probs == 3) & (bit_depth == 8) & (!phased)) {
    // A fast path when we know the ploidy is constant and the bit depth is 8,
    // this avoids the bit shifts/masks used in the variable bit_depth path.
    ref_dosage_fast(uncompressed.get(), idx, dose, n_samples);
  } else {
    if (!phased) {
      ref_dosage_slow_unphased(uncompressed.get(), idx, dose, n_samples);
    } else {
      ref_dosage_slow_phased(uncompressed.get(), idx, dose, n_samples);
    }
  }
  
  // mark the missing samples before finding the minor allele, since their zeroed dosages
  // would pull the frequency towards zero, and can flip which allele looks minor
  for (auto n: missing) {
    dose[n] = std::nanf("1");
  }
  
  minor_idx = find_minor_allele(dose);
  minor_known = true;
  if (use_alt | (use_minor & (minor_idx != 0))) {
    if (constant_ploidy & (max_ploidy == 2)) {
      swap_allele_dosage_simple(dose);
    } else {
      swap_allele_dosage_complex(dose);
    }
  } 
}

} //namespace bgen
