
#include <algorithm>
#include <cstdint>
#include <fstream>
#include <stdexcept>

#include "samples.h"
#include "utils.h"

namespace bgen {

/// most IDs to reserve space for up front, since a corrupt count reserves tens of bytes
/// per claimed sample. Growing past this only costs bgens over a million samples
const std::uint64_t MAX_ID_RESERVE = 1 << 20;

/// how many samples a block of this many bytes could describe. Each ID has a two byte
/// length prefix, after the two four byte counts, so this bound is generous but sound
static std::uint64_t samples_in_block(std::uint64_t block_length) {
  if (block_length <= 8) {
    return 0;
  }
  return (block_length - 8) / 2;
}

/// initialize sample list if present in the bgen file
///
///  @param handle stream positioned at the start of the sample block
///  @param n_samples sample count from the header, not yet trusted
///  @param file_size size of the bgen, or zero when it is not known
Samples::Samples(std::istream * handle, std::uint32_t n_samples,
                 std::uint64_t file_size) {
  // the sample block length bounds the sample count, and is checked against the IDs
  // at the end of this function
  std::uint32_t sample_header_length;
  std::uint32_t sample_n_check;
  if (!read_value(*handle, sample_header_length) || !read_value(*handle, sample_n_check)) {
    throw std::invalid_argument("bgen file is truncated inside the sample block");
  }
  if (n_samples != sample_n_check) {
    throw std::invalid_argument("inconsistent number of samples");
  }
  
  std::uint64_t capacity = samples_in_block(sample_header_length);
  if (file_size > 0) {
    capacity = std::min(capacity, samples_in_block(file_size));
  }
  if (n_samples > capacity) {
    // The data available in the sample block is not enough for the number of samples
    // declared in the header.
    throw std::invalid_argument("bgen has more samples than the sample block can hold");
  }

  // count the bytes the IDs occupy, to compare with sample_header_length afterwards
  std::uint64_t used = 8;
  
  // an ID costs far more in memory than the two bytes which justified counting it, so
  // reserve a bounded amount and grow as they arrive
  samples.reserve(std::min((std::uint64_t) n_samples, MAX_ID_RESERVE));
  for (std::uint32_t i=0; i<n_samples; i++) {
    samples.emplace_back();
    // read the IDs as raw bytes, so that IDs containing spaces survive intact
    if (!read_prefixed_string<std::uint16_t>(*handle, samples.back())) {
      throw std::invalid_argument("bgen file is truncated inside the sample block");
    }
    used += 2 + samples.back().size();
  }
  
  // The IDs should have taken the block right up to its declared length.
  if (used != sample_header_length) {
    throw std::invalid_argument("bgen sample block length does not match its sample IDs");
  }
}

/// initialize from external sample file
Samples::Samples(std::string path, std::uint32_t n_samples) {
  // binary, so a windows line ending stays on the line, and is dropped with the other
  // columns below
  std::ifstream handle(path, std::ios::in | std::ios::binary);
  if (!handle) {
    throw std::invalid_argument("error with sample file: '" + path + "'");
  }
  
  // skip the two header lines
  std::string line;
  std::getline(handle, line, '\n');
  std::getline(handle, line, '\n');
  
  // as above, the count cannot be trusted to size the allocation
  samples.reserve(std::min((std::uint64_t) n_samples, MAX_ID_RESERVE));
  
  // run through all lines and get the first column as sample_id
  while (std::getline(handle, line, '\n')) {
    // the .sample format is whitespace delimited, so take the first column, which also
    // drops any carriage return a windows file left on the end of the line
    std::string id = line.substr(0, line.find_first_of(" \t\r"));
    
    // skip empty lines, which includes the lone \r of a blank windows line, and lines
    // starting with a null character
    if ((id.size() == 0) || (id[0] == 0)) {
      continue;
    }
    
    if (samples.size() >= n_samples) {
      throw std::invalid_argument("inconsistent number of samples");
    }
    samples.push_back(std::move(id));
  }
  if (handle.bad()) {
    throw std::invalid_argument("error reading sample file: '" + path + "'");
  }
  
  if (samples.size() != n_samples) {
    throw std::invalid_argument("inconsistent number of samples");
  }
}

/// initialize with integer IDs if no sample list available
///
/// The IDs wait for get_samples(), which checks the count against the file size, since
/// nothing in the bgen can, and a corrupt count would otherwise allocate arbitrarily
Samples::Samples(std::uint32_t n_samples, std::uint64_t _file_size,
                 std::uint32_t _n_variants, int _compression) {
  n_placeholders = n_samples;
  file_size = _file_size;
  n_variants = _n_variants;
  compression = _compression;
}

/// sample IDs, generating placeholders on first use if the bgen had none. The count is
/// checked against the file size, to keep a small bgen from building billions of them
const std::vector<std::string> & Samples::get_samples() {
  if (n_placeholders > 0) {
    if (file_size > 0) {
      // 0-variant files hold only the header, while compressed variants can expand by
      // up to MAX_EXPANSION. Cap against impossible allocations while allowing real cohorts
      bool cannot_describe = false;
      if (n_variants == 0) {
        cannot_describe = (n_placeholders > 10000000);
      } else if (compression == 0) {
        cannot_describe = ((std::uint64_t) n_placeholders > file_size);
      } else {
        cannot_describe = ((std::uint64_t) n_placeholders > file_size * MAX_EXPANSION);
      }
      if (cannot_describe) {
        throw std::invalid_argument("bgen claims " + std::to_string(n_placeholders) +
                                    " samples, which a file of " +
                                    std::to_string(file_size) +
                                    " bytes cannot describe");
      }
    }
    samples.reserve(std::min((std::uint64_t) n_placeholders, MAX_ID_RESERVE));
    for (std::uint32_t i=0; i<n_placeholders; i++) {
      samples.push_back(std::to_string(i));
    }
    n_placeholders = 0;
  }
  return samples;
}

} // namespace bgen
