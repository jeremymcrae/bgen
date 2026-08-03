
#include <algorithm>
#include <cstdint>
#include <stdexcept>

#include "samples.h"
#include "utils.h"

namespace bgen {

/// most IDs to reserve space for up front when nothing bounds the count
///
/// Reserving the full count would ask the allocator for tens of bytes per claimed
/// sample, which a corrupt count can push past any address space limit even though
/// the pages are never touched. This is only needed where the bgen's own size is
/// unknown, so growing past it costs a few reallocations on a stream.
const std::uint64_t MAX_ID_RESERVE = 1 << 20;

/// how many samples a block of this many bytes could describe
///
/// Each ID carries a two byte length prefix, and the block also holds the two four
/// byte counts at its start, so this is a generous bound but a sound one.
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
  if (!read_value(*handle, sample_header_length)) {
    throw std::invalid_argument("bgen file is truncated inside the sample block");
  }
  
  std::uint32_t sample_n_check;
  if (!read_value(*handle, sample_n_check)) {
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
  
  if (file_size > 0) {
    samples.resize(n_samples);
    for (std::uint32_t i=0; i<n_samples; i++) {
      // read the IDs as raw bytes, so that IDs containing spaces survive intact
      if (!read_prefixed_string<std::uint16_t>(*handle, samples[i])) {
        throw std::invalid_argument("bgen file is truncated inside the sample block");
      }
      used += 2 + samples[i].size();
    }
  } else {
    // On a stream we cannot use the block length to validate sample count. Reserving
    // only takes address space, but could hit the process memory limit with large values,
    // so cap it and grow as the IDs actually arrive.
    samples.reserve(std::min((std::uint64_t) n_samples, MAX_ID_RESERVE));
    for (std::uint32_t i=0; i<n_samples; i++) {
      samples.emplace_back();
      if (!read_prefixed_string<std::uint16_t>(*handle, samples.back())) {
        throw std::invalid_argument("bgen file is truncated inside the sample block");
      }
      used += 2 + samples.back().size();
    }
  }
  
  // The IDs should have taken the block right up to its declared length.
  if (used != sample_header_length) {
    throw std::invalid_argument("bgen sample block length does not match its sample IDs");
  }
}

/// initialize from external sample file
Samples::Samples(std::string path, std::uint32_t n_samples) {
  std::ifstream handle(path, std::ios::in);
  if (!handle) {
    throw std::invalid_argument("error with sample file: '" + path + "'");
  }
  
  // read first two header lines
  std::string header;
  std::getline(handle, header, '\n');
  std::string types;
  std::getline(handle, types, '\n');
  
  // find the file length post header, then read it all in to memory
  std::string lines;
  if (handle) {
    // Only measure and read when the header lines were both there. On a file with
    // fewer than two lines the stream has already failed, and tellg() then
    // returns -1, so the size below would be meaningless. Leaving lines empty
    // lets the sample count check report the problem.
    auto pos = handle.tellg();
    handle.seekg(0, std::ios::end);
    auto end_pos = handle.tellg();
    if ((pos >= 0) && (end_pos > pos)) {
      std::uint64_t fsize = (std::uint64_t) (end_pos - pos);
      lines.resize(fsize);
      handle.seekg(pos);
      if (!handle.read(&lines[0], fsize)) {
        throw std::invalid_argument("error reading sample file: '" + path + "'");
      }
    }
  }
  
  // As above, the count cannot be trusted to size the allocation. Every ID needs at
  // least a character and a line ending, so the text read from the sample file
  // bounds how many it can hold
  samples.reserve(std::min((std::uint64_t) n_samples, lines.size() / 2));
  std::istringstream iss(lines);
  
  // run through all lines and gte the first column as sample_id
  std::uint32_t idx = 0;
  std::string line;
  while (std::getline(iss, line, '\n')) {
    // skip empty lines
    if ((line.size() == 0) || (line[0] == 0)) {
      // std::getline() on win32 at end of file can create string with null
      // characters. The null character indicates the line doesn't contain an
      // ID, even though string size might be > 0. Only affects win32.
      continue;
    }
    
    if (idx >= n_samples) {
      throw std::invalid_argument("inconsistent number of samples");
    }
    // the .sample format is whitespace delimited, so take the first column, and
    // drop any carriage return left on the end of a line by a windows file
    samples.push_back(line.substr(0, line.find_first_of(" \t\r")));
    idx += 1;
  }
  
  if (idx != n_samples) {
    throw std::invalid_argument("inconsistent number of samples");
  }
}

/// initialize with integer IDs if no sample list available
///
/// The IDs are not built here. Nothing has checked n_samples at this point, and
/// no per sample data exists in the bgen to check it against, so building them
/// now would let a corrupt count allocate arbitrarily. Reading a variant
/// validates the count against the genotype data, so the IDs are left until
/// something asks for them.
Samples::Samples(std::uint32_t n_samples) {
  n_placeholders = n_samples;
}

/// sample IDs, generating placeholders on first use if the bgen had none
const std::vector<std::string> & Samples::get_samples() {
  if (n_placeholders > 0) {
    samples.reserve(std::min((std::uint64_t) n_placeholders, MAX_ID_RESERVE));
    for (std::uint32_t i=0; i<n_placeholders; i++) {
      samples.push_back(std::to_string(i));
    }
    n_placeholders = 0;
  }
  return samples;
}

} // namespace bgen
