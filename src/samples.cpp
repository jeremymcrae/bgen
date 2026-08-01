
#include <cstdint>
#include <stdexcept>

#include "samples.h"
#include "utils.h"

namespace bgen {

Samples::Samples(std::istream * handle, std::uint32_t n_samples) {
  /* initialize sample list if present in the bgen file
  */
  // the sample block length is read only to step over it, the IDs follow it
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
  
  samples.resize(n_samples);
  for (std::uint32_t i=0; i<n_samples; i++) {
    // read the IDs as raw bytes, so that IDs containing spaces survive intact
    if (!read_prefixed_string<std::uint16_t>(*handle, samples[i])) {
      throw std::invalid_argument("bgen file is truncated inside the sample block");
    }
  }
}

Samples::Samples(std::string path, std::uint32_t n_samples) {
  /* initialize from external sample file
  */
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
  auto pos = handle.tellg();
  handle.seekg(0, std::ios::end);
  auto fsize = (std::uint64_t) handle.tellg() - pos;
  std::string lines(fsize, '\0');
  handle.seekg(pos);
  handle.read(&lines[0], fsize);
  
  samples.resize(n_samples);
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
    samples[idx] = line.substr(0, line.find(' '));
    idx += 1;
  }
  
  if (idx != n_samples) {
    throw std::invalid_argument("inconsistent number of samples");
  }
  
  if (n_samples != samples.size()) {
    throw std::invalid_argument("inconsistent number of samples");
  }
}

Samples::Samples(std::uint32_t n_samples) {
  /* initialize with integer IDs if no sample list available
  */
  samples.resize(n_samples);
  for (std::uint32_t i=0; i<n_samples; i++) {
    samples[i] = std::to_string(i);
  }
}

} // namespace bgen
