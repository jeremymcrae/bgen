
#include <cstring>
#include <stdexcept>

#include "header.h"
#include "utils.h"

namespace bgen {

Header::Header(std::istream * handle) {
  std::uint32_t header_length;
  char magic[4];
  if (!read_value(*handle, offset) ||
      !read_value(*handle, header_length) ||
      !read_value(*handle, nvariants) ||
      !read_value(*handle, nsamples) ||
      !handle->read(magic, 4)) {
    throw std::invalid_argument("doesn't appear to be a bgen file");
  }
  
  // the magic bytes are either 'bgen' or, for some files, four zero bytes
  if ((std::memcmp(magic, "bgen", 4) != 0) && (std::memcmp(magic, "\0\0\0\0", 4) != 0)) {
    throw std::invalid_argument("doesn't appear to be a bgen file");
  }
  
  // the header has to fit within the offset where the variants start, which stops a
  // corrupt header_length wrapping the subtraction below, or sizing a huge allocation
  if ((header_length < 20) || (header_length > offset)) {
    throw std::invalid_argument("bgen header is the wrong length");
  }
  
  // read any extra data contained in the header
  std::uint32_t size = header_length - 20;
  if (size > 0) {
    if (!read_string(*handle, extra, size)) {
      throw std::invalid_argument("bgen file is truncated inside the header");
    }
  }
  
  // read flags data
  std::uint32_t raw_flags;
  if (!read_value(*handle, raw_flags)) {
    throw std::invalid_argument("bgen file is truncated inside the header");
  }
  compression = raw_flags & 3;
  layout = (raw_flags >> 2) & 15;
  has_sample_ids = (raw_flags >> 31) & 1;
  
  // both fields have spare bits, so reject values the spec does not define
  if ((layout != 1) && (layout != 2)) {
    throw std::invalid_argument("unknown bgen layout version");
  }
  if (compression > 2) {
    throw std::invalid_argument("unknown bgen compression scheme");
  }
}

} // namespace bgen
