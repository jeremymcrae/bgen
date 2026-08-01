
#include <bitset>
#include <stdexcept>

#include "header.h"

namespace bgen {

Header::Header(std::istream * handle) {
  if (!read_value(*handle, offset) ||
      !read_value(*handle, header_length) ||
      !read_value(*handle, nvariants) ||
      !read_value(*handle, nsamples)) {
    throw std::invalid_argument("doesn't appear to be a bgen file");
  }
  magic.resize(4);
  if (!handle->read(&magic[0], 4)) {
    throw std::invalid_argument("doesn't appear to be a bgen file");
  }
  
  // the magic bytes are either 'bgen' or, for some files, four zero bytes
  if ((magic != "bgen") && (
      (int) magic[0] != 0 || 
      (int) magic[1] != 0 || 
      (int) magic[2] != 0 || 
      (int) magic[3] != 0)) {
    throw std::invalid_argument("doesn't appear to be a bgen file");
  }
  
  // The header runs from byte 4 to 4 + header_length, and the variant data
  // starts at 4 + offset, so the header has to fit within the offset. Checking
  // that stops a corrupt header_length from wrapping the unsigned subtraction
  // below, or sizing a huge allocation from a 4 byte field.
  if ((header_length < 20) || (header_length > offset)) {
    throw std::invalid_argument("bgen header is the wrong length");
  }
  
  // read any extra data contained in the header
  std::uint32_t size = header_length - 20;
  if (size > 0) {
    extra.resize(size);
    if (!handle->read(&extra[0], size)) {
      throw std::invalid_argument("bgen file is truncated inside the header");
    }
  }
  
  // read flags data
  std::uint32_t raw_flags;
  if (!read_value(*handle, raw_flags)) {
    throw std::invalid_argument("bgen file is truncated inside the header");
  }
  std::bitset<32> flags(raw_flags);
  
  std::bitset<32> compr_mask(0b11);
  std::bitset<32> layout_mask(0b111100);
  compression = (int) (flags & compr_mask).to_ulong();
  layout = (int) ((flags & layout_mask) >> 2).to_ulong();
  has_sample_ids = flags[31];
  
  // both fields have spare bits, so can hold values the spec does not define.
  // Reject those, as the genotype parsing has no branch for them, and would
  // otherwise return an unwritten buffer.
  if ((layout != 1) && (layout != 2)) {
    throw std::invalid_argument("unknown bgen layout version");
  }
  if (compression > 2) {
    throw std::invalid_argument("unknown bgen compression scheme");
  }
}

} // namespace bgen
