
#include <iterator>
#include <algorithm>
#include <cassert>
#include <bitset>

#include "header.h"

namespace bgen {

Header::Header(std::ifstream & handle) {
  handle.read(reinterpret_cast<char*>(&offset), sizeof(std::uint32_t));
  handle.read(reinterpret_cast<char*>(&header_length), sizeof(std::uint32_t));
  handle.read(reinterpret_cast<char*>(&nvariants), sizeof(std::uint32_t));
  handle.read(reinterpret_cast<char*>(&nsamples), sizeof(std::uint32_t));
  std::copy_n(std::istream_iterator<char>(handle), sizeof(std::uint32_t), std::back_inserter(magic));
  
  // make sure we are reading a bgen file
  if ((magic != "bgen") & (magic != "0000")) {
    throw std::invalid_argument("doesn't appear to be a bgen file");
  }
  
  // read any extra data contained in the header
  int size = header_length - 20;
  if (size > 0) {
    std::copy_n(std::istream_iterator<char>(handle), size, std::back_inserter(extra));
  }
  
  // read flags data
  std::bitset<32> flags;
  handle.read(reinterpret_cast<char*>(&flags), sizeof(std::uint32_t));
  
  std::bitset<32> compr_mask(0b000000000000000000000000000011);
  std::bitset<32> layout_mask(0b000000000000000000000000111100);
  compression = (int) (flags & compr_mask).to_ulong();
  layout = (int) ((flags & layout_mask) >> 2).to_ulong();
  has_sample_ids = flags[31];
}

} // namespace bgen
