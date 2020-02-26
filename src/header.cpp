
#include <iterator>
#include <algorithm>
#include <cassert>
#include <bitset>

#include "header.h"

namespace bgen {

unsigned char reverse(unsigned char b) {
  b = (b & 0xF0) >> 4 | (b & 0x0F) << 4;
  b = (b & 0xCC) >> 2 | (b & 0x33) << 2;
  b = (b & 0xAA) >> 1 | (b & 0x55) << 1;
  return b;
}

std::bitset<32> read_flags(std::ifstream & handle) {
  std::bitset<32> flags;
  for (int i=0; i < 4; i++) {
    char a;
    handle.read(&a, 1);
    std::bitset<8> tmp(reverse(a));
    for (int j = 0; j < 8; j++){
        flags[(i * 8) + j] = tmp[j];
      }
  }
  return flags;
}

Header::Header(std::ifstream & handle) {
  handle.read(reinterpret_cast<char*>(&offset), sizeof(std::uint32_t));
  handle.read(reinterpret_cast<char*>(&header_length), sizeof(std::uint32_t));
  handle.read(reinterpret_cast<char*>(&nvariants), sizeof(std::uint32_t));
  handle.read(reinterpret_cast<char*>(&nsamples), sizeof(std::uint32_t));
  std::copy_n(std::istream_iterator<char>(handle), sizeof(std::uint32_t), std::back_inserter(magic));
  
  // make sure we are reading a bgen file
  assert (magic == "bgen" | magic == "0000");
  
  // read any extra data contained in the header
  int size = header_length - 20;
  if (size > 0) {
    std::copy_n(std::istream_iterator<char>(handle), size, std::back_inserter(extra));
  }
  
  // read flags data
  std::bitset<32> flags = read_flags(handle);
  
  // compression, where 0=no_compression, 1=zlib, 2=zstd
  compression = flags[0] + flags[1];
  layout = flags[2] + flags[3] + flags[4] + flags[5];
  has_sample_ids = flags[31];
}

} // namespace bgen
