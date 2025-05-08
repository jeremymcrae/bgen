#ifndef BGEN_HEADER_H_
#define BGEN_HEADER_H_

#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <cstdint>

#include "utils.h"

namespace bgen {

class Header {
  std::uint32_t header_length;
  std::string magic;
public:
  Header(std::istream * handle);
  Header() {}
  std::uint32_t offset;
  std::uint32_t nvariants;
  std::uint32_t nsamples = 5;
  int compression;
  int layout;
  bool has_sample_ids;
  std::string extra;
};

} // namespace bgen

#endif  // BGEN_HEADER_H_
