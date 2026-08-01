#ifndef BGEN_HEADER_H_
#define BGEN_HEADER_H_

#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <cstdint>

#include "utils.h"

namespace bgen {

class Header {
  std::uint32_t header_length = 0;
  std::string magic;
public:
  Header(std::istream * handle);
  Header() {}
  std::uint32_t offset = 0;
  std::uint32_t nvariants = 0;
  std::uint32_t nsamples = 0;
  int compression = 0;
  int layout = 0;
  bool has_sample_ids = false;
  std::string extra;
};

} // namespace bgen

#endif  // BGEN_HEADER_H_
