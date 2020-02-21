
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <cstdint>

class BgenHeader {
  std::uint32_t offset;
  std::uint32_t header_length;
  std::string magic;
  std::string extra;
public:
  BgenHeader(std::ifstream & handle);
  BgenHeader() {};
  std::uint32_t nvariants;
  std::uint32_t nsamples = 5;
  std::string compression;
  int layout;
  bool has_sample_ids;
};
