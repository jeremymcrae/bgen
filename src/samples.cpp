
#include <algorithm>
#include <cstdint>
#include <iterator>
#include <stdexcept>

#include "samples.h"

namespace bgen {

BgenSamples::BgenSamples(std::ifstream & handle, int n_samples) {
  /* initialize sample list if present in the bgen file
  */
  std::uint32_t sample_header_length;
  handle.read(reinterpret_cast<char*>(&sample_header_length), sizeof(std::uint32_t));
  
  std::uint32_t sample_n_check;
  handle.read(reinterpret_cast<char*>(&sample_n_check), sizeof(std::uint32_t));
  if (n_samples != sample_n_check) {
    throw std::invalid_argument("inconsistent number of samples");
  }
  
  std::uint16_t id_len;
  for (int i=0; i<n_samples; i++){
    handle.read(reinterpret_cast<char*>(&id_len), sizeof(id_len));
    std::string sample_id;
    std::copy_n(std::istream_iterator<char>(handle), id_len, std::back_inserter(sample_id));
    samples.push_back(sample_id);
  }
}

std::vector<std::string> split(const std::string &s, char delim) {
  /* split a string by delimiter into a vector of elements
  */
  std::vector<std::string> elems;
  std::istringstream iss(s);
  std::string item;
  while (std::getline(iss, item, delim)) {
    elems.push_back(item);
  }
  return elems;
}

BgenSamples::BgenSamples(std::string path, int n_samples) {
  /* initialize from external sample file
  */
  std::ifstream handle(path, std::ios::in);
  
  std::string line;
  while (std::getline(handle, line)) {
    std::vector<std::string> elems = split(line, '\t');
    auto sample_id = elems[0];
    if ((sample_id == "FID") | (sample_id == "0")) { continue; } // skip header lines
    samples.push_back(sample_id);
  }
  if (n_samples != (int)samples.size()) {
    throw std::invalid_argument("inconsistent number of samples");
  }
}

} // namespace bgen
