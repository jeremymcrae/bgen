
#include <algorithm>
#include <cstdint>
#include <iterator>
#include <stdexcept>

#include "samples.h"

namespace bgen {

BgenSamples::BgenSamples(std::ifstream & handle, int n_samples) {
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

} // namespace bgen
