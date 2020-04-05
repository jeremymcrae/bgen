
#include <algorithm>
#include <cstdint>
#include <iterator>
#include <stdexcept>

#include "samples.h"
#include "utils.h"

namespace bgen {

Samples::Samples(std::ifstream & handle, int n_samples) {
  /* initialize sample list if present in the bgen file
  */
  std::uint32_t sample_header_length;
  handle.read(reinterpret_cast<char*>(&sample_header_length), sizeof(std::uint32_t));
  
  std::uint32_t sample_n_check;
  handle.read(reinterpret_cast<char*>(&sample_n_check), sizeof(std::uint32_t));
  if (n_samples != (int) sample_n_check) {
    throw std::invalid_argument("inconsistent number of samples");
  }
  
  std::uint16_t id_len;
  for (int i=0; i<n_samples; i++) {
    handle.read(reinterpret_cast<char*>(&id_len), sizeof(id_len));
    std::string sample_id;
    std::copy_n(std::istream_iterator<char>(handle), id_len, std::back_inserter(sample_id));
    samples.push_back(sample_id);
  }
}

Samples::Samples(std::string path, int n_samples) {
  /* initialize from external sample file
  */
  std::ifstream handle(path, std::ios::in);
  if (!handle) {
    throw std::invalid_argument("error with sample file: '" + path + "'");
  }
  
  std::string header;
  std::getline(handle, header);
  std::string types;
  std::getline(handle, types);
  
  std::string line;
  while (std::getline(handle, line)) {
    std::vector<std::string> elems = split(line, ' ');
    auto sample_id = elems[0];
    samples.push_back(sample_id);
  }
  if (n_samples != (int)samples.size()) {
    throw std::invalid_argument("inconsistent number of samples");
  }
}

Samples::Samples(int n_samples) {
  /* initialize with integer IDs if no sample list available
  */
  for (int i=0; i<n_samples; i++) {
    samples.push_back(std::to_string(i));
  }
}

} // namespace bgen
