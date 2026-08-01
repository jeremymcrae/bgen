#ifndef BGEN_SAMPLES_H_
#define BGEN_SAMPLES_H_

#include <cstdint>
#include <fstream>
#include <vector>
#include <string>

namespace bgen {

class Samples {
public:
  // the sample count is unsigned to match the bgen header field. Taking it as an
  // int would turn counts above INT_MAX negative, which then becomes a huge
  // size_t when it reaches resize(), and silently skips the loops below.
  Samples(std::istream * handle, std::uint32_t n_samples);
  Samples(std::string path, std::uint32_t n_samples);
  Samples(std::uint32_t n_samples);
  Samples() {}
  std::vector<std::string> samples;
};

} // namespace bgen

#endif  // BGEN_SAMPLES_H_
