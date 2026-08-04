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
  // file_size bounds the sample block, and is zero when the size is unknown
  Samples(std::istream * handle, std::uint32_t n_samples,
          std::uint64_t file_size=0);
  Samples(std::string path, std::uint32_t n_samples);
  // a bgen with no IDs in it. file_size limits how many samples the file could hold
  Samples(std::uint32_t n_samples, std::uint64_t file_size=0);
  Samples() {}
  const std::vector<std::string> & get_samples();
private:
  std::vector<std::string> samples;
  // number of placeholder IDs left to build, which is only non-zero for a bgen
  // with no IDs of its own. Those are generated on request rather than up front,
  // so the sample count cannot size an allocation before anything has checked it
  std::uint32_t n_placeholders = 0;
  // bytes the bgen holds, used to reject a placeholder count the file cannot
  // back. Zero means the size is unknown, so no bound can be applied
  std::uint64_t file_size = 0;
};

} // namespace bgen

#endif  // BGEN_SAMPLES_H_
