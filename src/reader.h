#ifndef BGEN_READER_H_
#define BGEN_READER_H_

#include <cstdint>
#include <istream>
#include <memory>
#include <string>
#include <vector>

#include "header.h"
#include "samples.h"
#include "variant.h"

namespace bgen {

class CppBgenReader {
  bool is_stdin = false;
  // size of the bgen in bytes, or zero if not known (stdin cannot be seeked)
  std::uint64_t file_size = 0;
public:
  CppBgenReader(std::string path, std::string sample_path = "", bool delay_parsing = false);
  void close_stream();
  std::uint64_t first_variant_offset();
  void parse_all_variants();
  void drop_variants(std::vector<int> indices);
  // the bgen stream, shared with every Variant opened from this reader, so the
  // file is closed once this reader and all of its variants are gone
  std::shared_ptr<std::istream> handle;
  std::vector<Variant> variants;
  Header header;
  Samples samples;
};

} // namespace bgen

#endif  // BGEN_READER_H_
