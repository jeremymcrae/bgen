#ifndef BGEN_READER_H_
#define BGEN_READER_H_

#include <fstream>
#include <memory>
#include <stdexcept>
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
  Variant next_var();
  void drop_variants(std::vector<int> indices);
  // the bgen stream, shared with every Variant opened from this reader, so the
  // file is closed once this reader and all of its variants are gone
  std::shared_ptr<std::istream> handle;
  std::vector<std::string> varids();
  std::vector<std::string> rsids();
  std::vector<std::string> chroms();
  std::vector<std::uint32_t> positions();
  Variant & operator[](std::size_t idx) { return variants[idx]; }
  Variant & get(std::size_t idx) { return variants[idx]; }
  std::vector<Variant> variants;
  Header header;
  Samples samples;
  std::uint64_t offset;
};

} // namespace bgen

#endif  // BGEN_READER_H_
