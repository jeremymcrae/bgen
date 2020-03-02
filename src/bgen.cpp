
#include "bgen.h"

namespace bgen {

Bgen::Bgen(std::string path, std::string sample_path) {
  std::ifstream handle(path, std::ios::binary);
  std::uint64_t fsize = handle.tellg();
  header = Header(handle);
  if (header.has_sample_ids) {
    samples = Samples(handle, header.nsamples);
  } else if (sample_path.size() > 0) {
    samples = Samples(sample_path, header.nsamples);
  } else {
    samples = Samples(header.nsamples);
  }
  
  // figure out the file length, so we don't go beyond it
  std::uint64_t current = handle.tellg();
  handle.seekg(0, std::ios::end);
  fsize = (std::uint64_t) handle.tellg() - fsize;
  handle.seekg(current);
  
  while (true) {
    if (handle.eof() | (std::uint64_t) handle.tellg() >= fsize) {
      break;
    }
    Variant variant(handle, header.layout, header.compression, header.nsamples);
    handle.seekg(variant.next_variant_offset());
    variants.push_back(variant);
  }
}

} // namespace bgen
