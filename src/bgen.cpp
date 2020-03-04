
#include "bgen.h"

namespace bgen {

Bgen::Bgen(std::string path, std::string sample_path) {
  handle.open(path, std::ios::binary);
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
    variants.push_back(variant);
  }
}

std::vector<std::string> Bgen::rsids() {
  /* get all the rsIDs for the variants in the bgen file
  */
  std::vector<std::string> rsid(variants.size());
  for (int x=0; x<variants.size(); x++) {
    rsid[x] = variants[x].rsid;
  }
  return rsid;
}

std::vector<std::string> Bgen::chroms() {
  /* get all the rsIDs for the variants in the bgen file
  */
  std::vector<std::string> chrom(variants.size());
  for (int x=0; x<variants.size(); x++) {
    chrom[x] = variants[x].chrom;
  }
  return chrom;
}

std::vector<std::uint32_t> Bgen::positions() {
  /* get all the rsIDs for the variants in the bgen file
  */
  std::vector<std::uint32_t> position(variants.size());
  for (int x=0; x<variants.size(); x++) {
    position[x] = variants[x].pos;
  }
  return position;
}

} // namespace bgen
