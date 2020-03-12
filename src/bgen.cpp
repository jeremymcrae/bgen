
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
  
  // seek forward to the first variant
  handle.seekg(header.offset + 4);
  
  while (true) {
    if (handle.eof() | ((std::uint64_t) handle.tellg() >= fsize)) {
      break;
    }
    Variant variant(handle, header.layout, header.compression, header.nsamples);
    variants.push_back(variant);
  }
}

void Bgen::drop_variants(std::vector<int> indices) {
  /* drop a subset of variants passed in by indexes
  */
  // sort indices in descending order, so dropping elemtns doesn't affect later items
  std::sort(indices.rbegin(), indices.rend());
  
  auto it = std::unique(indices.begin(), indices.end());
  if (it != indices.end()) {
    throw std::invalid_argument("can't drop variants with duplicate indices");
  }
  
  std::vector<Variant>::iterator begin = variants.begin();
  for (auto idx : indices) {
    variants.erase(begin + idx);
  }
}

std::vector<std::string> Bgen::varids() {
  /* get all the IDs for the variants in the bgen file
  */
  std::vector<std::string> varid(variants.size());
  for (uint x=0; x<variants.size(); x++) {
    varid[x] = variants[x].varid;
  }
  return varid;
}

std::vector<std::string> Bgen::rsids() {
  /* get all the rsIDs for the variants in the bgen file
  */
  std::vector<std::string> rsid(variants.size());
  for (uint x=0; x<variants.size(); x++) {
    rsid[x] = variants[x].rsid;
  }
  return rsid;
}

std::vector<std::string> Bgen::chroms() {
  /* get all the chroms for the variants in the bgen file
  */
  std::vector<std::string> chrom(variants.size());
  for (uint x=0; x<variants.size(); x++) {
    chrom[x] = variants[x].chrom;
  }
  return chrom;
}

std::vector<std::uint32_t> Bgen::positions() {
  /* get all the positions for the variants in the bgen file
  */
  std::vector<std::uint32_t> position(variants.size());
  for (uint x=0; x<variants.size(); x++) {
    position[x] = variants[x].pos;
  }
  return position;
}

} // namespace bgen
