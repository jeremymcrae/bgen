
#include <algorithm>
#include <string>

#include "reader.h"

namespace bgen {

CppBgenReader::CppBgenReader(std::string path, std::string sample_path, bool delay_parsing) {
  if (path != "/dev/stdin") {
    handle = std::make_shared<std::ifstream>(path, std::ios::in | std::ios::binary);
  } else {
    is_stdin = true;
    // std::cin is not ours to close, so hold it without owning it
    handle = borrowed_stream(&std::cin);
  }
  if (handle->fail()) {
    throw std::invalid_argument("error reading from '" + path + "'");
  }
  header = Header(handle.get());
  if (header.has_sample_ids) {
    samples = Samples(handle.get(), header.nsamples);
  } else if (sample_path.size() > 0) {
    samples = Samples(sample_path, header.nsamples);
  } else {
    samples = Samples(header.nsamples);
  }
  
  offset = first_variant_offset();
  if (!delay_parsing) {
    parse_all_variants();
  }
}

/// byte position of the first variant in the bgen
///
/// The cast matters, as header.offset is 32-bit, so the addition would wrap
/// before being widened, and a corrupt offset would then point back inside the
/// header rather than past the end of the file.
std::uint64_t CppBgenReader::first_variant_offset() {
  return (std::uint64_t) header.offset + 4;
}

Variant CppBgenReader::next_var() {
  if (handle->eof()) {
    throw std::out_of_range("reached end of file");
  }
  Variant var(handle, offset, header.layout, header.compression, header.nsamples, is_stdin);
  offset = var.next_variant_offset;
  return var;
}

/// load all variants into memory at once
void CppBgenReader::parse_all_variants() {
  if (variants.size() == header.nvariants) {
    return;
  }
  offset = first_variant_offset();
  variants.clear();
  variants.reserve(header.nvariants);
  try {
    for (std::uint32_t idx=0; idx < header.nvariants; idx++) {
      variants.push_back(next_var());
    }
  } catch (...) {
    // Drop the partial list, so a later call retries and raises again, rather
    // than finding a full sized vector and assuming the parse had finished. The
    // variants are pushed on as they parse for the same reason - resizing up
    // front would leave empty variants behind on failure.
    variants.clear();
    throw;
  }
  // finally reset the offset position to the first variant, so we can iterate
  // over variants more easily in python
  offset = first_variant_offset();
}

/// drop a subset of variants passed in by indexes
void CppBgenReader::drop_variants(std::vector<int> indices) {
  // Check every index before dropping any, so a bad index cannot leave the
  // variants half dropped. Without this, an out of range index builds an invalid
  // iterator and erases from outside the vector, which segfaults.
  for (auto idx : indices) {
    if ((idx < 0) || ((std::size_t) idx >= variants.size())) {
      throw std::out_of_range("variant index " + std::to_string(idx) +
                              " is out of range for a bgen with " +
                              std::to_string(variants.size()) + " variants");
    }
  }
  
  // sort indices in descending order, so dropping elemtns doesn't affect later items
  std::sort(indices.rbegin(), indices.rend());
  
  // adjacent_find checks for duplicates without editing indices, whereas
  // std::unique would shuffle them about before we then throw
  if (std::adjacent_find(indices.begin(), indices.end()) != indices.end()) {
    throw std::invalid_argument("can't drop variants with duplicate indices");
  }
  
  for (auto idx : indices) {
    variants.erase(variants.begin() + idx);
  }
  variants.shrink_to_fit();
  
  // and sort the variants again afterward
  std::sort(variants.begin(), variants.end(),
          [] (Variant const& a, Variant const& b) { return a.pos < b.pos; });
}

/// get all the IDs for the variants in the bgen file
std::vector<std::string> CppBgenReader::varids() {
  parse_all_variants();
  std::vector<std::string> varid(variants.size());
  for (std::uint32_t x=0; x<variants.size(); x++) {
    varid[x] = variants[x].varid;
  }
  return varid;
}

/// get all the rsIDs for the variants in the bgen file
std::vector<std::string> CppBgenReader::rsids() {
  parse_all_variants();
  std::vector<std::string> rsid(variants.size());
  for (std::uint32_t x=0; x<variants.size(); x++) {
    rsid[x] = variants[x].rsid;
  }
  return rsid;
}

/// get all the chroms for the variants in the bgen file
std::vector<std::string> CppBgenReader::chroms() {
  parse_all_variants();
  std::vector<std::string> chrom(variants.size());
  for (std::uint32_t x=0; x<variants.size(); x++) {
    chrom[x] = variants[x].chrom;
  }
  return chrom;
}

/// get all the positions for the variants in the bgen file
std::vector<std::uint32_t> CppBgenReader::positions() {
  parse_all_variants();
  std::vector<std::uint32_t> position(variants.size());
  for (std::uint32_t x=0; x<variants.size(); x++) {
    position[x] = variants[x].pos;
  }
  return position;
}

} // namespace bgen
