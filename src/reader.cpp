
#include <algorithm>
#include <fstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "reader.h"
#include "streaming.h"

namespace bgen {

/// smallest number of bytes a variant can occupy in a bgen. This underestimates the
/// per variant fields, which only makes the reserve it caps smaller than it could be
const std::uint64_t MIN_VARIANT_BYTES = 12;

/// most variants to reserve space up front for a stream, whose size cannot be
/// measured. A Variant is a few hundred bytes, so a corrupt count asks for terabytes
const std::uint64_t MAX_VARIANT_RESERVE = 1 << 16;

CppBgenReader::CppBgenReader(const std::string &path, const std::string &sample_path,
                             bool delay_parsing) {
  is_stdin = (path == "/dev/stdin");
  // stdin is read from its descriptor rather than held as std::cin, so that the error
  // state and the buffered read-ahead belong to this reader alone
  if (is_stdin) {
    handle = std::make_shared<DescriptorStream>(0);
  } else {
    handle = std::make_shared<BufferedFile>(path);
  }
  if (handle->fail()) {
    throw std::invalid_argument("error reading from '" + path + "'");
  }
  if (!is_stdin) {
    // Record the file size, so counts read from the header can be checked
    // against what the file could actually hold. stdin cannot be seeked, so its
    // size is left as zero, meaning unknown.
    handle->seekg(0, std::ios::end);
    if (handle->good()) {
      file_size = (std::uint64_t) handle->tellg();
    }
    handle->clear();
    handle->seekg(0);
  }
  header = Header(handle.get());
  if (header.has_sample_ids) {
    // Reject a wildly wrong count before the sample block is read, which keeps this
    // coarse failure separate from a block that merely disagrees with its IDs
    if ((file_size > 0) && (header.nsamples > file_size / 2)) {
      throw std::invalid_argument("bgen has more samples than the file can hold");
    }
    samples = Samples(handle.get(), header.nsamples, file_size);
  } else if (sample_path.size() > 0) {
    // the IDs come from outside the bgen, so the bgen cannot bound the count.
    // Samples() grows the list as it reads instead, and rejects a mismatch
    samples = Samples(sample_path, header.nsamples);
  } else {
    // no IDs anywhere, so nothing can bound the count here. Samples() defers
    // building the placeholder IDs until they are asked for, and checks the count
    // against the file size at that point
    samples = Samples(header.nsamples, file_size, header.nvariants, header.compression);
  }
  
  if (!delay_parsing) {
    parse_all_variants();
  }
}

/// release the operating system file handle, which windows needs before the file can
/// be deleted. Variants share the stream, so dropping our reference cannot do it
void CppBgenReader::close_stream() {
  DescriptorStream * stream = dynamic_cast<DescriptorStream *>(handle.get());
  if (stream != nullptr) {
    stream->close();
    return;
  }
  std::ifstream * file = dynamic_cast<std::ifstream *>(handle.get());
  if ((file != nullptr) && file->is_open()) {
    file->close();
  }
}

/// byte position of the first variant in the bgen, in 64 bits so that a corrupt
/// header.offset cannot wrap
std::uint64_t CppBgenReader::first_variant_offset() {
  return (std::uint64_t) header.offset + 4;
}

/// load all variants into memory at once
void CppBgenReader::parse_all_variants() {
  if (variants.size() == header.nvariants) {
    return;
  }
  variants.clear();
  // only reserve what the file could hold, since the header's count may be corrupt. A
  // stream cannot be measured, so its reserve is capped, and the vector grows from there
  std::uint64_t limit = (file_size > 0) ? file_size / MIN_VARIANT_BYTES + 1
                                        : MAX_VARIANT_RESERVE;
  variants.reserve(std::min<std::uint64_t>(header.nvariants, limit));
  std::uint64_t offset = first_variant_offset();
  try {
    for (std::uint32_t idx=0; idx < header.nvariants; idx++) {
      variants.emplace_back(handle, offset, header.layout, header.compression,
                            header.nsamples, is_stdin);
      offset = variants.back().next_variant_offset;
    }
  } catch (const std::out_of_range &) {
    // the bgen holds fewer variants than its header lists
    std::uint64_t parsed = variants.size();
    variants.clear();
    throw std::invalid_argument("bgen is truncated - the header lists " +
                                std::to_string(header.nvariants) +
                                " variants, but only " + std::to_string(parsed) +
                                " could be read");
  } catch (...) {
    // drop the partial list, so a later call parses again
    variants.clear();
    throw;
  }
}

/// drop a subset of variants passed in by indexes
void CppBgenReader::drop_variants(const std::vector<int> &indices) {
  // check every index before dropping any, so a bad index leaves the variants as
  // they were
  for (auto idx : indices) {
    if ((idx < 0) || ((std::size_t) idx >= variants.size())) {
      throw std::out_of_range("variant index " + std::to_string(idx) +
                              " is out of range for a bgen with " +
                              std::to_string(variants.size()) + " variants");
    }
  }
  
  // mark the variants to drop, checking for duplicates
  std::vector<bool> drop(variants.size(), false);
  for (auto idx : indices) {
    if (drop[idx]) {
      throw std::invalid_argument("can't drop variants with duplicate indices");
    }
    drop[idx] = true;
  }

  // shift the kept variants down in one pass, which keeps them in file order
  std::size_t kept = 0;
  for (std::size_t i = 0; i < variants.size(); i++) {
    if (drop[i]) {
      continue;
    }
    if (kept != i) {
      variants[kept] = std::move(variants[i]);
    }
    kept++;
  }
  variants.erase(variants.begin() + kept, variants.end());
  variants.shrink_to_fit();
}

} // namespace bgen
