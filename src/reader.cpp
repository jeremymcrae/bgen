
#include <algorithm>
#include <cerrno>
#include <climits>
#include <cstring>
#include <string>

#if defined(_WIN32)
  #include <fcntl.h>
  #include <io.h>
#else
  #include <unistd.h>
#endif

#include "reader.h"

namespace bgen {

/// bytes the bgen stream reads at a time (even on seek)
///
/// A seek discards the previous buffer, and the next read refills it. Every
/// variant seeks when starting, and by default would read 8 KB, despite only
/// using ~50 bytes for variant metadata. Resize buffer to 512 bytes instead.
///
/// Reading genotypes is unaffected, since that reads > 512 bytes.
const std::size_t STREAM_BUFFER = 512;

/// an ifstream which owns the buffer it reads through
///
/// The buffer has to outlive the stream, and the stream outlives the reader whenever a
/// Variant still holds it. StreamBuffer is a base class so that it is constructed
/// before, and destroyed after, the stream itself.
struct StreamBuffer {
  std::vector<char> data;
  StreamBuffer(std::size_t size) : data(size) {}
};

struct BufferedFile : private StreamBuffer, public std::ifstream {
  BufferedFile(const std::string & path, std::size_t size) : StreamBuffer(size) {
    // pubsetbuf only has an effect before the file is opened
    rdbuf()->pubsetbuf(data.data(), (std::streamsize) data.size());
    open(path, std::ios::in | std::ios::binary);
  }
};

/// bytes the stdin stream reads at a time (bigger than the file buffer, since a
/// stream never seeks, so the only cost of a large buffer is the read-ahead)
const std::size_t STDIN_BUFFER = 1 << 16;

/// duplicate a descriptor, in the binary mode the bgen has to be read in. Windows
/// opens stdin in text mode, which mangles CRLF pairs and stops at the first 0x1a
static int dup_binary(int source) {
#if defined(_WIN32)
  int fd = _dup(source);
  if (fd >= 0) {
    _setmode(fd, _O_BINARY);
  }
  return fd;
#else
  return ::dup(source);
#endif
}

static void close_descriptor(int fd) {
#if defined(_WIN32)
  _close(fd);
#else
  ::close(fd);
#endif
}

/// a streambuf which reads a bgen from a file descriptor, since stdin cannot be
/// opened by path everywhere, and std::cin belongs to the process, not to a reader
class DescriptorBuf : public std::streambuf {
public:
  DescriptorBuf(int source, std::size_t size) : data(size) {
    fd = dup_binary(source);
    // start with an empty get area, so the first read fills it
    setg(data.data(), data.data(), data.data());
  }
  ~DescriptorBuf() { close(); }
  bool is_open() const { return fd >= 0; }
  void close() {
    if (fd >= 0) {
      close_descriptor(fd);
      fd = -1;
    }
  }
protected:
  int_type underflow() override {
    if (gptr() >= egptr()) {
      std::streamsize taken = fill(data.data(), data.size());
      if (taken <= 0) {
        return traits_type::eof();
      }
      setg(data.data(), data.data(), data.data() + taken);
    }
    return traits_type::to_int_type(*gptr());
  }
  /// read a run of bytes, going straight to the descriptor for runs longer than the
  /// buffer, which the default would otherwise split into a read per bufferful
  std::streamsize xsgetn(char * dest, std::streamsize n) override {
    std::streamsize taken = 0;
    while (taken < n) {
      std::streamsize buffered = egptr() - gptr();
      if (buffered > 0) {
        std::streamsize used = std::min(buffered, n - taken);
        std::memcpy(dest + taken, gptr(), (std::size_t) used);
        gbump((int) used);
        taken += used;
      } else if ((n - taken) >= (std::streamsize) data.size()) {
        std::streamsize got = fill(dest + taken, (std::size_t) (n - taken));
        if (got <= 0) {
          break;
        }
        taken += got;
      } else if (traits_type::eq_int_type(underflow(), traits_type::eof())) {
        break;
      }
    }
    return taken;
  }
private:
  /// read from the descriptor, giving the bytes read, 0 at the end, or -1 on error. A
  /// pipe stops at what has been written, so callers ask again on a short read
  std::streamsize fill(char * dest, std::size_t n) {
    if (fd < 0) {
      return -1;
    }
    while (true) {
#if defined(_WIN32)
      int taken = _read(fd, dest, (unsigned int) std::min(n, (std::size_t) INT_MAX));
#else
      ssize_t taken = ::read(fd, dest, n);
#endif
      if (taken >= 0) {
        return (std::streamsize) taken;
      }
      if (errno != EINTR) {
        return -1;
      }
    }
  }
  int fd = -1;
  std::vector<char> data;
};

/// an istream owning the descriptor it reads and the buffer it reads through. The
/// buffer is a base class for the same reason as StreamBuffer above
struct DescriptorBuffer {
  DescriptorBuf buf;
  DescriptorBuffer(int fd, std::size_t size) : buf(fd, size) {}
};

struct DescriptorStream : private DescriptorBuffer, public std::istream {
  DescriptorStream(int fd, std::size_t size) : DescriptorBuffer(fd, size),
                                               std::istream(&buf) {
    if (!buf.is_open()) {
      setstate(std::ios::failbit);
    }
  }
  /// release the descriptor, while leaving whatever was buffered readable
  void close() { buf.close(); }
};

/// smallest number of bytes a variant can occupy in a bgen
///
/// This is a deliberate underestimate of the per variant fields (the identifiers
/// and alleles all have length prefixes, and layout 1 has fewer fields), since it
/// is only used to cap how much space is reserved up front. Undercounting just
/// makes the reserve smaller than it could be.
const std::uint64_t MIN_VARIANT_BYTES = 12;

/// most variants to reserve space for up front when the bgen's size is unknown
///
/// A Variant is a few hundred bytes, so reserving the claimed count asks the
/// allocator for over a terabyte when that count is corrupt, which fails even
/// though the pages are never touched. This is only needed for a stream, whose
/// size cannot be measured, so growing past it costs a few reallocations.
const std::uint64_t MAX_VARIANT_RESERVE = 1 << 16;

CppBgenReader::CppBgenReader(std::string path, std::string sample_path, bool delay_parsing) {
  is_stdin = (path == "/dev/stdin");
  // stdin is read from its descriptor rather than held as std::cin, so that the error
  // state and the buffered read-ahead belong to this reader alone
  if (is_stdin) {
    handle = std::shared_ptr<std::istream>(new DescriptorStream(0, STDIN_BUFFER));
  } else {
    handle = std::shared_ptr<std::istream>(new BufferedFile(path, STREAM_BUFFER));
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
    // Reject a wildly wrong count before the sample block is even read. Samples()
    // applies a tighter check of its own, so this is not load bearing, but it keeps
    // the coarse failure separate from a block that merely disagrees with its IDs.
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
    samples = Samples(header.nsamples, file_size);
  }
  
  offset = first_variant_offset();
  if (!delay_parsing) {
    parse_all_variants();
  }
}

/// release the operating system file handle
///
/// The stream is shared with every Variant taken from this reader, so dropping this
/// reader's reference is not enough to close the file while any of those survive.
/// Windows refuses to delete a file that is still open, so closing has to reach the
/// handle itself. Every accessor on a Variant checks the reader is open before it
/// reads, so none of them can reach the closed stream afterwards.
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
  // Only reserve what the file could hold. nvariants comes straight from the
  // header, so a corrupt count would otherwise ask for hundreds of GB of
  // Variants up front. Parsing still stops at whatever the file really contains.
  std::uint64_t n_reserve = header.nvariants;
  if (file_size > 0) {
    std::uint64_t possible = file_size / MIN_VARIANT_BYTES + 1;
    if (n_reserve > possible) {
      n_reserve = possible;
    }
  } else {
    // A stream cannot be seeked, so nothing here bounds the count. Cap the reserve
    // and let the vector grow as the variants actually arrive.
    n_reserve = std::min(n_reserve, MAX_VARIANT_RESERVE);
  }
  variants.reserve(n_reserve);
  try {
    for (std::uint32_t idx=0; idx < header.nvariants; idx++) {
      variants.push_back(next_var());
    }
  } catch (const std::out_of_range &) {
    // Running out of file part way through means the bgen holds fewer variants
    // than its header lists, so say that, rather than reporting the end of the
    // file for what looked like a perfectly valid variant index.
    std::uint64_t parsed = variants.size();
    variants.clear();
    throw std::invalid_argument("bgen is truncated - the header lists " +
                                std::to_string(header.nvariants) +
                                " variants, but only " + std::to_string(parsed) +
                                " could be read");
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
