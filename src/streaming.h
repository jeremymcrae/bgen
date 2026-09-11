#ifndef BGEN_STREAMING_H_
#define BGEN_STREAMING_H_

#include <cstddef>
#include <fstream>
#include <istream>
#include <streambuf>
#include <string>
#include <vector>

namespace bgen {

/// an ifstream which owns the buffer it reads through
///
/// The buffer has to outlive the stream, and the stream outlives the reader whenever a
/// Variant still holds it. FileBuffer is a base class so that it is constructed
/// before, and destroyed after, the stream itself.
struct FileBuffer {
  std::vector<char> data;
  FileBuffer(std::size_t size) : data(size) {}
};

/// a specialized system for reading files to speed up parsing variants metadata
///
/// Constructing a variant first seeks to the offset, then reads file data to
/// obtain variant metadata. By default this would read 8 kb, far more than the
/// ~50 bytes required for metadata. This reads into a buffer of 512 bytes instead.
struct BufferedFile : private FileBuffer, public std::ifstream {
  BufferedFile(const std::string & path);
};

/// a streambuf which reads a bgen from a file descriptor, since stdin cannot be
/// opened by path everywhere, and std::cin belongs to the process, not to a reader
class DescriptorBuf : public std::streambuf {
public:
  DescriptorBuf(int source);
  ~DescriptorBuf() { close(); }
  bool is_open() const { return fd >= 0; }
  void close();
protected:
  int_type underflow() override;
private:
  /// read from the descriptor, giving the bytes read, 0 at the end, or -1 on error. A
  /// pipe stops at what has been written, so callers ask again on a short read
  std::streamsize fill(char * dest, std::size_t n);
  int fd = -1;
  std::vector<char> data;
};

// struct for opening bgen from stdin. This is more complex than it would seem,
// since we can't use std::cin, and want it to work on windows/linux/macosx.
//
// The buffer is a base class for the same reason as FileBuffer above, so that it is
// constructed before, and destroyed after, the stream itself.
struct DescriptorStream : private DescriptorBuf, public std::istream {
  DescriptorStream(int fd) : DescriptorBuf(fd), std::istream(this) {
    if (!is_open()) {
      setstate(std::ios::failbit);
    }
  }
  /// release the descriptor, while leaving whatever was buffered readable
  using DescriptorBuf::close;
};

} // namespace bgen

#endif  // BGEN_STREAMING_H_
