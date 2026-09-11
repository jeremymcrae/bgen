
#include <algorithm>
#include <cerrno>
#include <climits>
#include <cstring>

#if defined(_WIN32)
  #include <fcntl.h>
  #include <io.h>
#else
  #include <unistd.h>
#endif

#include "streaming.h"

namespace bgen {

/// bytes the bgen file reads at a time (even on seek)
///
/// A seek discards the previous buffer, and the next read refills it. Every
/// variant seeks when starting, and by default would read 8 KB, despite only
/// using ~50 bytes for variant metadata. Resize buffer to 512 bytes instead.
///
/// Reading genotypes is unaffected, since that reads > 512 bytes.
const std::size_t FILE_BUFFER = 512;

/// bytes the stdin stream reads at a time (bigger than the file buffer, since a
/// stream never seeks, so the only cost of a large buffer is the read-ahead)
const std::size_t STDIN_BUFFER = 1 << 16;

BufferedFile::BufferedFile(const std::string & path) : FileBuffer(FILE_BUFFER) {
  // pubsetbuf only has an effect before the file is opened
  rdbuf()->pubsetbuf(data.data(), (std::streamsize) data.size());
  open(path, std::ios::in | std::ios::binary);
}

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

DescriptorBuf::DescriptorBuf(int source) : data(STDIN_BUFFER) {
  fd = dup_binary(source);
  // start with an empty get area, so the first read fills it
  setg(data.data(), data.data(), data.data());
}

void DescriptorBuf::close() {
  if (fd >= 0) {
    close_descriptor(fd);
    fd = -1;
  }
}

DescriptorBuf::int_type DescriptorBuf::underflow() {
  if (gptr() >= egptr()) {
    std::streamsize taken = fill(data.data(), data.size());
    if (taken <= 0) {
      return traits_type::eof();
    }
    setg(data.data(), data.data(), data.data() + taken);
  }
  return traits_type::to_int_type(*gptr());
}

std::streamsize DescriptorBuf::xsgetn(char * dest, std::streamsize n) {
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

std::streamsize DescriptorBuf::fill(char * dest, std::size_t n) {
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

} // namespace bgen
