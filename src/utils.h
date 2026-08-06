#ifndef BGEN_UTILS_H_
#define BGEN_UTILS_H_

#include <algorithm>
#include <bitset>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <istream>
#include <map>
#include <memory>
#include <string>
#include <sstream>
#include <vector>

namespace bgen {

/// wrap a stream which we do not own in a shared_ptr
///
/// The bgen stream is shared between a CppBgenReader and every Variant opened
/// from it, so that it stays open for as long as any of them still need it. A
/// couple of streams are not ours to close though: std::cin, and the stream
/// belonging to an unpickled variant (which is owned by whichever reader is
/// still holding it). Those get a shared_ptr with a no-op deleter, so they can
/// be stored in the same way without being closed from under their owner.
inline std::shared_ptr<std::istream> borrowed_stream(std::istream * handle) {
  return std::shared_ptr<std::istream>(handle, [](std::istream *) {});
}

// Mark a function as being compiled for a specific x86_64 instruction set, so
// that SIMD code can be built without those instructions being enabled for the
// translation unit as a whole. Without this, the compiler is free to use e.g.
// AVX2 in functions which are only reached after a runtime check has found that
// AVX2 is unavailable, which crashes with an illegal instruction.
//
// The SIMD code is all guarded by '#if defined(__x86_64__)', which GCC and clang
// define but MSVC does not, so MSVC never reaches these and needs no attribute.
#if defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__))
  #define BGEN_TARGET_AVX  __attribute__((target("avx")))
  #define BGEN_TARGET_AVX2 __attribute__((target("avx2")))
  #define BGEN_TARGET_SSE4 __attribute__((target("sse4.1")))
#else
  #define BGEN_TARGET_AVX
  #define BGEN_TARGET_AVX2
  #define BGEN_TARGET_SSE4
#endif

/// read a fixed width value from a stream, and report whether the read worked
///
/// A failed read leaves the value indeterminate, so reads must be checked before
/// the value gets used, otherwise we branch on (or size allocations from)
/// uninitialised memory. Reading into the value directly also avoids an
/// unaligned reinterpret_cast.
///
/// Callers throw, since they need different exceptions: a variant throws
/// out_of_range (IndexError, which ends iteration), whereas a truncated header
/// or sample block means the file is unusable, so raises ValueError.
template <typename T>
inline bool read_value(std::istream & handle, T & value) {
  return (bool) handle.read(reinterpret_cast<char *>(&value), sizeof(T));
}

/// max chunk size per file read at a time for a string
const std::size_t STRING_READ_CHUNK = 1 << 20;

/// read a string prefixed by its length, and report whether the reads worked
///
/// Reads exactly len bytes, defined by the length type. A std::istream_iterator
/// cannot be used here, as it does formatted input and so skips whitespace,
/// corrupting any ID containing a space and misaligning the stream after it.
///
/// Long strings are read in chunks, which costs a few extra resizes on genuinely
/// long alleles but means a corrupt length only allocates what the file can supply.
template <typename LenType>
inline bool read_prefixed_string(std::istream & handle, std::string & value) {
  LenType len;
  if (!read_value(handle, len)) {
    return false;
  }
  std::size_t remaining = (std::size_t) len;
  if (remaining <= STRING_READ_CHUNK) {
    // short enough that the claim cannot amplify much, so read it in one go
    value.resize(remaining);
    return (remaining == 0) || (bool) handle.read(&value[0], remaining);
  }
  value.clear();
  while (remaining > 0) {
    std::size_t chunk = std::min(remaining, STRING_READ_CHUNK);
    std::size_t start = value.size();
    value.resize(start + chunk);
    if (!handle.read(&value[start], chunk)) {
      return false;
    }
    remaining -= chunk;
  }
  return true;
}

struct Range {
    std::uint8_t _min;
    std::uint8_t _max;
};

std::uint32_t n_choose_k(int n, int k);
bool minor_certain(double freq, int n_checked, double z);
std::uint64_t fast_ploidy_sum(std::uint8_t * x, std::uint32_t & size);
Range fast_range(std::uint8_t * x, std::uint32_t & size);
void fast_missing_scan(const char * data, std::uint32_t size,
                       std::vector<std::uint32_t> & missing);
double fast_dosage_sum(float * dose, std::uint32_t size, std::uint64_t & n_called);
double fast_strided_sum(const float * dose, std::uint32_t n_samples, std::uint32_t start,
                        std::uint32_t stride, std::uint64_t & n_checked);

} // namespace bgen

#endif  // BGEN_UTILS_H_
