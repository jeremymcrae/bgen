#ifndef BGEN_UTILS_H_
#define BGEN_UTILS_H_

#include <bitset>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <map>
#include <string>
#include <sstream>
#include <vector>

namespace bgen {

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

struct Range {
    std::uint8_t _min;
    std::uint8_t _max;
};

std::uint32_t n_choose_k(int n, int k);
bool minor_certain(double freq, int n_checked, double z);
std::uint64_t fast_ploidy_sum(std::uint8_t * x, std::uint32_t & size);
Range fast_range(std::uint8_t * x, std::uint32_t & size);

} // namespace bgen

#endif  // BGEN_UTILS_H_
