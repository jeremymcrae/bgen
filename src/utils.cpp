
#include "utils.h"

namespace bgen {

unsigned char reverse(unsigned char b) {
  b = (b & 0xF0) >> 4 | (b & 0x0F) << 4;
  b = (b & 0xCC) >> 2 | (b & 0x33) << 2;
  b = (b & 0xAA) >> 1 | (b & 0x55) << 1;
  return b;
}

template<std::size_t N>
std::bitset<N> reverse(std::bitset<N> &b) {
  std::bitset<N> a;
  for(std::size_t i = 0; i < N/2; ++i) {
    bool t = b[i];
    a[i] = b[N-i-1];
    a[N-i-1] = t;
  }
  return a;
}

std::bitset<32> read_flags(std::ifstream & handle) {
  std::bitset<32> flags;
  for (int i=0; i < 4; i++) {
    char a;
    handle.read(&a, 1);
    std::bitset<8> tmp(reverse(a));
    for (int j = 0; j < 8; j++){
        flags[(i * 8) + j] = tmp[j];
      }
  }
  return reverse(flags);
}

std::vector<std::string> split(const std::string &s, char delim) {
  /* split a string by delimiter into a vector of elements
  */
  std::vector<std::string> elems;
  std::istringstream iss(s);
  std::string item;
  while (std::getline(iss, item, delim)) {
    elems.push_back(item);
  }
  return elems;
}

int BinomialCoefficient::n_choose_k(int n, int k) {
  /* calculate value of binomial coefficient, but only if it's not cached already
  
  This functtion depends on the ploidy and allele number, which will mostly be
  < 3, so this just caches the result after the first time it's computed for
  quick lookup.
  */
  if ( k > n - k ) {
    k = n - k;
  }
  std::vector<int> pair = {n, k};
  if (cached.count(pair) == 1) {
    return cached[pair];
  }
  
  int res = 1;
  
  // Calculate value of
  // [n * (n-1) *---* (n-k+1)] / [k * (k-1) *----* 1]
  for (int i = 0; i < k; ++i) {
    res *= (n - i);
    res /= (i + 1);
  }
  cached[pair] = res;
  
  return res;
}

// Returns value of Binomial Coefficient C(n, k)
int n_choose_k(int n, int k) {
  int res = 1;
  
  // Since C(n, k) = C(n, n-k)
  if ( k > n - k ) {
    k = n - k;
  }
  
  // Calculate value of
  // [n * (n-1) *---* (n-k+1)] / [k * (k-1) *----* 1]
  for (int i = 0; i < k; ++i) {
    res *= (n - i);
    res /= (i + 1);
  }
  
  return res;
}

} // namespace bgen
