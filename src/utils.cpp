
#include "utils.h"

namespace bgen {

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
  
  This function depends on the ploidy and allele number, which will mostly be
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

} // namespace bgen
