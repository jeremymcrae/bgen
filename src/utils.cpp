
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
