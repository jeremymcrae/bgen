
#include "utils.h"

namespace bgen {

// Returns value of Binomial Coefficient C(n, k)
uint n_choose_k(int n, int k) {
  uint res = 1;

  // Since C(n, k) = C(n, n-k)
  if ( k > n - k ) {
    k = n - k;
  }

  // Calculate value of
  // [n * (n-1) *---* (n-k+1)] / [k * (k-1) *----* 1]
  for (uint i = 0; i < (uint) k; ++i) {
    res *= (n - i);
    res /= (i + 1);
  }

  return res;
}

} // namespace bgen
