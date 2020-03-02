#ifndef BGEN_UTILS_H_
#define BGEN_UTILS_H_

#include <bitset>
#include <fstream>
#include <map>
#include <string>
#include <sstream>
#include <vector>

namespace bgen {

unsigned char reverse(unsigned char b);
std::bitset<32> read_flags(std::ifstream & handle);
std::vector<std::string> split(const std::string &s, char delim);
int n_choose_k(int n, int k);

class BinomialCoefficient {
  std::map<std::vector<int>, int> cached;
public:
  BinomialCoefficient() {};
  int n_choose_k(int n, int k);
};

} // namespace bgen

#endif  // BGEN_UTILS_H_
