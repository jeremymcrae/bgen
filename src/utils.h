#ifndef BGEN_UTILS_H_
#define BGEN_UTILS_H_

#include <bitset>
#include <fstream>
#include <map>
#include <string>
#include <sstream>
#include <vector>

namespace bgen {

std::vector<std::string> split(const std::string &s, char delim);
int n_choose_k(int n, int k);

} // namespace bgen

#endif  // BGEN_UTILS_H_
