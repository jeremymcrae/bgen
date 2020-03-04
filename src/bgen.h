#ifndef BGEN_BGEN_H_
#define BGEN_BGEN_H_

#include <fstream>
#include <vector>

#include "header.h"
#include "samples.h"
#include "variant.h"

namespace bgen {

class Bgen {
  std::ifstream handle;
public:
  Bgen(std::string path, std::string sample_path="");
  Variant next_var();
  std::vector<std::string> rsids();
  std::vector<std::string> chroms();
  std::vector<std::uint32_t> positions();
  std::vector<Variant> variants;
  Header header;
  Samples samples;
};

} // namespace bgen

#endif  // BGEN_BGEN_H_
