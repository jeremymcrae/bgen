
#include <fstream>
#include <vector>
#include <string>
#include <sstream>

namespace bgen {

class Samples {
public:
  Samples(std::ifstream & handle, int n_samples);
  Samples(std::string path, int n_samples);
  Samples() {};
  std::vector<std::string> samples;
};

} // namespace bgen
