
#include <fstream>
#include <vector>
#include <string>

namespace bgen {

class BgenSamples {
public:
  BgenSamples(std::ifstream & handle, int n_samples);
  BgenSamples() {};
  std::vector<std::string> samples;
};

} // namespace bgen
