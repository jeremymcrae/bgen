
#include <fstream>
#include <vector>
#include <string>
#include <sstream>

namespace bgen {

class BgenSamples {
public:
  BgenSamples(std::ifstream & handle, int n_samples);
  BgenSamples(std::string path, int n_samples);
  BgenSamples() {};
  std::vector<std::string> samples;
};

} // namespace bgen
