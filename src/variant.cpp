
#include <stdexcept>

#include "variant.h"

namespace bgen {

Variant::Variant(std::ifstream & handle, int layout, int compression) {
  /* initialise a single variant with chrom, pos, rsID identifiers
  
  This starts a Genotypes object, but this doesn't parse the genotypes until
  required, just starts it so we can get the offset of the next variant, so as
  to parse the bgen variants at speed.
  */
  offset = handle.tellg();
  if (layout == 1) {
    handle.read(reinterpret_cast<char*>(&n_samples), sizeof(n_samples));
  }
  
  // get the variant ID (first need to know how long the field is)
  std::uint16_t item_len;
  handle.read(reinterpret_cast<char*>(&item_len), sizeof(std::uint16_t));
  std::copy_n(std::istream_iterator<char>(handle), item_len, std::back_inserter(varid));
  
  // get the rsID (first need to know how long the field is)
  handle.read(reinterpret_cast<char*>(&item_len), sizeof(std::uint16_t));
  
  std::copy_n(std::istream_iterator<char>(handle), item_len, std::back_inserter(rsid));
  
  // get the chromosome (first need to know how long the field is)
  handle.read(reinterpret_cast<char*>(&item_len), sizeof(std::uint16_t));
  std::copy_n(std::istream_iterator<char>(handle), item_len, std::back_inserter(chrom));
  
  handle.read(reinterpret_cast<char*>(&pos), sizeof(std::uint32_t));
  if (layout == 1) {
    n_alleles = 2;
  } else {
    handle.read(reinterpret_cast<char*>(&n_alleles), sizeof(std::uint16_t));
  }
  
  for (int x=0; x < n_alleles; x++) {
    std::uint32_t allele_len;
    std::string allele;
    handle.read(reinterpret_cast<char*>(&allele_len), sizeof(std::uint32_t));
    std::copy_n(std::istream_iterator<char>(handle), allele_len, std::back_inserter(allele));
    alleles.push_back(allele);
  }
  
  geno = Genotypes(&handle, layout, compression, n_alleles, n_samples);
}

std::uint64_t Variant::next_variant_offset() {
  /* uses the genotypes object to find the offset of the next variant
  */
  return geno.next_var_offset;
}

std::string Variant::alt() {
  throw std::invalid_argument("haven't completed alt() fucntion yet");
}

std::vector<std::vector<double>> Variant::genotypes() {
  return geno.genotypes();
}

std::vector<double> Variant::alt_dosage() {
  throw std::invalid_argument("haven't completed alt_dosage() fucntion yet");
}

} // namespace bgen
