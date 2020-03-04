
#include <stdexcept>

#include "variant.h"

namespace bgen {

Variant::Variant(std::ifstream & handle, int layout, int compression, int expected_n) {
  /* initialise a single variant with chrom, pos, rsID identifiers
  
  This starts a Genotypes object, but this doesn't parse the genotypes until
  required, just starts it so we can get the offset of the next variant, so as
  to parse the bgen variants at speed.
  */
  offset = handle.tellg();
  if (layout == 1) {
    handle.read(reinterpret_cast<char*>(&n_samples), sizeof(n_samples));
  } else {
    n_samples = expected_n;
  }
  
  if (n_samples != expected_n) {
    throw std::invalid_argument("number of samples doesn't match");
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

int Variant::alt_index(std::vector<std::vector<float>> & dose) {
  /* find the index offset for the alt allele (second most frequent allele)
  */
  // get the per column sums. Assumes only two alleles, probably ok in practise
  std::vector<double> sums(2);
  for (int n=0; n<n_samples; n++) {
    sums[0] += dose[n][0];
    sums[1] += dose[n][1];
  }
  
  if (sums[0] > sums[1]) {
    return 0;
  } else if (sums[1] > sums[0]) {
    return 1;
  } else {
    return 1; // pick the first if the alelles are 50:50
  }
}

std::vector<std::vector<float>> & Variant::probabilities() {
  return geno.probabilities();
}

void Variant::dosages(std::vector<std::vector<float>> & dose) {
  /* get allele dosages (assumes biallelic variant)
  */
  if (n_alleles != 2) {
    throw std::invalid_argument("can't get allele dosages for non-biallelic var.");
  }
  
  auto & probs = probabilities();
  
  std::uint8_t ploidy;
  for (int n=0; n<n_samples; n++) {
    ploidy = geno.ploidy[n];
    float halved = probs[n][1] * ((float) ploidy / 2);
    dose[n][0] = (probs[n][0] * ploidy) + halved;
    dose[n][1] = (probs[n][2] * ploidy) + halved;
  }
}

std::vector<float> & Variant::alt_dosage() {
  dose = std::vector<std::vector<float>>(n_samples, std::vector<float>(2, 0));
  dosages(dose);
  int idx = alt_index(dose);
  
  alt_dose = std::vector<float>(n_samples);
  for (int n=0; n<n_samples; n++) {
    alt_dose[n] = dose[n][idx];
  }
  
  return alt_dose;
}

} // namespace bgen
