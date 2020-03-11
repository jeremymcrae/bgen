
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
  
  if ((int) n_samples != expected_n) {
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

std::vector<std::vector<float>> Variant::probabilities() {
  /* get genotype probabilites for the variant
  
  Note this converts the probability arrays to a 2D vector, which is much slower
  than the initial parsing.
  */
  float * float_probs = geno.probabilities();
  std::vector<std::vector<float>> probs(n_samples, std::vector<float>(geno.max_probs));
  int offset;
  for (uint n=0; n<n_samples; n++) {
    offset = n * geno.max_probs;
    for (int x=0; x<geno.max_probs; x++){
      probs[n][x] = float_probs[offset + x];
    }
  }
  return probs;
}

void Variant::dosages(float * first, float * second) {
  /* get allele dosages (assumes biallelic variant)
  */
  if (n_alleles != 2) {
    throw std::invalid_argument("can't get allele dosages for non-biallelic var.");
  }
  
  float * probs = geno.probabilities();
  
  int offset;
  float sums[2] = {0, 0};
  float ploidy = geno.max_ploidy;
  float half_ploidy = ploidy / 2;
  for (uint n=0; n<n_samples; n++) {
    offset = n * geno.max_probs;
    if (!geno.constant_ploidy) {
      ploidy = (float) geno.ploidy[n];
      half_ploidy = ploidy / 2;
    }
    float halved = probs[offset + 1] * half_ploidy;
    first[n] = (probs[offset] * ploidy) + halved;
    second[n] = (probs[offset + 2] * ploidy) + halved;
    
    sums[0] += first[n];
    sums[1] += second[n];
  }
  
  if (sums[0] < sums[1]) {
    minor_idx = 0;
  } else if (sums[1] < sums[0]) {
    minor_idx = 1;
  } else {
    minor_idx = 0; // pick the first if the alelles are 50:50
  }
}

std::vector<float> & Variant::minor_allele_dosage() {
  /* get dosage of the minor allele (only works for biallelic variants)
  */
  float * first = new float[n_samples];
  float * second = new float[n_samples];
  dosages(first, second);
  
  minor_allele = alleles[minor_idx];
  if (minor_idx == 0) {
    minor_allele_dose = std::vector<float>(first, first + n_samples);
  } else {
    minor_allele_dose = std::vector<float>(second, second + n_samples);
  }
  
  delete []first;
  delete []second;
  
  return minor_allele_dose;
}

} // namespace bgen
