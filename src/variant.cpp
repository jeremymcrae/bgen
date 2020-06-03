
#include <stdexcept>
#include <cmath>

#include "variant.h"

namespace bgen {

Variant::Variant(std::ifstream & handle, std::uint64_t & varoffset, int layout, int compression, int expected_n) {
  /* initialise a single variant with chrom, pos, rsID identifiers
  
  This starts a Genotypes object, but this doesn't parse the genotypes until
  required, just starts it so we can get the offset of the next variant, so as
  to parse the bgen variants at speed.
  */
  offset = varoffset;
  handle.seekg(offset);
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
  if (item_len > 0) {
    std::copy_n(std::istream_iterator<char>(handle), item_len, std::back_inserter(varid));
  }
  
  // get the rsID (first need to know how long the field is)
  handle.read(reinterpret_cast<char*>(&item_len), sizeof(std::uint16_t));
  if (item_len > 0) {
    std::copy_n(std::istream_iterator<char>(handle), item_len, std::back_inserter(rsid));
  }
  
  // get the chromosome (first need to know how long the field is)
  handle.read(reinterpret_cast<char*>(&item_len), sizeof(std::uint16_t));
  if (item_len > 0) {
    std::copy_n(std::istream_iterator<char>(handle), item_len, std::back_inserter(chrom));
  }
  
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

int Variant::probs_per_sample() {
  return geno.max_probs;
}

bool Variant::phased() {
  if (geno.max_probs == 0) {
    throw std::invalid_argument("unknown phase, run variant.probabilities() first");
  }
  return geno.phased;
}

std::uint8_t * Variant::ploidy() {
  if (geno.max_probs == 0) {
    throw std::invalid_argument("unknown ploidy, run variant.probabilities() first");
  }
  return geno.ploidy;
}

float * Variant::probs_1d() {
  /* get genotype probabilities for the variant as a 1-dimensional vector
  
  This makes it easy to pass the data via cython into a numpy array, which can
  be reshaped to a 2-D array.
  */
  return geno.probabilities();
}

std::vector<std::vector<float>> & Variant::probabilities() {
  /* get genotype probabilities for the variant
  
  Note this converts the probability arrays to a 2D vector, which is much slower
  than the initial parsing.
  */
  float * float_probs = geno.probabilities();
  int nrows = 0;
  if (geno.phased) {
    for (uint n=0; n<n_samples; n++) {
      nrows += geno.ploidy[n];
    }
  } else {
    nrows = n_samples;
  }
  
  probs2d = std::vector<std::vector<float>>(nrows, std::vector<float>(geno.max_probs));
  int offset;
  for (uint n=0; n<n_samples; n++) {
    offset = n * geno.max_probs;
    for (int x=0; x<geno.max_probs; x++){
      probs2d[n][x] = float_probs[offset + x];
    }
  }
  return probs2d;
}

bool minor_certain(double freq, int n_checked, double z) {
    /** check if the minor allele is certain (to 99.9999999999999& confidence)
    
    Take the frequency, and number of individuals checked so far, and see if the
    99.99..(fifteen nines) confidence interval overlaps 0.5. If not, then we can
    be sure we've identified the minor allele, even without checking the full
    population.
    
      @freq estimated minor allele frequency
      @n_checked number of individsuals checked so far
      @z standard normal deviate (eg 1.96 for 95% CI, here we use 10.0 for
        stronger confidence, and the fact the normal approximation for confidence
        intervals isn't perfect)
      @return True/False for whether to halt the permuations
    */
    double delta = (z * std::sqrt((freq * (1 - freq))))/n_checked;
    
    // check if the confidence interval overlaps 0.5
    return !((freq - delta < 0.5) & (freq + delta > 0.5));
}

void Variant::dosages() {
  /* get allele dosages (assumes biallelic variant)
  */
  if (n_alleles != 2) {
    throw std::invalid_argument("can't get allele dosages for non-biallelic var.");
  }
  
  dose = new float[n_samples];
  float * probs = geno.probabilities();
  
  int offset;
  int batchsize = 1000;
  int increment = n_samples / batchsize;
  float sums[2] = {0, 0};
  float ploidy = geno.max_ploidy;
  float half_ploidy = ploidy / 2;
  float halved;
  int total;
  
  // rather than checking every individual to see which is the minor allele, we
  // check subsets, in batches of 1000. We obtain alleles for individuals in the
  // batch, then check if a confidence interval for the frequency of the less
  // frequent allele could overlap 0.5. If not, we can be reasonably certain the
  // less frequent allele is the true minor allele, without having to check the
  // full cohort. This can be 60X faster than checking the full cohort in larger
  // populations.
  
  // To make sure we don't hit weird groupings of alleles in individuals, this
  // picks samples uniformly thoughout the population, by using an appropriate
  // step size.
  for (int idx=0; idx<increment; idx++) {
    for (uint n=idx; n<n_samples; n += increment) {
      offset = n * geno.max_probs;
      if (!geno.constant_ploidy) {
        ploidy = (float) geno.ploidy[n];
        half_ploidy = ploidy / 2;
      }
      halved = probs[offset + 1] * half_ploidy;
      sums[0] += (probs[offset] * ploidy) + halved;
      sums[1] += (probs[offset + 2] * ploidy) + halved;
    }
    total = sums[0] + sums[1];
    double freq = (double) std::min(sums[0], sums[1]) / total;
    if (minor_certain(freq, batchsize * (idx + 1), 10.0)) {
      break;
    }
  }
  
  int geno_idx = 0;
  if (sums[0] < sums[1]) {
    minor_idx = 0;
  } else if (sums[1] < sums[0]) {
    minor_idx = 1;
    geno_idx = 2;
  } else {
    minor_idx = 0; // pick the first if the alelles are 50:50
  }
  
  // now that we know which allele to use, calculate dosage for all samples
  if (geno.constant_ploidy) {
    for (uint n=0; n<n_samples; n++) {
      offset = n * geno.max_probs;
      dose[n] = (probs[offset + geno_idx] * ploidy) + probs[offset + 1] * half_ploidy;
    }
  } else {
    for (uint n=0; n<n_samples; n++) {
      offset = n * geno.max_probs;
      ploidy = (float) geno.ploidy[n];
      half_ploidy = ploidy / 2;
      dose[n] = (probs[offset + geno_idx] * ploidy) + probs[offset + 1] * half_ploidy;
    }
  }
}

float * Variant::minor_allele_dosage() {
  /* get dosage of the minor allele (only works for biallelic variants)
  */
  clear_probs(); // clean up so repeated calls don't leak memory
  
  dosages();
  
  minor_allele = alleles[minor_idx];
  return dose;
}

void Variant::clear_probs() {
  if (minor_idx != -1) {
    delete[] dose;
  }
  minor_idx = -1;
}

} // namespace bgen
