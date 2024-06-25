### Another bgen parser
![bgen](https://github.com/jeremymcrae/bgen/workflows/bgen/badge.svg)

This is a package for reading and writing [bgen files](https://www.well.ox.ac.uk/~gav/bgen_format).

This package uses cython to wrap c++ code for parsing bgen files. It's fairly
quick, it can parse genotypes from 500,000 individuals at ~300 variants per
second within a single python process (~450 million probabilities per second
with a 3GHz CPU). Decompressing the genotype probabilities can be the slow step,
zlib decompression takes 80% of the total time, using zstd compressed genotypes
is ~2X faster.

This has been optimized for UKBiobank bgen files (i.e. bgen version 1.2 with
zlib compressed 8-bit genotype probabilities, but the other bgen versions and
zstd compression have also been tested using example bgen files).

#### Install
`pip install bgen`

#### Usage
```python
from bgen import BgenReader, BgenWriter

bfile = BgenReader(BGEN_PATH)
rsids = bfile.rsids()

# select a variant by indexing
var = bfile[1000]

# pull out genotype probabilities
probs = var.probabilities  # returns 2D numpy array
dosage = var.minor_allele_dosage  # returns 1D numpy array for biallelic variant

# iterate through every variant in the file
with BgenReader(BGEN_PATH, delay_parsing=True) as bfile:
  for var in bfile:
      dosage = var.minor_allele_dosage

# get all variants in a genomic region
variants = bfile.fetch('21', 10000, 5000000)

# or for writing bgen files
import numpy as np
from bgen import BgenWriter

geno = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]]).astype(np.float64)
with BgenWriter(BGEN_PATH, n_samples=3) as bfile:
  bfile.add_variant(varid='var1', rsid='rs1', chrom='chr1', pos=1,
                    alleles=['A', 'G'], genotypes=geno)
```

#### API documentation

``` py
class BgenReader(path, sample_path='', delay_parsing=False)
    # opens a bgen file. If a bgenix index exists for the file, the index file
    # will be opened automatically for quicker access of specific variants.
    Arguments:
      path: path to bgen file
      sample_path: optional path to sample file. Samples will be given integer IDs
          if sample file is not given and sample IDs not found in the bgen file
      delay_parsing: True/False option to allow for not loading all variants into
          memory when the BgenFile is opened. This can save time when iterating
          across variants in the file
  
  Attributes:
    samples: list of sample IDs
    header: BgenHeader with info about the bgen version and compression.
  
  Methods:
    slicing: BgenVars can be accessed by slicing the BgenFile e.g. bfile[1000]
    iteration: variants in a BgenFile can be looped over e.g. for x in bfile: print(x)
    fetch(chrom, start=None, stop=None): get all variants within a genomic region
    drop_variants(list[int]): drops variants by index from being used in analyses
    with_rsid(rsid): returns list of BgenVars with given rsid
    at_position(pos): returns list of BgenVars at a given position
    varids(): returns list of varids for variants in the bgen file.
    rsids(): returns list of rsids for variants in the bgen file.
    chroms(): returns list of chromosomes for variants in the bgen file.
    positions(): returns list of positions for variants in the bgen file.

class BgenVar(handle, offset, layout, compression, n_samples):
  # Note: this isn't called directly, but instead returned from BgenFile methods
  Attributes:
    varid: ID for variant
    rsid: reference SNP ID for variant
    chrom: chromosome variant is on
    pos: nucleotide position variant is at
    alleles: list of alleles for variant
    is_phased: True/False for whether variant has phased genotype data
    ploidy: list of ploidy for each sample. Samples are ordered as per BgenFile.samples
    minor_allele: the least common allele (for biallelic variants)
    minor_allele_dosage: 1D numpy array of minor allele dosages for each sample
    alt_dosage: 1D numpy array of alt allele dosages for each sample
    probabilities:  2D numpy array of genotype probabilities, one sample per row
      These are most likely for biallelic diploid variants. In that scenario
      unphased probabilities have three columns, for homozygous first allele 
      (AA), heterozygous (Aa), homozygous second allele (aa).
      In contrast, phased probabilities (for a biallelic diploid variant) would
      have four columns, first two for haplotype 1 (hap1-allele1, hap1-allele2), 
      last two for haplotype 2 (hap2-allele1, hap2-allele2).
  
  BgenVars can be pickled e.g. pickle.dumps(var)


class BgenWriter(path, n_samples, samples=[], compression='zstd' layout=2, metadata=None)
    # opens a bgen file to write variants to. Automatically makes a bgenix index file
    Arguments:
      path: path to write data to
      n_samples: number of samples that you have data for
      samples: list of sample IDs (same length as n_samples)
      compression: compression type: None, 'zstd', or 'zlib' (default='zstd')
      layout: bgen layout format (default=2)
      metadata: any additional metadata you want o include in the file (as str)
    
    Methods:
      add_variant_direct(variant)
        Arguments:
            variant: BgenVar, to be directly copied from one begn file to 
                another. This can be done when the new bgen file is for the same
                set of samples as the one being read from. This is much faster
                due to not having to decode and re-encode the genotype data.
      add_variant(varid, rsid, chrom, pos, alleles, genotypes, ploidy=2, 
                  phased=False, bit_depth=8)
        Arguments:
            varid: variant ID e.g. 'var1'
            rsid: reference SNP ID e.g. 'rs1'
            chrom: chromosome the variant is on e.g 'chr1'
            pos: nucleotide position of the variant e.g. 100
            alleles: list of allele strings e.g. ['A', 'C']
            genotypes: numpy array of genotype probabilities, ordered as per the
                bgen samples e.g. np.array([[0, 0, 1], [0.5, 0.5, 0]])
            ploidy: ploidy state, either as integer to indicate constant ploidy
                (e.g. 2), or numpy array of ploidy values per sample, e.g. np.array([1, 2, 2])
            phased: whether the genotypes are for phased data or not (default=False)
            bit_depth: how many bits to store each genotype as (1-32, default=8)

```
