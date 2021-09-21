### Another bgen reader
![bgen](https://github.com/jeremymcrae/bgen/workflows/bgen/badge.svg)

This is a package for reading [bgen files](https://www.well.ox.ac.uk/~gav/bgen_format).

This package uses cython to wrap c++ code for parsing bgen files. It's fairly
quick, it can parse genotypes from 500,000 individuals at ~300 variants per
second within a single python process (~450 million probabilities per second
with a 3GHz CPU). Decompressing the genotype probabilities is the slow step,
zlib decompression takes 80% of the total time, using zstd compressed genotypes
would be much faster, maybe 2-3X faster?

This has been optimized for UKBiobank bgen files (i.e. bgen version 1.2 with
zlib compressed 8-bit genotype probabilities, but the other bgen versions and
zstd compression have also been tested using example bgen files).

#### Install
`pip install bgen`

#### Usage
```python
from bgen.reader import BgenFile

bfile = BgenFile(BGEN_PATH)
rsids = bfile.rsids()

# select a variant by indexing
var = bfile[1000]

# pull out genotype probabilities
probs = var.probabilities  # returns 2D numpy array
dosage = var.minor_allele_dosage  # returns 1D numpy array for biallelic variant

# iterate through every variant in the file
with BgenFile(BGEN_PATH, delay_parsing=True) as bfile:
  for var in bfile:
      dosage = var.minor_allele_dosage

# get all variants in a genomic region
variants = bfile.fetch('21', 10000, 5000000)
```

#### API documentation

``` py
class BgenFile(path, sample_path='', delay_parsing=False)
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
    with_rsid(pos): returns BgenVar with given position
    at_position(rsid): returns BgenVar with given rsid
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
    probabilitiies:  2D numpy array of genotype probabilities, one sample per row
  
  BgenVars can be pickled e.g. pickle.dumps(var)
```
