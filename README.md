### Another bgen reader
![travis](https://travis-ci.org/jeremymcrae/bgen.svg?branch=master)

This is a package for reading [bgen files](https://www.well.ox.ac.uk/~gav/).

This package uses cython to wrap c++ code for parsing bgen files. It's not too
slow, it can parse genotypes from 500,000 individuals at >100 variants per
second within python.

This has been primarily been designed around UKBiobank bgen files (i.e. bgen
version 1.2 with zlib compressed genotype probabilities, but the other versions
and compression schemes have also been tested using example bgen files).

#### Install
`pip install bgen`

#### Usage
```python
from bgen import BgenFile
bfile = BgenFile(BGEN_PATH, sample_path=None)
rsids = bfile.rsids()

# select a variant by indexing
var = bfile[1000]

# pull out genotype probabilities
probs = var.probabilities  # returns 2D numpy array
dosage = var.minor_allele_dosage  # requires biallelic variant, returns numpy array

# exclude variants from analyses by passing in indices
to_drop = [1, 3, 500]
bfile.drop_variants(to_drop)

# pickle variants for easy message passing
import pickle
dumped = pickle.dumps(var)
var = pickle.loads(dumped)

# iterate through every variant in the file, without preloading every variant
with BgenFile(BGEN_PATH, sample_path=None, delay_parsing=True) as bfile:
  for var in bfile:
      probs = var.probabilities
      dosage = var.minor_allele_dosage
      ploidy = var.ploidy
```
