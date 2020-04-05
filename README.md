### Another bgen reader
This is a package for reading [bgen files](https://www.well.ox.ac.uk/~gav/).

This package uses cython to wrap c++ code for parsing bgen files. It's not too
slow, it can parse genotypes from 500,000 individuals at >100 variants per
second within python.

This has been tested with UKBiobank bgen files (i.e. bgen version 1.2 with zlib
compressed genotype probabilities, but the other versions and compression
schemes should also work).

#### Install
`pip install bgen` (possibly needs `pip install cython` beforehand)

#### Usage
```python
from bgen import BgenFile
bfile = BgenFile(BGEN_PATH, SAMPLE_PATH=None)
rsids = var.rsids()

# iterate through every variant in the file
with BgenFile(BGEN_PATH, SAMPLE_PATH=None) as bfile:
  for var in bfile:
      probs = var.probabilities()  # returns 2D numpy array
      dosage = var.alt_dosage()  # requires biallelic variant, returns numpy array

# select a variant by indexing
var = bfile[1000]

# exclude variants from analyses by passing in indices
to_drop = [1, 3, 500]
bfile.drop_variants(to_drop)

# pickle variants for easy message passing
import pickle
dumped = pickle.dumps(var)
var = pickle.loads(dumped)
```
