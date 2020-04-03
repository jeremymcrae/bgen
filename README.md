### Another bgen reader
This is a package for reading [bgen files](https://www.well.ox.ac.uk/~gav/).

This package uses cython to wrap c++ code for bgen parsing. It's fairly quick,
it can parse genotypes from 500,000 individuals at ~90 variants per second
within python.

#### Install
`pip install bgen`

#### Usage
```python
from bgen import BgenFile
bfile = BgenFile(BGEN_PATH, SAMPLE_PATH=None)
rsids = var.rsids()

# iterate through every variant in the file
for var in bfile:
    probs = var.probabilities()
    dosage = var.alt_dosage()  # requires biallelic variant, returns numpy array

# get one variant by indexing
var = bfile[1000]

# exclude variants from analyses by passing in indices
to_drop = [1, 3, 500]
bfile.drop_variants(to_drop)

# pickle variants for easy message passing
import pickle
dumped = pickle.dumps(var)
var = pickle.loads(dumped)
```
