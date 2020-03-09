### Yet another bgen reader

#### Install
pip install bgen

#### Usage
```python
from bgen import BgenFile
bfile = BgenFile(BGEN_PATH, SAMPLE_PATH=None)
rsids = var.rsids()

for var in bfile:
    probs = var.probabilities()
    dosage = var.alt_dosage()  # requires biallelic variant
```
