

from pathlib import Path
import unittest

import numpy as np

from bgen import BgenReader

def get_alt_dose(variant):
    ''' compute alt dosage from variant probabilities
    
    This accounts for ploidy and phased status.
    '''
    probs = variant.probabilities
    haploid = variant.ploidy == 1
    
    alt_dose = np.empty(len(probs))
    if variant.is_phased:
        alt_dose[~haploid] = probs[~haploid, 1] + probs[~haploid, 3]
        alt_dose[haploid] = probs[haploid, 1]
    else:
        alt_dose[~haploid] = 2 * probs[~haploid, 2] + probs[~haploid, 1]
        alt_dose[haploid] = probs[haploid, 1]
    
    return alt_dose

class TestAltDosage(unittest.TestCase):
    ''' class to make sure alt dosage is computed correctly
    '''
    
    def setUp(self):
        ''' set path to folder with test data
        '''
        self.folder = Path(__file__).parent /  "data"
    
    def test_alt_dosage_nonstandard(self):
        ''' variant.alt_dosage is correct with variable ploidy and with phased data
        '''
        path = self.folder / 'alt_dosage_check.bgen'
        with BgenReader(path) as bfile:
            for variant in bfile:
                dose = variant.alt_dosage
                self.assertTrue((dose >= 0).all())
                self.assertTrue((dose == get_alt_dose(variant)).all())
