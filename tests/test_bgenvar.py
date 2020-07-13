
from pathlib import Path
import unittest
import pickle

import numpy as np

from bgen.reader import BgenFile

class TestBgenVar(unittest.TestCase):
    ''' class to make sure BgenVar works correctly
    '''
    
    def setUp(self):
        ''' set path to folder with test data
        '''
        self.folder = Path(__file__).parent /  "data"
      
    def test_minor_allele_dosage(self):
          ''' test we calculate minor_allele_dosage correctly
          '''
          path = self.folder / 'example.16bits.zstd.bgen'
          with BgenFile(path) as bfile:
              for var in bfile:
                  dose = var.minor_allele_dosage
                  probs = var.probabilities
                  
                  # calculate dosages for each allele
                  a1 = (probs[:, 0] * 2 + probs[:, 1])
                  a2 = (probs[:, 2] * 2 + probs[:, 1])
                  
                  # get delta between var.minor_allele_dosage and values calculated here
                  recomputed = a2 if np.nansum(a1) >= np.nansum(a2) else a1
                  delta = abs(dose - recomputed)
                  
                  # check difference between the two estimates is sufficiently low
                  self.assertTrue(np.nanmax(delta) < 2e-7)
    
    def test_minor_allele_dosage_fast(self):
        ''' test we calculate minor_allele_dosage correctly with the fast path
        '''
        path = self.folder / 'example.8bits.bgen'
        with BgenFile(path) as bfile:
            for var in bfile:
                dose = var.minor_allele_dosage
                probs = var.probabilities
                
                # calculate dosages for each allele
                a1 = (probs[:, 0] * 2 + probs[:, 1])
                a2 = (probs[:, 2] * 2 + probs[:, 1])
                
                # get delta between var.minor_allele_dosage and values calculated here
                recomputed = a2 if np.nansum(a1) >= np.nansum(a2) else a1
                delta = abs(dose - recomputed)
                
                # check difference between the two estimates is sufficiently low
                self.assertTrue(np.nanmax(delta) < 3e-7)
    
    def test_pickling(self):
        ''' BgenVar should pickle and unpickle
        '''
        path = self.folder / 'example.16bits.zstd.bgen'
        with BgenFile(path) as bfile:
            for var in bfile:
                # this checks that we can pickle and unpickle a BgenVar
                pickled = pickle.dumps(var)
                unpickled = pickle.loads(pickled)
                
                # check attributes of the original and unpickled are identical
                self.assertEqual(var.varid, unpickled.varid)
                self.assertEqual(var.rsid, unpickled.rsid)
                self.assertEqual(var.chrom, unpickled.chrom)
                self.assertEqual(var.pos, unpickled.pos)
                self.assertEqual(var.alleles, unpickled.alleles)
                
