
from pathlib import Path
import unittest

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
                  a2 = (probs[:, 2] * 2 + probs[0, 1])
                  
                  # get delta between var.minor_allele_dosage and values calculated here
                  if np.nansum(a1) >= np.nansum(a2):
                      delta = abs(dose - a1)
                  else:
                      delta = abs(dose - a1)
                  
                  # check difference between the two estimates is sufficiently low
                  self.assertTrue(np.nanmax(delta) < 1e-15)
