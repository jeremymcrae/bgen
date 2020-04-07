
from pathlib import Path
import unittest

import numpy as np

from bgen.reader import BgenFile

from tests.utils import load_gen_data, arrays_equal

class TestExampleBgens(unittest.TestCase):
    ''' class to make sure we can load bgen files
    '''
    @classmethod
    def setUpClass(cls):
        cls.gen_data = load_gen_data()
    
    def setUp(self):
        ''' set path to folder with test data
        '''
        self.folder = Path(__file__).parent /  "data"
    
    def test_load_example_files(self):
        ''' check we can open the example files
        '''
        for path in self.folder.glob('example*.bgen'):
            bfile = BgenFile(str(path))
    
    def test_load_example_genotypes(self):
        ''' check parsing genotypes from the example files
        '''
        for path in self.folder.glob('example*.bgen'):
            print(f'testing {path}')
            try:
                bit_depth = int(path.stem.split('.')[1].strip('bits'))
            except ValueError:
                bit_depth = 16
            bfile = BgenFile(str(path))
            for var, g in zip(bfile, self.gen_data):
                self.assertEqual(g, var)
                self.assertTrue(arrays_equal(g.probabilities, var.probabilities, bit_depth))
    
    def test_zstd_compressed(self):
        ''' check we can parse genotypes from zstd compressed geno probabilities
        '''
        path = self.folder / 'example.16bits.zstd.bgen'
        bfile = BgenFile(str(path))
        for var, g in zip(bfile, self.gen_data):
            self.assertEqual(g, var)
            self.assertTrue(arrays_equal(g.probabilities, var.probabilities, 16))
    
    def test_load_complex_files(self):
        ''' make sure we can open the complex bgen files
        '''
        for path in self.folder.glob('complex*.bgen'):
            # bfile = BgenFile(str(path))
            pass
    
    def test_load_missing_file(self):
        ''' check passing in a path to a missing file fails gracefully
        '''
        with self.assertRaises(ValueError):
            BgenFile('/zzz/jjj/qqq.bgen')
    
    # def test_load_missing_sample_file(self):
    #     path = str(self.folder / 'example.8bits.bgen')
    #     bfile = BgenFile(path, '/zzz/jjj/qqq.sample')
