
from pathlib import Path
import unittest

import numpy as np

from bgen.reader import BgenFile

from tests.utils import load_gen_data, load_vcf_data, load_haps_data, arrays_equal

class TestExampleBgens(unittest.TestCase):
    ''' class to make sure we can load bgen files
    '''
    @classmethod
    def setUpClass(cls):
        cls.gen_data = load_gen_data()
        cls.vcf_data = load_vcf_data()
        cls.haps_data = load_haps_data()
    
    def setUp(self):
        ''' set path to folder with test data
        '''
        self.folder = Path(__file__).parent /  "data"
    
    def test_load_example_genotypes_bit_depths(self):
        ''' check parsing genotypes from the example files with different bit depths
        '''
        for path in self.folder.glob('example.*bits.bgen'):
            bit_depth = int(path.stem.split('.')[1].strip('bits'))
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
    
    def test_v11(self):
        ''' check we can open a bgen in v1.1 format, and parse genotypes correctly
        '''
        path = self.folder / 'example.v11.bgen'
        bfile = BgenFile(str(path))
        bit_depth = 16
        for var, g in zip(bfile, self.gen_data):
            self.assertEqual(g, var)
            self.assertTrue(arrays_equal(g.probabilities, var.probabilities, bit_depth))
    
    def test_Path(self):
        ''' check we can open bgen files from Path objects
        '''
        path = self.folder / 'example.v11.bgen'
        bfile = BgenFile(path)
    
    def test_load_haplotypes_bgen(self):
        ''' check we can open a bgen with haplotypes, and parse genotypes correctly
        '''
        path = self.folder / 'haplotypes.bgen'
        bfile = BgenFile(str(path))
        bit_depth = 16
        for var, g in zip(bfile, self.haps_data):
            self.assertEqual(g, var)
            self.assertTrue(arrays_equal(g.probabilities, var.probabilities, bit_depth))
    
    def test_load_complex_file(self):
        ''' make sure we can open a complex bgen file
        '''
        path = self.folder / 'complex.bgen'
        bfile = BgenFile(path)
        bit_depth = 16
        for var, g in zip(bfile, self.vcf_data):
            self.assertEqual(g, var)
            self.assertTrue(arrays_equal(g.probabilities, var.probabilities, bit_depth))
            self.assertTrue(all(x == y for x, y in zip(g.ploidy, var.ploidy)))
    
    def test_load_complex_files(self):
        ''' make sure we can open the complex bgen files
        '''
        
        for path in self.folder.glob('complex.*.bgen'):
            bit_depth = int(path.stem.split('.')[1].strip('bits'))
            bfile = BgenFile(path)
            for var, g in zip(bfile, self.vcf_data):
                self.assertEqual(g, var)
                self.assertTrue(arrays_equal(g.probabilities, var.probabilities, bit_depth))
    
    def test_load_missing_file(self):
        ''' check passing in a path to a missing file fails gracefully
        '''
        with self.assertRaises(ValueError):
            BgenFile('/zzz/jjj/qqq.bgen')
    
    # def test_load_missing_sample_file(self):
    #     path = str(self.folder / 'example.8bits.bgen')
    #     bfile = BgenFile(path, '/zzz/jjj/qqq.sample')
