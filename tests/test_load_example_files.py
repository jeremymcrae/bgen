
from pathlib import Path
import unittest
import tempfile

import numpy as np

from bgen import BgenReader, BgenWriter

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
            bfile = BgenReader(str(path))
            for var, g in zip(bfile, self.gen_data):
                self.assertEqual(g, var)
                self.assertTrue(arrays_equal(g.probabilities, var.probabilities, bit_depth))
    
    def test_zstd_compressed(self):
        ''' check we can parse genotypes from zstd compressed geno probabilities
        '''
        path = self.folder / 'example.16bits.zstd.bgen'
        bfile = BgenReader(str(path))
        for var, g in zip(bfile, self.gen_data):
            self.assertEqual(g, var)
            self.assertTrue(arrays_equal(g.probabilities, var.probabilities, 16))
    
    def test_v11(self):
        ''' check we can open a bgen in v1.1 format, and parse genotypes correctly
        '''
        path = self.folder / 'example.v11.bgen'
        bfile = BgenReader(str(path))
        bit_depth = 16
        for var, g in zip(bfile, self.gen_data):
            self.assertEqual(g, var)
            self.assertTrue(arrays_equal(g.probabilities, var.probabilities, bit_depth))
    
    def test_Path(self):
        ''' check we can open bgen files from Path objects
        '''
        path = self.folder / 'example.v11.bgen'
        bfile = BgenReader(path)
    
    def test_load_haplotypes_bgen(self):
        ''' check we can open a bgen with haplotypes, and parse genotypes correctly
        '''
        path = self.folder / 'haplotypes.bgen'
        bfile = BgenReader(str(path))
        bit_depth = 16
        for var, g in zip(bfile, self.haps_data):
            self.assertEqual(g, var)
            self.assertTrue(arrays_equal(g.probabilities, var.probabilities, bit_depth))
    
    def test_load_complex_file(self):
        ''' make sure we can open a complex bgen file
        '''
        path = self.folder / 'complex.bgen'
        bfile = BgenReader(path)
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
            bfile = BgenReader(path)
            for var, g in zip(bfile, self.vcf_data):
                self.assertEqual(g, var)
                self.assertTrue(arrays_equal(g.probabilities, var.probabilities, bit_depth))
    
    def test_load_sample_file(self):
        ''' make sure we can open a sample file for a bgen
        '''
        bgen_path = self.folder / 'complex.bgen'
        sample_path = self.folder / 'complex.sample'
        orig = BgenReader(bgen_path)
        orig_samples = orig.samples
        
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            # construct a bgen file without any sample IDs inside (otherwise 
            # the internal IDs take priority).
            bgen_path = tmp / 'temp.bgen'
            with BgenWriter(bgen_path, n_samples=len(orig.samples)) as bfile:
                for var in orig:
                    bfile.add_variant_direct(var)
            
            # check bgen files without internal IDs or an external file instead
            # use numeric IDs (converted to strings)
            bfile = BgenReader(bgen_path)
            numeric_ids = [f'{x}' for x in range(len(orig_samples))]
            self.assertEqual(numeric_ids, bfile.samples)
            bfile.close()
            
            # reading sample IDs from the corresponding sample file should give
            # identical IDs
            with BgenReader(bgen_path, sample_path) as bfile:
                self.assertEqual(orig_samples, bfile.samples)
            
            # check we raise an error with too few sample IDs
            missing_path = tmp / 'empty.sample'
            missing = open(missing_path, 'wt')
            missing.close()
            
            with self.assertRaises(ValueError):
                BgenReader(bgen_path, missing_path)
            
            # check we raise an error with too many sample IDs
            extra_path = tmp / 'empty.sample'
            with open(extra_path, 'wt') as extra:
                extra.write('id\n0\nsample_0\nsample_1\nsample_2\nsample_3\nsample_4\n')
            
            with self.assertRaises(ValueError):
                BgenReader(bgen_path, extra_path)
    
    def test_load_missing_file(self):
        ''' check passing in a path to a missing file fails gracefully
        '''
        with self.assertRaises(ValueError):
            BgenReader('/zzz/jjj/qqq.bgen')
    
    # def test_load_missing_sample_file(self):
    #     path = str(self.folder / 'example.8bits.bgen')
    #     bfile = BgenReader(path, '/zzz/jjj/qqq.sample')
