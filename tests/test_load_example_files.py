
from pathlib import Path
import unittest

from bgen.reader import BgenFile

class TestExampleBgens(unittest.TestCase):
    ''' class to make sure we can load bgen files
    '''
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
            bfile = BgenFile(str(path))
            for var in bfile:
                geno = var.probabilities
    
    def test_zstd_compressed(self):
        ''' check we can parse genotypes from zstd compressed geno probabilities
        '''
        path = self.folder / 'example.16bits.zstd.bgen'
        bfile = BgenFile(str(path))
        for var in bfile:
            geno = var.probabilities
    
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
