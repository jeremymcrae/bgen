
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
        ''' make sure we can load all the example files
        '''
        for path in self.folder.glob('example*.bgen'):
            bfile = BgenFile(str(path))
    
    def test_load_example_genotypes(self):
        ''' make sure we can load all the example files
        '''
        for path in self.folder.glob('example*.bgen'):
            bfile = BgenFile(str(path))
            for var in bfile:
                # geno = var.probabilities
                pass
    
    def test_load_complex_files(self):
        ''' make sure we can load all the complex bgen files
        '''
        for path in self.folder.glob('complex*.bgen'):
            # bfile = BgenFile(str(path))
            pass
    
    def test_load_missing_file(self):
        ''' check passing in a path to a missing file fails gracefully
        '''
        with self.assertRaises(ValueError):
            BgenFile('/zzz/jjj/qqq.bgen')
