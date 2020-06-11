
from pathlib import Path
import unittest

import numpy as np

from bgen.reader import BgenFile
from bgen.index import Index

from tests.utils import load_gen_data

class TestBgenIndex(unittest.TestCase):
    ''' class to make sure bgen.index.Index works correctly
    '''
    @classmethod
    def setUpClass(cls):
        cls.gen_data = load_gen_data()
    
    def setUp(self):
        ''' set path to folder with test data
        '''
        self.folder = Path(__file__).parent /  "data"
    
    def test_index_opens(self):
        ''' loads index when available
        '''
        bfile = BgenFile(self.folder / 'example.15bits.bgen')
        self.assertFalse(bfile._check_for_index(str(self.folder / 'example.15bits.bgen')))
        
        bfile = BgenFile(self.folder / 'example.16bits.bgen')
        self.assertTrue(bfile._check_for_index(str(self.folder / 'example.16bits.bgen')))
    
    def test_index_fetch(self):
        ''' fetches file offsets
        '''
        chrom = '01'
        start = 5000
        stop = 50000
        
        index = Index(self.folder / 'example.16bits.bgen.bgi')
        self.assertTrue(len(list(index.fetch(chrom))) == len(self.gen_data))
        self.assertTrue(len(list(index.fetch('02'))) == 0)
        self.assertTrue(len(list(index.fetch(chrom, start * 100, stop * 100))) == 0)
        
        # check for a whole chromosome
        chrom_offsets = list(index.fetch(chrom))
        self.assertTrue(len(chrom_offsets) > 0)
        self.assertTrue(len(chrom_offsets) == len(self.gen_data))
        
        # check for all variants following a position
        after_pos_offsets = list(index.fetch(chrom, start))
        self.assertTrue(len(after_pos_offsets) > 0)
        self.assertTrue(len(after_pos_offsets) == len([x for x in self.gen_data if start <= x.pos]))
        
        # check for all variants within a region
        in_region_offsets = list(index.fetch(chrom, start, stop))
        self.assertTrue(len(in_region_offsets) > 0)
        self.assertTrue(len(in_region_offsets) == len([x for x in self.gen_data if start <= x.pos <= stop]))
        
        # make sure the queries return different lists
        self.assertTrue(len(chrom_offsets) != len(after_pos_offsets))
        self.assertTrue(len(chrom_offsets) != len(in_region_offsets))
        self.assertTrue(len(after_pos_offsets) != len(in_region_offsets))
        
