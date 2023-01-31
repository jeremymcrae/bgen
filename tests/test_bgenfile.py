
from pathlib import Path
import unittest

import numpy as np

from bgen.reader import BgenFile

from tests.utils import load_gen_data

class TestBgenFile(unittest.TestCase):
    ''' class to make sure BgenFile works correctly
    '''
    
    @classmethod
    def setUpClass(cls):
        cls.gen_data = load_gen_data()
    
    def setUp(self):
        ''' set path to folder with test data
        '''
        self.folder = Path(__file__).parent /  "data"
    
    def test_context_handler_closed_bgen_samples(self):
        ''' no samples available from exited BgenFile
        '''
        path = self.folder / 'example.16bits.zstd.bgen'
        with BgenFile(path) as bfile:
            self.assertTrue(len(bfile.samples) > 0)
        
        with self.assertRaises(ValueError):
            bfile.samples
    
    def test_context_handler_closed_bgen_varids(self):
        ''' no varids available from exited BgenFile
        '''
        path = self.folder / 'example.16bits.zstd.bgen'
        with BgenFile(path) as bfile:
            self.assertTrue(len(bfile.varids()) > 0)
        
        with self.assertRaises(ValueError):
            bfile.varids()
    
    def test_context_handler_closed_bgen_rsids(self):
        ''' no rsids available from exited BgenFile
        '''
        path = self.folder / 'example.16bits.zstd.bgen'
        with BgenFile(path) as bfile:
            self.assertTrue(len(bfile.rsids()) > 0)
        
        with self.assertRaises(ValueError):
            bfile.rsids()
    
    def test_context_handler_closed_bgen_positions(self):
        ''' no positions available from exited BgenFile
        '''
        path = self.folder / 'example.16bits.zstd.bgen'
        with BgenFile(path) as bfile:
            self.assertTrue(len(bfile.positions()) > 0)
        
        with self.assertRaises(ValueError):
            bfile.positions()
    
    def test_context_handler_closed_bgen_length(self):
        ''' error raised if accessing length of exited BgenFile
        '''
        path = self.folder / 'example.16bits.zstd.bgen'
        with BgenFile(path) as bfile:
            self.assertTrue(len(bfile) > 0)
        
        with self.assertRaises(ValueError):
             len(bfile)
    
    def test_context_handler_closed_bgen_slice(self):
        ''' error raised if slicing variant from exited BgenFile
        '''
        path = self.folder / 'example.16bits.zstd.bgen'
        with BgenFile(path) as bfile:
            self.assertTrue(len(bfile) > 0)
        
        with self.assertRaises(ValueError):
             var = bfile[0]
    
    def test_context_handler_closed_bgen_at_position(self):
        ''' error raised if getting variant at position from exited BgenFile
        '''
        path = self.folder / 'example.16bits.zstd.bgen'
        with BgenFile(path) as bfile:
            self.assertTrue(len(bfile) > 0)
        
        with self.assertRaises(ValueError):
             var = bfile.at_position(100)
    
    def test_context_handler_closed_bgen_with_rsid(self):
        ''' error raised if getting variant with rsid from exited BgenFile
        '''
        path = self.folder / 'example.16bits.zstd.bgen'
        with BgenFile(path) as bfile:
            self.assertTrue(len(bfile) > 0)
        
        with self.assertRaises(ValueError):
             var = bfile.with_rsid('rs111')
    
    def test_fetch(self):
        ''' can fetch variants within a genomic region
        '''
        chrom, start, stop = '01', 5000, 50000
        bfile = BgenFile(self.folder / 'example.16bits.bgen')
        self.assertTrue(bfile._check_for_index(str(self.folder / 'example.16bits.bgen')))
        
        self.assertTrue(list(bfile.fetch('02')) == [])
    
    def test_fetch_whole_chrom(self):
        ''' fetching just with chrom gives all variants on chromosome
        '''
        chrom, start, stop = '01', 5000, 50000
        bfile = BgenFile(self.folder / 'example.16bits.bgen')
        
        # test fetching a whole chromosome
        sortkey = lambda x: (x.chrom, x.pos)
        for x, y in zip(sorted(bfile.fetch(chrom), key=sortkey), sorted(self.gen_data, key=sortkey)):
            self.assertEqual(x.rsid, y.rsid)
            self.assertEqual(x.chrom, y.chrom)
            self.assertEqual(x.pos, y.pos)
    
    def test_fetch_after_position(self):
        ''' fetching variants with chrom and start gives all variants after pos
        '''
        chrom, start, stop = '01', 5000, 50000
        print('testing fetch after position')
        bfile = BgenFile(self.folder / 'example.16bits.bgen')
        
        print('opened gen')
        sortkey = lambda x: (x.chrom, x.pos)
        gen_vars = [x for x in sorted(self.gen_data, key=sortkey) if start <= x.pos]
        print('have expected data')
        for x, y in zip(sorted(bfile.fetch(chrom, start), key=sortkey), gen_vars):
            print(x, y)
            self.assertEqual(x.rsid, y.rsid)
            self.assertEqual(x.chrom, y.chrom)
            self.assertEqual(x.pos, y.pos)
    
    def test_fetch_in_region(self):
        ''' fetching variants with chrom, start, stop gives variants in region
        '''
        chrom, start, stop = '01', 5000, 50000
        bfile = BgenFile(self.folder / 'example.16bits.bgen')
        
        sortkey = lambda x: (x.chrom, x.pos)
        gen_vars = [x for x in sorted(self.gen_data, key=sortkey) if start <= x.pos <= stop]
        for x, y in zip(sorted(bfile.fetch(chrom, start, stop), key=sortkey), gen_vars):
            self.assertEqual(x.rsid, y.rsid)
            self.assertEqual(x.chrom, y.chrom)
            self.assertEqual(x.pos, y.pos)
        
        # check that we don't get any variants in a region without any
        self.assertEqual(list(bfile.fetch(chrom, start * 1000, stop * 1000)), [])
