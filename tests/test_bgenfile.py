
from pathlib import Path
import unittest

import numpy as np

from bgen.reader import BgenFile

class TestBgenFile(unittest.TestCase):
    ''' class to make sure BgenFile works correctly
    '''
    
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
