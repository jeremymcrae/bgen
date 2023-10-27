

from pathlib import Path
import unittest
import tempfile
import os
import time
import sys

import numpy as np

from bgen import BgenReader, BgenWriter
from bgen.index import Index

def probs_close(orig, updat, bit_depth):
    ''' check if the genotype probabilities are near the original values

    We lose some precision when storing genotype proabilites, dependent on the
    bit depth we select. This accounts for the variable bit depth.
    '''
    # fix the nans from differing ploidies
    orig = orig.copy()
    orig[np.isnan(orig)] = 1
    updat[np.isnan(updat)] = 1

    # check abolute error (accounting for bit depth)
    max_error = 1 / ((2 ** bit_depth) - 1)
    matched = (abs(orig - updat) < max_error)

    # get relative difference between original and stored values
    delta = abs(updat / orig)
    delta[delta < 1] = 1 / delta[delta < 1]
    delta -= 1
    delta[np.isinf(delta)] = 0

    int_vals = np.floor((2 ** bit_depth - 1) * orig)
    max_delta = 1.0001 / int_vals

    # allow it through if the abolute error is sufficiently low, or the relative
    # error is sufficiently low, of if they differ by less than 1 part in 10 million
    matched = matched | (delta <= max_delta) | (delta < 1e-7)
    return matched.all()

class TestBgenWriter(unittest.TestCase):
    ''' class to make sure bgen.BgenWriter works correctly
    '''

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.tmpdir = Path(self.tmp.name)
    
    def tearDown(self):
        try:
            self.tmp.cleanup()
        except:
            pass
    
    def test_writing(self):
        ''' test basic BgenWriter file without variants
        '''
        path = self.tmpdir / 'temp.bgen'

        with BgenWriter(path, n_samples=3, samples=['a', 'b', 'c']) as bfile:
            pass
        
        with BgenReader(path, delay_parsing=True) as bfile:
            # check all the header attributes look ok
            self.assertEqual(bfile.samples, ['a', 'b', 'c'])
            self.assertEqual(bfile.header.offset, 37)
            self.assertEqual(bfile.header.nsamples, 3)
            self.assertEqual(bfile.header.nvariants, 0)
            self.assertEqual(bfile.header.compression, 'zstd')
            self.assertEqual(bfile.header.layout, 2)
            self.assertTrue(bfile.header.has_sample_ids)
            self.assertEqual(bfile.header.metadata, '')

        path2 = self.tmpdir / 'temp2.bgen'

        # check if we change the attributes, then we get the right data
        # bfile = BgenWriter(tmp.name, n_samples=4, samples=['a', 'b', 'c', 'd'], 
        with BgenWriter(path2, n_samples=4, samples=[], 
                           compression=None, layout=1, metadata='1234') as bfile:
            pass
        
        with BgenReader(path2, delay_parsing=True) as bfile:
            # check the new header attributes look ok
            self.assertEqual(bfile.samples, ['0', '1', '2', '3'])
            self.assertEqual(bfile.header.offset, 24)
            self.assertEqual(bfile.header.nsamples, 4)
            self.assertEqual(bfile.header.nvariants, 0)
            self.assertIsNone(bfile.header.compression)
            self.assertEqual(bfile.header.layout, 1)
            self.assertFalse(bfile.header.has_sample_ids)
            self.assertEqual(bfile.header.metadata, '1234')

    def test_writing_to_closed_file(self):
        ''' check we can't write to a closed file
        '''
        path = self.tmpdir / 'temp.bgen'
        sample_ids = ['a', 'b', 'c']
        with BgenWriter(path, n_samples=3, samples=sample_ids) as bfile:
            geno = np.array([[0.1, 0.8, 0.1], [0.5, 0.25, 0.25], [0.1, 0.2, 0.7]])
            bfile.add_variant('var1', 'rs1', 'chr1', 10, ['A', 'C'], geno,)

        with self.assertRaises(ValueError):
            bfile.add_variant('var1', 'rs1', 'chr1', 10, ['A', 'C'], geno,)
    
    def test_wrong_sample_number(self):
        ''' check we can't write variants with the wrong number of samples
        '''
        path = self.tmpdir / 'temp.bgen'
        sample_ids = ['a', 'b', 'c']
        with BgenWriter(path, n_samples=3, samples=sample_ids) as bfile:
            geno = np.array([[0.1, 0.8, 0.1], [0.5, 0.25, 0.25]])

            with self.assertRaises(ValueError):
                bfile.add_variant('var1', 'rs1', 'chr1', 10, ['A', 'C'], geno,)

    def test_writing_variant_attributes(self):
        ''' check we write variant attributes correctly (Aside from genotype)
        '''
        path = self.tmpdir / 'temp.bgen'
        sample_ids = ['a', 'b', 'c']
        with BgenWriter(path, n_samples=3, samples=sample_ids) as bfile:
            geno = np.array([[0.1, 0.8, 0.1], [0.5, 0.25, 0.25], [0.1, 0.2, 0.7]])
            bfile.add_variant('var1', 'rs1', 'chr1', 10, ['A', 'C'], geno,)
            bfile.add_variant('var2', 'rs2', 'chr1', 11, ['C', 'T'], geno / 2)
        
        with BgenReader(path, delay_parsing=True) as bfile:
            self.assertEqual(bfile.samples, sample_ids)

            # check that we've written all the attributes correctly
            var = next(bfile)
            self.assertEqual(var.varid, 'var1')
            self.assertEqual(var.rsid, 'rs1')
            self.assertEqual(var.chrom, 'chr1')
            self.assertEqual(var.pos, 10)
            self.assertEqual(var.alleles, ['A', 'C'])

            # and check the second variant
            var = next(bfile)
            self.assertEqual(var.varid, 'var2')
            self.assertEqual(var.rsid, 'rs2')
            self.assertEqual(var.chrom, 'chr1')
            self.assertEqual(var.pos, 11)
            self.assertEqual(var.alleles, ['C', 'T'])

    def test_writing_genotypes(self):
        ''' test BgenWriter
        '''
        path = self.tmpdir / 'temp.bgen'
        with BgenWriter(path, n_samples=3, samples=['a', 'b', 'c'],
                           compression=None, layout=2, metadata='1234') as bfile:
            geno = np.array([[0.1, 0.8, 0.1], 
                            [0.5, 0.25, 0.25], 
                            [float('nan'), float('nan'), float('nan')],
                            ])
            bit_depth = 16
            bfile.add_variant('var1', 'rs1', 'chr1', 10, ['A', 'C'], geno, bit_depth=bit_depth)
            bfile.add_variant('var2', 'rs2', 'chr1', 11, ['C', 'T'], geno / 2, bit_depth=bit_depth)
        
        with BgenReader(path, delay_parsing=True) as bfile:
            self.assertEqual(bfile.samples, ['a', 'b', 'c'])

            # check the first variant works
            # check the genotype probabilities are correct, but last sample is all nan
            var = next(bfile)
            self.assertTrue(np.isnan(var.probabilities[2, :]).all())
            self.assertTrue(probs_close(geno[:, :-1], var.probabilities[:, :-1], bit_depth))
            
            # and check the second variant works
            var = next(bfile)
            # adjust the genotype probabilities for the second variant
            geno = geno / 2
            geno[:, -1] = 1 - geno[:, :-1].sum(axis=1)
            self.assertTrue(probs_close(geno[:, :-1], var.probabilities[:, :-1], bit_depth))

            # we hit the end after two variants
            with self.assertRaises(StopIteration):
                next(bfile)

    def test_writing_many_genotypes(self):
        ''' test BgenWriter with a few thousand samples
        '''
        path = self.tmpdir / 'temp.bgen'
        n_samples = 2000
        a = np.linspace(0, 0.3, n_samples)
        b = np.linspace(0.3, 0.6, n_samples)
        geno = np.vstack([a, b, 1 - (a + b)]).T
        with BgenWriter(path, n_samples=n_samples, layout=1, compression=None) as bfile:
            bfile.add_variant('var1', 'rs1', 'chr1', 10, ['A', 'C'], geno)
        
        with BgenReader(path, delay_parsing=True) as bfile:
            var = next(bfile)
            self.assertTrue(probs_close(geno[:, :-1], var.probabilities[:, :-1], bit_depth=8))

    def test_compression_and_layouts(self):
        compressions = [None, 'zlib', 'zstd']
        layouts = [1, 2]
        geno = np.array([[0.1, 0.8, 0.1], 
                        [0.5, 0.25, 0.25], 
                        [float('nan'), float('nan'), float('nan')],
                        ])

        idx = 1
        for compression in compressions:
            for layout in layouts:
                if compression == 'zstd' and layout == 1:
                    continue
                path = path = self.tmpdir / f'temp_{idx}.bgen'
                idx += 1
                with BgenWriter(path, 3, samples=['a', 'b', 'c'],
                                compression=compression, layout=layout) as bfile:
                    bfile.add_variant('var1', 'rs1', 'chr1', 10, ['A', 'C'], geno)

                with BgenReader(path, delay_parsing=True) as bfile:
                    self.assertEqual(bfile.header.compression, compression)
                    self.assertEqual(bfile.header.layout, layout)
                    for x in bfile:
                        probs = x.probabilities
                        self.assertTrue(probs_close(geno[:, :-1], probs[:, :-1], 8))

    def test_bit_depths(self):
        ''' check writing to different bit depths works
        '''
        geno = np.array([[0.1, 0.8, 0.1], 
                        [0.5, 0.25, 0.25], 
                        [float('nan'), float('nan'), float('nan')],
                        ])
        
        idx = 1
        for bit_depth in range(1, 33):
            path = path = self.tmpdir / f'temp_{idx}.bgen'
            idx += 1
            with BgenWriter(path, 3, samples=['a', 'b', 'c']) as bfile:
                bfile.add_variant('var1', 'rs1', 'chr1', 10, ['A', 'C'], geno,
                                bit_depth=bit_depth)

            with BgenReader(path, delay_parsing=True) as bfile:
                for x in bfile:
                    probs = x.probabilities
                    self.assertTrue(probs_close(geno[:, :-1], probs[:, :-1], bit_depth))
    
    def test_more_alleles(self):
        ''' check writing to different bit depths works
        '''
        path = self.tmpdir / 'temp.bgen'
        geno1 = np.array([[0.1, 0.8, 0.1],
                        [0.5, 0.25, 0.25],
                        [float('nan'), float('nan'), float('nan')],
                        ])
        geno2 = np.array([[0.1, 0.6, 0, 0, 0.1, 0.2],
                        [0.1, 0.2, 0.1, 0.2, 0.1, 0.3],
                        [float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan')],
                        ])
        bit_depth = 8
        bfile = BgenWriter(path, 3, samples=['a', 'b', 'c'])
        bfile.add_variant('var1', 'rs1', 'chr1', 10, ['A', 'C'], geno1,
                        bit_depth=bit_depth)
        bfile.add_variant('var1', 'rs1', 'chr1', 10, ['A', 'C', 'T'], geno2,
                        bit_depth=bit_depth)
        bfile.close()

        bfile = BgenReader(path, delay_parsing=True)
        x = next(bfile)
        self.assertTrue(probs_close(geno1[:, :-1], x.probabilities[:, :-1], bit_depth))

        x = next(bfile)
        self.assertEqual(len(x.alleles), 3)
        self.assertTrue(probs_close(geno2[:, :-1], x.probabilities[:, :-1], bit_depth))

    def test_phased_data(self):
        '''checking writing phased data'''
        path = self.tmpdir / 'temp.bgen'
        bfile = BgenWriter(path, 3, samples=['a', 'b', 'c'])
        geno = np.array([[0.1, 0.9, 0.5, 0.5], 
                        [0.2, 0.8, 0.4, 0.6],
                        [float('nan'), float('nan'), float('nan'), float('nan')]])
        bfile.add_variant('var1', 'rs1', 'chr1', 10, ['A', 'C'], geno,
                            phased=1, bit_depth=8)
        bfile.close()

        bfile = BgenReader(path, delay_parsing=True)
        for x in bfile:
            probs = x.probabilities
            self.assertTrue(probs_close(geno[:, :-1], probs[:, :-1], bit_depth=8))
    
    def test_phased_data_many_samples(self):
        '''checking writing phased data with many samples.
        
        This also hits a fast path for parsing phased data'''
        path = self.tmpdir / 'temp.bgen'
        n_samples = 1000
        bfile = BgenWriter(path, n_samples)
        # construct a genotype array where the values
        a = np.linspace(0, 0.3, n_samples)
        b = np.linspace(0.7, 1, n_samples)
        geno = np.vstack([a, 1-a, b, 1-b]).T
        bfile.add_variant('var1', 'rs1', 'chr1', 10, ['A', 'C'], geno,
                            phased=1, bit_depth=8)
        bfile.close()

        bfile = BgenReader(path, delay_parsing=True)
        for x in bfile:
            probs = x.probabilities
            self.assertTrue(probs_close(geno[:, :-1], probs[:, :-1], bit_depth=8))
    
    def test_phased_data_different_sizes(self):
        '''checking writing phased data with a range of sample sizes'''
        path = self.tmpdir / 'temp.bgen'
        for n_samples in range(1, 100):
            bfile = BgenWriter(path, n_samples)
            # construct a genotype array where the values
            a = np.linspace(0, 0.3, n_samples)
            b = np.linspace(0.7, 1, n_samples)
            geno = np.vstack([a, 1-a, b, 1-b]).T
            bfile.add_variant('var1', 'rs1', 'chr1', 10, ['A', 'C'], geno,
                                phased=1, bit_depth=8)
            bfile.close()

            bfile = BgenReader(path, delay_parsing=True)
            for x in bfile:
                probs = x.probabilities
                self.assertTrue(probs_close(geno[:, :-1], probs[:, :-1], bit_depth=8))
    
    def test_ploidy_unphased(self):
        ''' check we can write unphased variants with variable ploidy per sample
        '''
        path = self.tmpdir / 'temp.bgen'
        bfile = BgenWriter(path, 3, samples=['a', 'b', 'c'])
        ploidy = np.array([1, 2, 3], dtype=np.uint8)
        geno = np.array([[0.1, 0.9, float('nan'), float('nan')], 
                        [0.2, 0.4, 0.4, float('nan')],
                        [float('nan'), float('nan'), float('nan'), float('nan')]])
        bfile.add_variant('var1', 'rs1', 'chr1', 10, ['A', 'C'], geno,
                            ploidy=ploidy)
        bfile.close()

        bfile = BgenReader(path, delay_parsing=True)
        for x in bfile:
            probs = x.probabilities
            self.assertTrue(probs_close(geno, probs, bit_depth=8))
        
    def test_ploidy_phased(self):
        ''' check we can write phased variants with variable ploidy per sample
        '''
        path = self.tmpdir / 'temp.bgen'
        bfile = BgenWriter(path, 4, samples=['a', 'b', 'c', 'd'])
        ploidy = np.array([1, 2, 3, 3], dtype=np.uint8)
        geno = np.array([[0.1, 0.9, float('nan'), float('nan'), float('nan'), float('nan')], 
                         [0.2, 0.8, 0.5, 0.5, float('nan'), float('nan')],
                         [float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan')],
                         [0.3, 0.7, 0.2, 0.8, 1, 0],
                         ])
        bfile.add_variant('var1', 'rs1', 'chr1', 10, ['A', 'C'], geno,
                            ploidy=ploidy, phased=1)
        bfile.close()

        bfile = BgenReader(path, delay_parsing=True)
        for x in bfile:
            probs = x.probabilities
            self.assertTrue(probs_close(geno, probs, bit_depth=8))
