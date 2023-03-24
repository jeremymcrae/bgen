

from pathlib import Path
import unittest
import tempfile
import os

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
    max_delta = 1 / int_vals

    # allow it through if the abolute error is sufficiently low, or the relative
    # error is sufficiently low, of if they differ by less than 1 part in 10 million
    matched = matched | (delta <= max_delta) | (delta < 1e-7)
    return matched.all()

class TestBgenWriter(unittest.TestCase):
    ''' class to make sure bgen.BgenWriter works correctly
    '''

    def setUp(self):
        self.tmp = tempfile.NamedTemporaryFile(delete=False)
        self.tmp.close()
        self.path = self.tmp.name
    
    def tearDown(self):
        os.unlink(self.path)
    
    def test_writing(self):
        ''' test basic BgenWriter file without variants
        '''
        bfile = BgenWriter(self.path, n_samples=3, samples=['a', 'b', 'c'])
        bfile = BgenReader(self.path, delay_parsing=True)

        # check all the header attributes look ok
        self.assertEqual(bfile.samples, ['a', 'b', 'c'])
        self.assertEqual(bfile.header.offset, 37)
        self.assertEqual(bfile.header.nsamples, 3)
        self.assertEqual(bfile.header.nvariants, 0)
        self.assertEqual(bfile.header.compression, 'zstd')
        self.assertEqual(bfile.header.layout, 2)
        self.assertTrue(bfile.header.has_sample_ids)
        self.assertEqual(bfile.header.metadata, '')

        # check if we change the attributes, then we get the right data
        # bfile = BgenWriter(self.path, n_samples=4, samples=['a', 'b', 'c', 'd'], 
        bfile = BgenWriter(self.path, n_samples=4, samples=[], 
                           compression=None, layout=1, metadata='1234')
        bfile.close()
        bfile = BgenReader(self.path, delay_parsing=True)

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
        sample_ids = ['a', 'b', 'c']
        bfile = BgenWriter(self.path, n_samples=3, samples=sample_ids)
        geno = np.array([[0.1, 0.8, 0.1], [0.5, 0.25, 0.25], [0.1, 0.2, 0.7]])
        bfile.add_variant('var1', 'rs1', 'chr1', 10, ['A', 'C'], 3, geno,)
        bfile.close()

        with self.assertRaises(ValueError):
            bfile.add_variant('var1', 'rs1', 'chr1', 10, ['A', 'C'], 3, geno,)
    
    def test_wrong_sample_number(self):
        ''' check we can't write variants with the wrong number of samples
        '''
        sample_ids = ['a', 'b', 'c']
        bfile = BgenWriter(self.path, n_samples=3, samples=sample_ids)
        geno = np.array([[0.1, 0.8, 0.1], [0.5, 0.25, 0.25]])

        with self.assertRaises(ValueError):
            bfile.add_variant('var1', 'rs1', 'chr1', 10, ['A', 'C'], 3, geno,)

    def test_writing_variant_attributes(self):
        ''' check we write variant attributes correctly (Aside from genotype)
        '''
        sample_ids = ['a', 'b', 'c']
        bfile = BgenWriter(self.path, n_samples=3, samples=sample_ids)
        geno = np.array([[0.1, 0.8, 0.1], [0.5, 0.25, 0.25], [0.1, 0.2, 0.7]])
        bfile.add_variant('var1', 'rs1', 'chr1', 10, ['A', 'C'], 3, geno,)
        bfile.add_variant('var2', 'rs2', 'chr1', 11, ['C', 'T'], 3, geno / 2)
        
        bfile = BgenReader(self.path, delay_parsing=True)
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
        bfile = BgenWriter(self.path, n_samples=3, samples=['a', 'b', 'c'],
                           compression=None, layout=2, metadata='1234')

        geno = np.array([[0.1, 0.8, 0.1], 
                        [0.5, 0.25, 0.25], 
                        [float('nan'), float('nan'), float('nan')],
                        ])
        bit_depth = 16
        bfile.add_variant('var1', 'rs1', 'chr1', 10, ['A', 'C'], 3, geno, bit_depth=bit_depth)
        bfile.add_variant('var2', 'rs2', 'chr1', 11, ['C', 'T'], 3, geno / 2, bit_depth=bit_depth)
        
        bfile = BgenReader(self.path, delay_parsing=True)
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
        bfile._close()


    def test_compression_and_layouts(self):
        compressions = [None, 'zlib', 'zstd']
        layouts = [1, 2]
        geno = np.array([[0.1, 0.8, 0.1], 
                        [0.5, 0.25, 0.25], 
                        [float('nan'), float('nan'), float('nan')],
                        ])

        for compression in compressions:
            for layout in layouts:
                if compression == 'zstd' and layout == 1:
                    continue
                bfile = BgenWriter(self.path, 3, samples=['a', 'b', 'c'],
                                compression=compression, layout=layout)
                bfile.add_variant('var1', 'rs1', 'chr1', 10, ['A', 'C'], 3, geno)
                bfile.close()

                bfile = BgenReader(self.path, delay_parsing=True)
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
        for bit_depth in range(1, 33):
            bfile = BgenWriter(self.path, 3, samples=['a', 'b', 'c'])
            bfile.add_variant('var1', 'rs1', 'chr1', 10, ['A', 'C'], 3, geno,
                            bit_depth=bit_depth)
            bfile.close()

            bfile = BgenReader(self.path, delay_parsing=True)
            for x in bfile:
                probs = x.probabilities
                self.assertTrue(probs_close(geno[:, :-1], probs[:, :-1], bit_depth))
    
    def test_more_alleles(self):
        ''' check writing to different bit depths works
        '''
        geno1 = np.array([[0.1, 0.8, 0.1],
                        [0.5, 0.25, 0.25],
                        [float('nan'), float('nan'), float('nan')],
                        ])
        geno2 = np.array([[0.1, 0.6, 0, 0, 0.1, 0.2],
                        [0.1, 0.2, 0.1, 0.2, 0.1, 0.3],
                        [float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan')],
                        ])
        bit_depth = 8
        bfile = BgenWriter(self.path, 3, samples=['a', 'b', 'c'])
        bfile.add_variant('var1', 'rs1', 'chr1', 10, ['A', 'C'], 3, geno1,
                        bit_depth=bit_depth)
        bfile.add_variant('var1', 'rs1', 'chr1', 10, ['A', 'C', 'T'], 3, geno2,
                        bit_depth=bit_depth)
        bfile.close()

        bfile = BgenReader(self.path, delay_parsing=True)
        x = next(bfile)
        self.assertTrue(probs_close(geno1[:, :-1], x.probabilities[:, :-1], bit_depth))

        x = next(bfile)
        self.assertEqual(len(x.alleles), 3)
        self.assertTrue(probs_close(geno2[:, :-1], x.probabilities[:, :-1], bit_depth))

    def test_phased_data(self):
        '''checking writing phased data'''
        bfile = BgenWriter(self.path, 3, samples=['a', 'b', 'c'])
        geno = np.array([[0.1, 0.9, 0.5, 0.5], 
                        [0.2, 0.8, 0.4, 0.6],
                        [float('nan'), float('nan'), float('nan'), float('nan')]])
        bfile.add_variant('var1', 'rs1', 'chr1', 10, ['A', 'C'], 3, geno, 
                            phased=1, bit_depth=8)
        bfile.close()

        bfile = BgenReader(self.path, delay_parsing=True)
        for x in bfile:
            probs = x.probabilities
            self.assertTrue(probs_close(geno[:, :-1], probs[:, :-1], bit_depth=8))
    
    def test_ploidy_unphased(self):
        ''' check we can write unphased variants with variable ploidy per sample
        '''
        bfile = BgenWriter(self.path, 3, samples=['a', 'b', 'c'])
        ploidy = np.array([1, 2, 3], dtype=np.uint8)
        geno = np.array([[0.1, 0.9, float('nan'), float('nan')], 
                        [0.2, 0.4, 0.4, float('nan')],
                        [float('nan'), float('nan'), float('nan'), float('nan')]])
        bfile.add_variant('var1', 'rs1', 'chr1', 10, ['A', 'C'], 3, geno, 
                            ploidy=ploidy)
        bfile.close()

        bfile = BgenReader(self.path, delay_parsing=True)
        for x in bfile:
            probs = x.probabilities
            self.assertTrue(probs_close(geno, probs, bit_depth=8))
        
    def test_ploidy_phased(self):
        ''' check we can write phased variants with variable ploidy per sample
        '''
        bfile = BgenWriter(self.path, 4, samples=['a', 'b', 'c', 'd'])
        ploidy = np.array([1, 2, 3, 3], dtype=np.uint8)
        geno = np.array([[0.1, 0.9, float('nan'), float('nan'), float('nan'), float('nan')], 
                         [0.2, 0.8, 0.5, 0.5, float('nan'), float('nan')],
                         [float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan')],
                         [0.3, 0.7, 0.2, 0.8, 1, 0],
                         ])
        bfile.add_variant('var1', 'rs1', 'chr1', 10, ['A', 'C'], 4, geno, 
                            ploidy=ploidy, phased=1)
        bfile.close()

        bfile = BgenReader(self.path, delay_parsing=True)
        for x in bfile:
            probs = x.probabilities
            self.assertTrue(probs_close(geno, probs, bit_depth=8))
