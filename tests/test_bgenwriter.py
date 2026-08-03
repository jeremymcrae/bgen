

from pathlib import Path
import unittest
import tempfile
import os
import time
import sys
import math

import numpy as np

from bgen import BgenReader, BgenWriter
from bgen.index import Index

def probs_close(orig, updat, bit_depth):
    ''' check if the genotype probabilities are near the original values

    We lose some precision when storing genotype probabilities, dependent on the
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
    
    def test_writing_metadata_with_newline(self):
        ''' check we can write metadata (and read back!)
        '''
        sample_ids = ['a', 'b', 'c']
        metadata = 'a\nbc'  # previously had errors with metadata containing newlines
        path = self.tmpdir / 'temp.bgen'
        
        bit_depth = 16
        geno = np.array([[0.1, 0.8, 0.1], [0.5, 0.25, 0.25], [0.1, 0.2, 0.7]])
        with BgenWriter(path, n_samples=3, samples=sample_ids, metadata=metadata) as bfile:
            bfile.add_variant('var1', 'rs1', 'chr1', 10, ['A', 'C'], geno, bit_depth=bit_depth)
        
        bfile = BgenReader(path)
        self.assertEqual(metadata, bfile.header.metadata)
        var = next(bfile)
        self.assertTrue(probs_close(geno[:, :-1], var.probabilities[:, :-1], bit_depth))
    
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
    
    def test_probabilities_stay_in_range(self):
        ''' check the inferred final probability is never out of range

        Layout 2 only stores all but the last probability per sample, and the
        reader infers the last as max - sum(stored). Rounding each probability
        independently let the stored values sum to more than max (e.g. 0.5 and
        0.5 both round to 128 at a bit depth of 8, but the max is 255), which
        made the reader infer a negative probability.
        '''
        geno = np.array([[0.5, 0.5, 0.0],
                         [0.0, 0.5, 0.5],
                         [0.25, 0.5, 0.25],
                         ])

        for bit_depth in range(1, 33):
            path = self.tmpdir / f'range_{bit_depth}.bgen'
            with BgenWriter(path, 3, samples=['a', 'b', 'c']) as bfile:
                bfile.add_variant('var1', 'rs1', 'chr1', 10, ['A', 'C'], geno,
                                  bit_depth=bit_depth)

            with BgenReader(path, delay_parsing=True) as bfile:
                for x in bfile:
                    probs = x.probabilities
                    # allow for float32 error in the reader's remainder
                    self.assertGreaterEqual(probs.min(), -1e-7,
                                            f'negative probability at bit_depth={bit_depth}')
                    self.assertLessEqual(probs.max(), 1 + 1e-7,
                                         f'probability above one at bit_depth={bit_depth}')
                    self.assertTrue(np.allclose(probs.sum(axis=1), 1.0, atol=1e-6),
                                    f'probabilities do not sum to one at bit_depth={bit_depth}')

    def test_probabilities_match_across_read_paths(self):
        ''' check the constant-ploidy fast path agrees with the generic path

        Reading 8-bit unphased diploid data uses a lookup table keyed on
        255 - first - second, which read outside the table when the stored
        values summed above 255, and so disagreed with the generic path.
        '''
        geno = np.array([[0.5, 0.5, 0.0],
                         [0.0, 0.5, 0.5],
                         [0.25, 0.5, 0.25],
                         ])

        fast_path = self.tmpdir / 'fast.bgen'
        with BgenWriter(fast_path, 3, samples=['a', 'b', 'c']) as bfile:
            bfile.add_variant('var1', 'rs1', 'chr1', 10, ['A', 'C'], geno,
                              ploidy=2, bit_depth=8)

        # a varying ploidy skips the fast path, but sample 0 and 1 stay diploid
        generic_path = self.tmpdir / 'generic.bgen'
        with BgenWriter(generic_path, 3, samples=['a', 'b', 'c']) as bfile:
            bfile.add_variant('var1', 'rs1', 'chr1', 10, ['A', 'C'], geno,
                              ploidy=np.array([2, 2, 1], dtype=np.uint8),
                              bit_depth=8)

        with BgenReader(fast_path, delay_parsing=True) as bfile:
            fast = bfile[0].probabilities
        with BgenReader(generic_path, delay_parsing=True) as bfile:
            generic = bfile[0].probabilities

        for i in range(2):
            self.assertTrue(np.allclose(fast[i], generic[i, :3], atol=1e-7),
                            f'sample {i}: {fast[i]} != {generic[i, :3]}')

    def test_reading_probabilities_above_maximum(self):
        ''' check a bgen storing probabilities summing above the maximum is safe

        Older versions of this package wrote such files, and nothing stops a
        third party from doing so. The reader used to index outside its lookup
        table for these, which read whatever was adjacent in memory.
        '''
        path = self.tmpdir / 'temp.bgen'
        geno = np.array([[1.0, 0.0, 0.0],
                         [1.0, 0.0, 0.0],
                         [1.0, 0.0, 0.0],
                         ])
        with BgenWriter(path, 3, samples=['a', 'b', 'c'],
                        compression=None) as bfile:
            bfile.add_variant('var1', 'rs1', 'chr1', 10, ['A', 'C'], geno,
                              bit_depth=8)

        # each sample stores (255, 0), so rewrite them to (255, 255), which sums
        # to 510 and makes the inferred final probability 255 - 510 = -255
        raw = bytearray(path.read_bytes())
        patched = 0
        for i in range(len(raw) - 1):
            if raw[i] == 255 and raw[i + 1] == 0:
                raw[i + 1] = 255
                patched += 1
        self.assertEqual(patched, 3)
        path.write_bytes(bytes(raw))

        with BgenReader(path, delay_parsing=True) as bfile:
            for x in bfile:
                probs = x.probabilities
                self.assertTrue(np.all((probs >= 0) & (probs <= 1)),
                                f'probabilities out of range: {probs}')
                # the inferred probability is clamped to zero rather than read
                # from outside the lookup table, which returned adjacent memory
                self.assertTrue(np.all(probs[:, 2] == 0),
                                f'expected zeroed final probability: {probs}')
                dose = x.alt_dosage
                self.assertTrue(np.all((dose >= 0) & (dose <= 2)),
                                f'dosages out of range: {dose}')

    def test_probabilities_out_of_range(self):
        ''' check we reject genotype probabilities outside 0-1

        Layout 1 always checked this, but layout 2 used to cast whatever it was
        given into an unsigned integer, so a negative or too-large probability
        wrapped around into an unrelated value.
        '''
        for layout, compression in [(1, 'zlib'), (2, 'zstd')]:
            for label, geno in [
                    ('negative', np.array([[-0.5, 0.7, 0.8]] * 3)),
                    ('above one', np.array([[2.0, 0.0, -1.0]] * 3)),
                    ('infinite', np.array([[np.inf, 0.0, 0.0]] * 3)),
                    ]:
                path = self.tmpdir / f'oob_{layout}_{label.replace(" ", "")}.bgen'
                with BgenWriter(path, 3, samples=['a', 'b', 'c'], layout=layout,
                                compression=compression) as bfile:
                    with self.assertRaisesRegex(ValueError, 'between 0 and 1'):
                        bfile.add_variant('var1', 'rs1', 'chr1', 10, ['A', 'C'],
                                          geno)

    def test_probabilities_sum_above_one(self):
        ''' check we reject genotype probabilities which sum above one

        Only the first probabilities are stored, and the reader infers the last
        as one minus their sum, so a sample summing above one cannot be encoded.
        '''
        for layout, compression in [(1, 'zlib'), (2, 'zstd')]:
            path = self.tmpdir / f'sum_{layout}.bgen'
            geno = np.array([[0.5, 0.9, 0.1]] * 3)
            with BgenWriter(path, 3, samples=['a', 'b', 'c'], layout=layout,
                            compression=compression) as bfile:
                with self.assertRaisesRegex(ValueError, 'sum to more than 1'):
                    bfile.add_variant('var1', 'rs1', 'chr1', 10, ['A', 'C'], geno)

    def test_phased_probabilities_out_of_range(self):
        ''' check phased probabilities are checked within each haplotype
        '''
        path = self.tmpdir / 'temp.bgen'
        with BgenWriter(path, 3, samples=['a', 'b', 'c']) as bfile:
            # the second haplotype of each sample sums above one
            geno = np.array([[0.3, 0.7, 0.8, 0.8]] * 3)
            with self.assertRaisesRegex(ValueError, 'sum to more than 1'):
                bfile.add_variant('var1', 'rs1', 'chr1', 10, ['A', 'C'], geno,
                                  phased=True)

            geno = np.array([[0.3, 0.7, -0.2, 1.2]] * 3)
            with self.assertRaisesRegex(ValueError, 'between 0 and 1'):
                bfile.add_variant('var2', 'rs2', 'chr1', 11, ['A', 'C'], geno,
                                  phased=True)

    def test_partly_missing_phased_sample(self):
        ''' check a sample with only some haplotypes missing is rejected

        A fully missing sample is stored with a missingness flag, but there is
        nowhere to record that just one haplotype of a sample is missing.
        '''
        path = self.tmpdir / 'temp.bgen'
        with BgenWriter(path, 2, samples=['a', 'b']) as bfile:
            geno = np.array([[0.5, 0.5, float('nan'), float('nan')],
                             [0.2, 0.8, 0.5, 0.5]])
            with self.assertRaisesRegex(ValueError, 'must encode all as missing'):
                bfile.add_variant('var1', 'rs1', 'chr1', 10, ['A', 'C'], geno,
                                  phased=True)

    def test_probabilities_within_rounding_tolerance(self):
        ''' check probabilities just outside 0-1 are still accepted

        The reader hands back float32, so probabilities which have been through
        a read/write round trip can sit a few float32 epsilons outside the legal
        range, and those must not be rejected.
        '''
        path = self.tmpdir / 'temp.bgen'
        geno = np.array([[-8.940697e-08, 0.5, 0.5],
                         [0.5, 0.5, 2.98023e-08],
                         [1 + 5.96046e-08, 0.0, 0.0],
                         ])
        with BgenWriter(path, 3, samples=['a', 'b', 'c']) as bfile:
            bfile.add_variant('var1', 'rs1', 'chr1', 10, ['A', 'C'], geno)

        with BgenReader(path, delay_parsing=True) as bfile:
            probs = bfile[0].probabilities
            self.assertTrue(np.all((probs >= -1e-7) & (probs <= 1 + 1e-7)))

    def test_probabilities_summing_below_one(self):
        ''' check probabilities summing to less than one are still written

        The bgen spec only requires the stored values to sum to at most the
        maximum, so a caller storing a partial distribution is legal.
        '''
        path = self.tmpdir / 'temp.bgen'
        geno = np.array([[0.1, 0.4, 0.05], [0.2, 0.1, 0.1], [0.0, 0.0, 0.0]])
        with BgenWriter(path, 3, samples=['a', 'b', 'c']) as bfile:
            bfile.add_variant('var1', 'rs1', 'chr1', 10, ['A', 'C'], geno,
                              bit_depth=16)

        with BgenReader(path, delay_parsing=True) as bfile:
            probs = bfile[0].probabilities
            # the first two probabilities are stored, so they round trip
            self.assertTrue(np.allclose(probs[:, :2], geno[:, :2], atol=1e-4))

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
        with BgenWriter(path, n_samples) as bfile:
            # construct a genotype array where the values
            a = np.linspace(0, 0.3, n_samples)
            b = np.linspace(0.7, 1, n_samples)
            geno = np.vstack([a, 1-a, b, 1-b]).T
            bfile.add_variant('var1', 'rs1', 'chr1', 10, ['A', 'C'], geno,
                                phased=1, bit_depth=8)

        time.sleep(1)

        with BgenReader(path, delay_parsing=True) as bfile:
            for x in bfile:
                probs = x.probabilities
                self.assertTrue(probs_close(geno[:, :-1], probs[:, :-1], bit_depth=8))
    
    def test_phased_data_different_sizes(self):
        '''checking writing phased data with a range of sample sizes'''
        for n_samples in range(1, 100):
            path = self.tmpdir / f'temp_{n_samples}.bgen'
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
            bfile.close()
    
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
    
    def test_add_variant_direct(self):
        ''' test adding variant data directly from BgenVar
        '''
        path_1 = self.tmpdir / 'temp.bgen'
        bfile = BgenWriter(path_1, 4, samples=['a', 'b', 'c', 'd'])
        ploidy = np.array([1, 2, 3, 3], dtype=np.uint8)
        geno = np.array([[0.1, 0.9, float('nan'), float('nan'), float('nan'), float('nan')],
                         [0.2, 0.8, 0.5, 0.5, float('nan'), float('nan')],
                         [float('nan'), float('nan'), float('nan'),
                          float('nan'), float('nan'), float('nan')],
                         [0.3, 0.7, 0.2, 0.8, 1, 0],
                         ])
        bfile.add_variant('var1', 'rs1', 'chr1', 10, ['A', 'C'], geno,
                          ploidy=ploidy, phased=1)

        # and add another, slightly different, variant
        geno = np.array([[0.1, 0.7, float('nan'), float('nan'), 0.3, 0.4],
                         [0.2, 0.7, 0.5, 0.5, float('nan'), float('nan')],
                         [float('nan'), float('nan'), float('nan'),
                          float('nan'), float('nan'), float('nan')],
                         [0.3, 0.1, 0.2, 0.8, 0.5, 0],
                         ])
        bfile.add_variant('var2', 'rs2', 'chr2', 20, ['G', 'TT'], geno,
                          ploidy=ploidy, phased=1)
        bfile.close()

        path_2 = self.tmpdir / 'temp2.bgen'
        with BgenReader(path_1) as bfile:
            samples = bfile.samples
            with BgenWriter(path_2, len(samples), samples) as output:
                for var in bfile:
                    output.add_variant_direct(var)

        with BgenReader(path_1) as bfile_1, BgenReader(path_2) as bfile_2:
            for var1, var2 in zip(bfile_1, bfile_2):
                self.assertEqual(var1.rsid, var2.rsid)
                self.assertEqual(var1.chrom, var2.chrom)
                self.assertEqual(var1.pos, var2.pos)
                self.assertEqual(var1.alleles, var2.alleles)
                self.assertEqual(var1.fileoffset, var2.fileoffset)
                self.assertEqual(var1.next_variant_offset, var2.next_variant_offset)

                # check all the nan values match
                self.assertTrue((np.isnan(var1.probabilities) == np.isnan(var2.probabilities)).all())

                # check all the non-nan values match
                mask = np.isfinite(var1.probabilities)
                self.assertTrue((var1.probabilities[mask] == var2.probabilities[mask]).all())

    def test_add_variant_direct_mismatched_file(self):
        ''' check we reject copying variant data into an incompatible bgen

        add_variant_direct copies the genotype block across as raw bytes, so the
        destination has to interpret those bytes the same way as the source. If
        the sample count, layout or compression differ, the copy used to be
        written anyway, and only failed later when the new bgen was read.
        '''
        source = self.tmpdir / 'source.bgen'
        geno = np.array([[0.1, 0.8, 0.1], [0.5, 0.25, 0.25], [0.1, 0.2, 0.7]])
        with BgenWriter(source, 3, samples=['a', 'b', 'c'], layout=2,
                        compression='zstd') as bfile:
            bfile.add_variant('var1', 'rs1', 'chr1', 10, ['A', 'C'], geno)

        mismatches = [
            ('samples', dict(n_samples=5, layout=2, compression='zstd'), 'samples'),
            ('layout', dict(n_samples=3, layout=1, compression='zlib'), 'layout'),
            ('compression', dict(n_samples=3, layout=2, compression='zlib'),
             'compression'),
            ('uncompressed', dict(n_samples=3, layout=2, compression=None),
             'compression'),
            ]
        with BgenReader(source) as bfile:
            var = bfile[0]
            for label, kwargs, message in mismatches:
                dest = self.tmpdir / f'dest_{label}.bgen'
                n = kwargs.pop('n_samples')
                with BgenWriter(dest, n, samples=[f'{i}' for i in range(n)],
                                **kwargs) as output:
                    with self.assertRaisesRegex(ValueError, message):
                        output.add_variant_direct(var)

    def test_add_variant_direct_all_layouts(self):
        ''' check copying variant data works for each layout and compression
        '''
        geno = np.array([[0.1, 0.8, 0.1], [0.5, 0.25, 0.25], [0.1, 0.2, 0.7]])
        for layout, compression in [(1, None), (1, 'zlib'), (2, None),
                                    (2, 'zlib'), (2, 'zstd')]:
            source = self.tmpdir / f'src_{layout}_{compression}.bgen'
            with BgenWriter(source, 3, samples=['a', 'b', 'c'], layout=layout,
                            compression=compression) as bfile:
                bfile.add_variant('var1', 'rs1', 'chr1', 10, ['A', 'C'], geno)

            dest = self.tmpdir / f'dst_{layout}_{compression}.bgen'
            with BgenReader(source) as bfile:
                with BgenWriter(dest, 3, samples=['a', 'b', 'c'], layout=layout,
                                compression=compression) as output:
                    for var in bfile:
                        output.add_variant_direct(var)

            # the copy has to hold the same genotypes as the source
            with BgenReader(source) as first, BgenReader(dest) as second:
                for var1, var2 in zip(first, second):
                    self.assertTrue((var1.probabilities == var2.probabilities).all(),
                                    f'layout={layout} compression={compression}')

    def test_variant_file_properties(self):
        ''' check a BgenVar reports the layout and compression it was read with
        '''
        geno = np.array([[0.1, 0.8, 0.1], [0.5, 0.25, 0.25], [0.1, 0.2, 0.7]])
        for layout, compression in [(1, None), (1, 'zlib'), (2, None),
                                    (2, 'zlib'), (2, 'zstd')]:
            path = self.tmpdir / f'props_{layout}_{compression}.bgen'
            with BgenWriter(path, 3, samples=['a', 'b', 'c'], layout=layout,
                            compression=compression) as bfile:
                bfile.add_variant('var1', 'rs1', 'chr1', 10, ['A', 'C'], geno)

            with BgenReader(path) as bfile:
                var = bfile[0]
                self.assertEqual(var.layout, layout)
                self.assertEqual(var.compression, compression)
                self.assertEqual(var.n_samples, 3)
                # and these have to agree with the file header
                self.assertEqual(var.layout, bfile.header.layout)
                self.assertEqual(var.compression, bfile.header.compression)
                self.assertEqual(var.n_samples, bfile.header.nsamples)

    def test_multiple_read_writes(self):
        ''' check values pass correctly through multiple rounds of read/writes
        '''
        first_path = self.tmpdir / 'temp1.bgen'
        second_path = self.tmpdir / 'temp2.bgen'
        
        bit_depth = 9
        max_val = (2 ** bit_depth) - 1
        half = max_val // 2
        integer_values = [[0,                max_val,     0],
                          [half,             half + 1,    0],
                          [half + 1,         half,        0],
                          [(max_val - 1),    1,           0],
                         ]
        integer_values = np.array(integer_values)
        geno = integer_values / max_val
        
        with BgenWriter(first_path, n_samples=len(geno)) as bfile:
            bfile.add_variant('var1', 'rs1', 'chr1', 10, ['A', 'C'], geno, bit_depth=bit_depth)
        
        with BgenWriter(second_path, n_samples=len(geno)) as out_bfile:
            with BgenReader(first_path, delay_parsing=True) as bfile:
                probs = next(bfile).probabilities
                as_integers = (probs * max_val).round()
                self.assertTrue(probs_close(geno, probs, bit_depth=bit_depth))
                self.assertTrue((integer_values == as_integers).all())
                out_bfile.add_variant(
                    'var1', 'rs1', 'chr1', 10, ['A', 'C'], probs, bit_depth=bit_depth)
        
        with BgenReader(second_path, delay_parsing=True) as bfile:
                probs = next(bfile).probabilities
                as_integers = (probs * max_val).round()
                self.assertTrue(probs_close(geno, probs, bit_depth=bit_depth))
                self.assertTrue((integer_values == as_integers).all())
    
    def test_all_possible_genotypes(self):
        ''' check all possible values in the range available to the bit depths
        '''
        max_samples = 10000000
        for bit_depth in range(1, 24):
            first_path = self.tmpdir / f'temp_{bit_depth}.v1.bgen'
            second_path = self.tmpdir / f'temp2_{bit_depth}.v2.bgen'
            
            max_val = (2 ** bit_depth) - 1
            increment = 1
            if (max_val / 2) > max_samples:
                increment = ((max_val + 1) / 2) / max_samples
            
            integer_values = np.arange(0, max_val + 1, increment)
            
            remainder = max_val - (integer_values[::2] + integer_values[::-2])
            
            integer_values = np.array([integer_values[::2], integer_values[::-2],
                                       remainder], dtype=np.uint32)
            integer_values = np.ascontiguousarray(integer_values.T)
            geno = integer_values / max_val
            
            # write first round
            with BgenWriter(first_path, n_samples=len(geno)) as bfile:
                bfile.add_variant('var1', 'rs1', 'chr1', 10, ['A', 'C'], geno,
                                  bit_depth=bit_depth)
            
            # check the first write
            with BgenWriter(second_path, n_samples=len(geno)) as out_bfile:
                with BgenReader(first_path, delay_parsing=True) as bfile:
                    probs = next(bfile).probabilities
                    as_integers = (probs * max_val).round()
                    self.assertTrue(probs_close(geno, probs, bit_depth=bit_depth))
                    self.assertTrue((integer_values == as_integers).all())
                    out_bfile.add_variant(
                        'var1', 'rs1', 'chr1', 10, ['A', 'C'], probs, bit_depth=bit_depth)
            
            # check the re-written bgen
            with BgenReader(second_path, delay_parsing=True) as bfile:
                    probs = next(bfile).probabilities
                    as_integers = (probs * max_val).round()
                    self.assertTrue(probs_close(geno, probs, bit_depth=bit_depth))
                    self.assertTrue((integer_values == as_integers).all())

    def test_long_sample_id(self):
        ''' check we reject a sample ID too long for the bgen length field

        Sample IDs are stored behind a two byte length, so a longer ID used to
        wrap that field to a small value while still writing every character.
        The sample block length stayed self-consistent, so the file looked fine,
        but each later ID was read from the wrong offset and the whole bgen
        became unreadable.
        '''
        geno = np.array([[0.1, 0.8, 0.1], [0.5, 0.25, 0.25]])

        # the longest ID which still fits must be written and read back intact
        path = self.tmpdir / 'longest.bgen'
        longest = 'x' * 65535
        with BgenWriter(path, 2, samples=[longest, 'b']) as bfile:
            bfile.add_variant('var1', 'rs1', 'chr1', 10, ['A', 'C'], geno)
        with BgenReader(path) as bfile:
            self.assertEqual(bfile.samples, [longest, 'b'])

        for length in [65536, 100000]:
            path = self.tmpdir / f'long_{length}.bgen'
            with self.assertRaisesRegex(ValueError, 'sample ID is too long'):
                BgenWriter(path, 2, samples=['x' * length, 'b'])

    def test_long_variant_identifiers(self):
        ''' check we reject variant IDs, rsIDs and chromosomes which are too long

        These are all stored behind a two byte length, so anything longer used to
        wrap the length field and leave a file which could not be read back.
        '''
        geno = np.array([[0.1, 0.8, 0.1], [0.5, 0.25, 0.25]])
        fields = {'varid': 'variant ID', 'rsid': 'rsID', 'chrom': 'chromosome'}

        for field, message in fields.items():
            # the longest value which still fits must round trip
            path = self.tmpdir / f'longest_{field}.bgen'
            longest = 'x' * 65535
            args = {'varid': 'var1', 'rsid': 'rs1', 'chrom': 'chr1'}
            args[field] = longest
            with BgenWriter(path, 2, samples=['a', 'b']) as bfile:
                bfile.add_variant(pos=10, alleles=['A', 'C'], genotypes=geno,
                                  **args)
            with BgenReader(path) as bfile:
                self.assertEqual(getattr(next(iter(bfile)), field), longest)

            for length in [65536, 100000]:
                path = self.tmpdir / f'long_{field}_{length}.bgen'
                args[field] = 'x' * length
                with BgenWriter(path, 2, samples=['a', 'b']) as bfile:
                    with self.assertRaisesRegex(ValueError, f'{message} is too long'):
                        bfile.add_variant(pos=10, alleles=['A', 'C'],
                                          genotypes=geno, **args)

    def test_long_identifier_leaves_file_usable(self):
        ''' check rejecting an over-long identifier doesn't corrupt the file

        The identifiers are checked before anything is written, so the rejected
        variant must not appear in the file, must not be counted in the header,
        and must not disturb the variants written either side of it.
        '''
        geno = np.array([[0.1, 0.8, 0.1], [0.5, 0.25, 0.25]])
        path = self.tmpdir / 'temp.bgen'
        with BgenWriter(path, 2, samples=['a', 'b']) as bfile:
            bfile.add_variant('var1', 'rs1', 'chr1', 10, ['A', 'C'], geno)
            with self.assertRaisesRegex(ValueError, 'variant ID is too long'):
                bfile.add_variant('x' * 70000, 'rs2', 'chr1', 20, ['A', 'C'],
                                  geno)
            bfile.add_variant('var3', 'rs3', 'chr1', 30, ['A', 'C'], geno)

        with BgenReader(path) as bfile:
            self.assertEqual(len(bfile), 2)
            self.assertEqual([x.varid for x in bfile], ['var1', 'var3'])
            for var in bfile:
                self.assertTrue(probs_close(geno, var.probabilities,
                                            bit_depth=8))

    def test_duplicate_and_empty_sample_ids(self):
        ''' check duplicate and empty sample IDs are still accepted

        The bgen spec doesn't require sample IDs to be unique or non-empty, so
        these are written as given rather than rejected.
        '''
        geno = np.array([[0.1, 0.8, 0.1], [0.5, 0.25, 0.25]])
        for samples in [['a', 'a'], ['', ''], ['', 'b']]:
            path = self.tmpdir / f'dup_{"_".join(samples)}.bgen'
            with BgenWriter(path, 2, samples=samples) as bfile:
                bfile.add_variant('var1', 'rs1', 'chr1', 10, ['A', 'C'], geno)
            with BgenReader(path) as bfile:
                self.assertEqual(bfile.samples, samples)
