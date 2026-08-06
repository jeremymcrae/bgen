

from pathlib import Path
import unittest
import tempfile
import logging
import os
import sqlite3
import threading
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
    
    def _collect_warnings(self, geno, bit_depth, **kwargs):
        ''' write a variant and return the bit_depth warnings it logged
        '''
        messages = []

        class Collect(logging.Handler):
            def emit(self, record):
                messages.append(record.getMessage())

        path = self.tmpdir / f'warn_{bit_depth}_{len(messages)}_{id(geno)}.bgen'
        handler = Collect()
        logger = logging.getLogger()
        logger.addHandler(handler)
        try:
            with BgenWriter(path, geno.shape[0],
                            samples=[f's{i}' for i in range(geno.shape[0])]) as bfile:
                bfile.add_variant('var1', 'rs1', 'chr1', 10, ['A', 'C'], geno,
                                  bit_depth=bit_depth, **kwargs)
        finally:
            logger.removeHandler(handler)
        return [x for x in messages if 'bit_depth' in x]

    def test_coarse_bit_depth_warns(self):
        ''' a bit depth too coarse for the probabilities has to say so

        Probabilities are stored as integers out of 2**bit_depth - 1, so a low depth does
        not round so much as replace the values: [0.3, 0.3, 0.4] comes back from a one bit
        file as [0, 1, 0]. That used to happen silently.

        Which depths warn depends on the values, since the test is whether the loss beats
        what the default depth would cost. These probabilities are coarse enough to warn
        up to a depth of 5, and by 6 they already round to within the tolerance.
        '''
        geno = np.array([[0.3, 0.3, 0.4]])
        for bit_depth in range(1, 6):
            with self.subTest(bit_depth=bit_depth):
                messages = self._collect_warnings(geno, bit_depth)
                self.assertEqual(len(messages), 1, messages)
                self.assertIn(f'bit_depth={bit_depth}', messages[0])

        # and the values really are changed, which is what the warning is about
        path = self.tmpdir / 'coarse.bgen'
        with BgenWriter(path, 1, samples=['a']) as bfile:
            bfile.add_variant('var1', 'rs1', 'chr1', 10, ['A', 'C'], geno, bit_depth=1)
        with BgenReader(path) as bfile:
            probs = np.asarray(next(iter(bfile)).probabilities)
        self.assertTrue(np.allclose(probs, [[0, 1, 0]]), probs)

    def test_adequate_bit_depth_stays_quiet(self):
        ''' the default depth and above must not warn, or the warning gets ignored

        The default of 8 gives a step of 1/255, which keeps two decimal places, so
        ordinary imputed probabilities are stored well enough. Warning merely because
        they are not exactly representable would fire on nearly every file written.
        '''
        rng = np.random.default_rng(9)
        geno = rng.random((20, 3))
        geno /= geno.sum(axis=1)[:, None]
        for bit_depth in [8, 12, 16, 24, 32]:
            with self.subTest(bit_depth=bit_depth):
                self.assertEqual(self._collect_warnings(geno, bit_depth), [])

    def test_exactly_representable_values_stay_quiet(self):
        ''' data a low depth can hold exactly must not warn

        Hard called genotypes are exact at any depth, so storing them in a single bit is
        a legitimate way to save space. Keying the warning on the depth alone would scold
        people doing that.
        '''
        hard = np.zeros((6, 3))
        hard[np.arange(6), np.arange(6) % 3] = 1.0
        for bit_depth in [1, 2, 4, 8]:
            with self.subTest(bit_depth=bit_depth, data='hard calls'):
                self.assertEqual(self._collect_warnings(hard, bit_depth), [])

        # values that are exact multiples for their depth are equally fine
        for bit_depth in [2, 3, 4]:
            factor = 2 ** bit_depth - 1
            geno = np.array([[1 / factor, 2 / factor, 1 - 3 / factor]])
            with self.subTest(bit_depth=bit_depth, data='exact multiples'):
                self.assertEqual(self._collect_warnings(geno, bit_depth), [])

    def test_missing_genotypes_do_not_warn(self):
        ''' samples stored as missing lose nothing, whatever the depth
        '''
        geno = np.full((3, 3), np.nan)
        for bit_depth in [1, 4, 8]:
            with self.subTest(bit_depth=bit_depth):
                self.assertEqual(self._collect_warnings(geno, bit_depth), [])

    def _assert_missing_layout(self, base, probs, missing, bit_depth, width=None,
                               note=''):
        ''' every missing sample reads back missing, and every other one reads back

        The tolerance covers both sources of error: the bit depth the values were
        quantised to, and the float32 the reader returns them in. Above about 22 bits the
        float32 dominates, so a bit depth term on its own rejects correct output.
        '''
        tol = 1.0 / (2 ** bit_depth - 1) + 4 * float(np.spacing(np.float32(1.0)))
        for i in range(len(base)):
            k = width[i] if width is not None else base.shape[1]
            got = probs[i][:k]
            want = base[i][:k]
            where = f'sample {i}, bit_depth={bit_depth}, missing={missing}{note}'
            if i in missing:
                self.assertTrue(np.all(np.isnan(got)),
                                f'{where} should read back as missing')
            else:
                self.assertFalse(np.any(np.isnan(got)),
                                 f'{where} should not read back as missing')
                np.testing.assert_allclose(got, want, atol=tol, rtol=0,
                                           err_msg=where)

    def test_missing_samples_roundtrip_at_every_bit_depth(self):
        ''' missing samples must read back as missing, wherever they sit

        Their probabilities are not written value by value, since every one of them
        encodes as zero and the buffer already holds zeroes. That means the encoder
        advances the bit offset by hand, so a sample written after a missing one lands
        at an offset that was never checked against a stored value. Read the file back
        and confirm both the missing rows and their neighbours survive.
        '''
        n = 9
        rng = np.random.default_rng(42)
        base = rng.random((n, 3))
        base /= base.sum(axis=1, keepdims=True)
        samples = [f's{i}' for i in range(n)]

        placements = [[0], [n - 1], [0, n - 1], [4], [3, 4, 5],
                      list(range(0, n, 2)), list(range(n))]
        for bit_depth in [1, 2, 5, 8, 11, 16, 17, 24, 32]:
            for missing in placements:
                geno = base.copy()
                geno[missing] = np.nan
                path = self.tmpdir / 'miss.bgen'
                with BgenWriter(path, n, samples=samples) as bfile:
                    bfile.add_variant('v', 'rs', '01', 10, ['A', 'C'], geno,
                                      bit_depth=bit_depth)

                with BgenReader(path) as bfile:
                    probs = next(iter(bfile)).probabilities
                self._assert_missing_layout(base, probs, missing, bit_depth)

    def test_missing_phased_samples_roundtrip(self):
        ''' the same for phased data, which skips per haplotype

        The phased encoder walks haplotypes within a sample, so skipping a missing one
        has to advance by ploidy * (n_probs - 1) rather than by a single row.
        '''
        n = 8
        rng = np.random.default_rng(7)
        samples = [f's{i}' for i in range(n)]
        for ploidy in [1, 2, 3]:
            base = rng.random((n, 2 * ploidy))
            for hap in range(ploidy):
                cols = slice(hap * 2, (hap + 1) * 2)
                base[:, cols] /= base[:, cols].sum(axis=1, keepdims=True)
            for bit_depth in [1, 8, 16, 32]:
                for missing in [[0], [n - 1], [2, 3], list(range(0, n, 2))]:
                    geno = base.copy()
                    geno[missing] = np.nan
                    path = self.tmpdir / 'missp.bgen'
                    with BgenWriter(path, n, samples=samples) as bfile:
                        bfile.add_variant(
                            'v', 'rs', '01', 10, ['A', 'C'], geno,
                            ploidy=np.full(n, ploidy, dtype=np.uint8),
                            phased=True, bit_depth=bit_depth)

                    with BgenReader(path) as bfile:
                        var = next(iter(bfile))
                        probs = var.probabilities
                        self.assertTrue(var.is_phased)
                    self._assert_missing_layout(base, probs, missing, bit_depth,
                                                note=f', ploidy={ploidy}')

    def test_missing_samples_with_varying_ploidy_roundtrip(self):
        ''' and with a ploidy that differs per sample

        Each sample then stores a different number of probabilities, so the skip has to
        use that sample's own width rather than the widest.
        '''
        n = 9
        rng = np.random.default_rng(19)
        samples = [f's{i}' for i in range(n)]
        ploidy = np.tile([1, 2, 3], n // 3).astype(np.uint8)
        widths = [math.comb(int(p) + 1, 1) for p in ploidy]
        base = np.full((n, max(widths)), np.nan)
        for i, ploid in enumerate(ploidy):
            row = rng.random(widths[i])
            base[i, :widths[i]] = row / row.sum()

        for bit_depth in [1, 8, 15, 32]:
            for missing in [[0], [n - 1], [1, 2], list(range(0, n, 3))]:
                geno = base.copy()
                geno[missing] = np.nan
                path = self.tmpdir / 'missv.bgen'
                with BgenWriter(path, n, samples=samples) as bfile:
                    bfile.add_variant('v', 'rs', '01', 10, ['A', 'C'], geno,
                                      ploidy=ploidy, bit_depth=bit_depth)

                with BgenReader(path) as bfile:
                    var = next(iter(bfile))
                    probs = var.probabilities
                    self.assertEqual(list(var.ploidy), list(ploidy))
                self._assert_missing_layout(base, probs, missing, bit_depth,
                                            width=widths)

    def test_biallelic_8bit_matches_neighbouring_shapes(self):
        ''' the common shape is encoded by a specialised loop, which must agree

        Biallelic, unphased, constant ploidy 2, bit depth 8 is what almost all real data
        looks like, so it skips the generic encoder. Nothing about the stored bytes may
        change as a result. There is no way to route the same variant through both
        encoders from Python, so instead check the parts that are observable: the values
        read back, and that the shapes either side of the dispatch still behave.
        '''
        n = 50
        rng = np.random.default_rng(101)
        geno = rng.random((n, 3))
        geno /= geno.sum(axis=1, keepdims=True)
        samples = [f's{i}' for i in range(n)]

        # the dispatched shape, and its immediate neighbours which take the generic path
        shapes = [
            ('dispatched', dict(bit_depth=8), geno),
            ('bit_depth 7', dict(bit_depth=7), geno),
            ('bit_depth 9', dict(bit_depth=9), geno),
            ('explicit ploidy', dict(bit_depth=8,
                                     ploidy=np.full(n, 2, dtype=np.uint8)), geno),
        ]
        for label, kwargs, probs in shapes:
            path = self.tmpdir / 'shape.bgen'
            with BgenWriter(path, n, samples=samples) as bfile:
                bfile.add_variant('v', 'rs', '01', 10, ['A', 'C'], probs,
                                  **kwargs)
            with BgenReader(path) as bfile:
                var = next(iter(bfile))
                got = var.probabilities
                self.assertEqual(list(var.ploidy), [2] * n, label)
                self.assertFalse(var.is_phased, label)
            depth = kwargs['bit_depth']
            tol = 1 / (2 ** depth - 1) + 4 * float(np.spacing(np.float32(1.0)))
            np.testing.assert_allclose(got, probs, atol=tol, rtol=0,
                                       err_msg=label)

    def test_biallelic_8bit_stores_every_quantised_value(self):
        ''' the specialised encoder must cover the whole 0 to 255 range

        It scales a running total rather than each probability, and clamps so the total
        never decreases, so walk the first probability across every step the bit depth
        can represent and confirm each one round trips.
        '''
        n = 256
        first = np.arange(n) / 255.0
        geno = np.zeros((n, 3))
        geno[:, 0] = first
        geno[:, 1] = 1.0 - first
        samples = [f's{i}' for i in range(n)]

        path = self.tmpdir / 'steps.bgen'
        with BgenWriter(path, n, samples=samples) as bfile:
            bfile.add_variant('v', 'rs', '01', 10, ['A', 'C'], geno,
                              bit_depth=8)

        with BgenReader(path) as bfile:
            got = next(iter(bfile)).probabilities
        # each stored value is exact at this bit depth, so only float32 error remains
        np.testing.assert_allclose(got, geno, atol=1e-6, rtol=0)

    def test_biallelic_8bit_rejects_bad_rows(self):
        ''' the specialised encoder keeps the same rejections as the generic one

        It writes the probability checks out by hand rather than looping, so each one has
        to still be there: a partially missing row, a row summing above one, and a
        negative probability.
        '''
        n = 8
        samples = [f's{i}' for i in range(n)]
        base = np.full((n, 3), 1 / 3)

        partial = base.copy()
        partial[3, 1] = np.nan
        cases = [
            (partial, 'must encode all as missing'),
            (np.full((n, 3), 0.5), 'sum to more than 1'),
            (np.array([[-0.5, 0.75, 0.75]] * n), 'must be between 0 and 1'),
        ]
        for geno, message in cases:
            path = self.tmpdir / 'bad.bgen'
            with BgenWriter(path, n, samples=samples) as bfile:
                with self.assertRaisesRegex(ValueError, message):
                    bfile.add_variant('v', 'rs', '01', 10, ['A', 'C'], geno,
                                      bit_depth=8)

    def test_biallelic_8bit_covers_every_batch_boundary(self):
        ''' the specialised encoder vectorises, so sample counts must not matter

        It encodes several samples at a time where the hardware allows (four under avx2,
        two under NEON) and finishes the remainder one at a time. A sample count that is
        not a whole number of batches therefore takes both paths, and an off by one in
        either would corrupt or drop the samples at the join. Walk every count across two
        batch widths so no alignment is missed.
        '''
        rng = np.random.default_rng(202)
        for n in range(1, 18):
            with self.subTest(n_samples=n):
                geno = rng.random((n, 3))
                geno /= geno.sum(axis=1, keepdims=True)
                samples = [f's{i}' for i in range(n)]
                path = self.tmpdir / f'batch{n}.bgen'
                with BgenWriter(path, n, samples=samples) as bfile:
                    bfile.add_variant('v', 'rs', '01', 10, ['A', 'C'], geno,
                                      bit_depth=8)
                with BgenReader(path) as bfile:
                    var = next(iter(bfile))
                    got = var.probabilities
                self.assertEqual(got.shape, (n, 3))
                np.testing.assert_allclose(got, geno, atol=1 / 255 + 1e-6, rtol=0)

    def test_biallelic_8bit_handles_missing_in_every_lane(self):
        ''' a batch mixing missing and present samples must encode both correctly

        The vectorised encoder tests a whole batch at once, and missing samples fail that
        test, so it has to separate a batch that is partly missing from one that is
        unencodable. Placing missing samples at every position within and across batches
        covers each combination of lanes, which is where a mask error would show up as a
        neighbour being written as missing, or a missing sample being written as data.
        '''
        rng = np.random.default_rng(303)
        n = 16
        base = rng.random((n, 3))
        base /= base.sum(axis=1, keepdims=True)
        samples = [f's{i}' for i in range(n)]
        # every subset of the first four samples, so each lane pattern of a batch occurs,
        # plus patterns that straddle a batch edge
        patterns = [[i for i in range(4) if mask & (1 << i)] for mask in range(16)]
        patterns += [[3, 4], [2, 3, 4, 5], [n - 1], [0, n - 1],
                     list(range(0, n, 2)), list(range(n))]
        for missing in patterns:
            with self.subTest(missing=tuple(missing)):
                geno = base.copy()
                geno[missing] = np.nan
                path = self.tmpdir / 'lanes.bgen'
                with BgenWriter(path, n, samples=samples) as bfile:
                    bfile.add_variant('v', 'rs', '01', 10, ['A', 'C'], geno,
                                      bit_depth=8)
                with BgenReader(path) as bfile:
                    var = next(iter(bfile))
                    got = var.probabilities
                    ploidy = np.array(var.ploidy)
                absent = np.zeros(n, dtype=bool)
                absent[missing] = True
                # a missing sample reads back as nan, and a present one must not
                self.assertTrue(np.isnan(got[absent]).all())
                self.assertFalse(np.isnan(got[~absent]).any())
                np.testing.assert_allclose(got[~absent], base[~absent],
                                           atol=1 / 255 + 1e-6, rtol=0)
                # the ploidy byte still has to read as 2 for every sample
                self.assertEqual(list(ploidy), [2] * n)

    def test_biallelic_8bit_rounds_ties_away_from_zero(self):
        ''' the vectorised encoder must round halves the same way as the scalar one

        The obvious vector rounding instruction rounds halves to even, whereas the scalar
        encoder rounds them away from zero, and probabilities that are multiples of 1/510
        land exactly on a half after scaling. Those are the values where the two rules
        disagree, so store them and check what comes back matches the away from zero rule.
        '''
        # (k + 0.5) / 255 scales to exactly k + 0.5, the tie for every stored value
        ties = np.array([(k + 0.5) / 255.0 for k in range(255)])
        n = len(ties)
        geno = np.zeros((n, 3))
        geno[:, 0] = ties
        geno[:, 1] = 1.0 - ties
        samples = [f's{i}' for i in range(n)]

        path = self.tmpdir / 'ties.bgen'
        with BgenWriter(path, n, samples=samples) as bfile:
            bfile.add_variant('v', 'rs', '01', 10, ['A', 'C'], geno, bit_depth=8)

        with BgenReader(path) as bfile:
            got = next(iter(bfile)).probabilities
        # away from zero sends k + 0.5 up to k + 1, so the stored value is one higher than
        # truncation would give. Half to even would send every other tie down instead
        expected = (np.arange(255) + 1) / 255.0
        np.testing.assert_allclose(got[:, 0], expected, atol=1e-6, rtol=0)

    def test_biallelic_8bit_rejects_running_total_over_tolerance(self):
        ''' the running total after two probabilities is checked, not just the final one

        check_probability allows a probability down to -PROB_TOLERANCE, so a negative third
        value can pull a row back under the limit after the first two have already gone
        over it. Such a row must still be rejected, because the encoder stores the running
        total and a value above the maximum cannot be represented.

        The vectorised encoder tests the intermediate total separately from the final one,
        and dropping that test is invisible to any row whose probabilities are all
        non-negative: the final total is then the larger of the two, so checking it alone
        appears to be enough. These rows are the only ones that tell the two apart, and
        the window is about 1e-6 wide, so they have to be constructed rather than sampled.
        '''
        tolerance = 1e-6
        n = 8
        samples = [f's{i}' for i in range(n)]
        # each row has p0 + p1 above 1 + tolerance, but p0 + p1 + p2 back inside it
        rows = [
            (0.5, 0.5 + 1.5 * tolerance, -tolerance),
            (0.5, 0.5 + 2 * tolerance, -tolerance),
            (0.25, 0.75 + 1.5 * tolerance, -tolerance / 2),
        ]
        for p0, p1, p2 in rows:
            self.assertGreater(p0 + p1, 1 + tolerance)
            self.assertLessEqual(p0 + p1 + p2, 1 + tolerance)
            # place the row in each position of a batch, so no lane is left untested
            for lane in range(4):
                with self.subTest(row=(p0, p1, p2), lane=lane):
                    geno = np.full((n, 3), 1 / 3)
                    geno[lane] = [p0, p1, p2]
                    path = self.tmpdir / f'running{lane}.bgen'
                    with BgenWriter(path, n, samples=samples) as bfile:
                        with self.assertRaisesRegex(ValueError, 'sum to more than 1'):
                            bfile.add_variant('v', 'rs', '01', 10, ['A', 'C'], geno,
                                              bit_depth=8)

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

    def test_reading_8bit_matches_across_vector_boundary(self):
        ''' the vectorised 8-bit read must agree with its own scalar tail

        The unphased diploid fast path decodes eight samples at a time and leaves the
        remainder to a scalar loop, so a sample's probabilities are produced by different
        code depending on its position in the file. Give every sample identical genotypes,
        so anything other than one distinct value per column means the two halves disagree,
        and vary the sample count around multiples of eight so the split falls in a
        different place each time.

        Probabilities are chosen as multiples of 1/255 which are stored exactly, so the
        comparison is of bit patterns rather than of rounding. Values that sum above one at
        8 bits are included because the third probability is inferred as 255 - first -
        second, which goes negative there and has to clamp identically in both halves.
        '''
        cases = [(1 / 255, 1 / 255), (0.2, 0.3), (1.0, 0.0), (0.0, 1.0),
                 (254 / 255, 1 / 255), (0.5, 0.5)]
        for first, second in cases:
            for n in [1, 2, 7, 8, 9, 15, 16, 17, 23, 24, 25, 100]:
                with self.subTest(probs=(first, second), n=n):
                    geno = np.tile(np.array([[first, second, 1 - first - second]]), (n, 1))
                    path = self.tmpdir / f'boundary_{n}.bgen'
                    with BgenWriter(path, n, samples=[f's{i}' for i in range(n)]) as bfile:
                        bfile.add_variant('v', 'rs', '01', 10, ['A', 'C'], geno,
                                          ploidy=2, bit_depth=8)
                    with BgenReader(path) as bfile:
                        probs = next(iter(bfile)).probabilities

                    self.assertEqual(probs.shape, (n, 3))
                    for col in range(3):
                        values = np.unique(probs[:, col])
                        self.assertEqual(len(values), 1,
                                         f'column {col} of {n} identical samples read back '
                                         f'as {len(values)} distinct values, {values}, so '
                                         f'the vectorised and scalar reads disagree')
                    # and the values themselves must still be the ones stored
                    np.testing.assert_allclose(probs[0], geno[0], atol=1 / 255 + 1e-6,
                                               rtol=0)

    def test_reading_8bit_dosage_matches_across_vector_boundary(self):
        ''' the vectorised 8-bit dosage must agree with its own scalar tail

        test_reading_8bit_matches_across_vector_boundary covers probabilities. Dosage takes
        a separate route, ref_dosage_fast, which has its own vector loop and its own scalar
        tail, so it needs its own check. Give every sample identical genotypes: more than
        one distinct dosage can only mean the two halves disagree.

        The dosage count is first * 2 + second, and the stored first probability is swept
        over its whole range, because the two halves used to differ only for counts whose
        exact value is not representable as a float.
        '''
        for first in [1, 23, 51, 77, 128, 200, 254, 255]:
            for n in [1, 2, 7, 8, 9, 15, 16, 17, 31, 33, 100]:
                with self.subTest(first=first, n=n):
                    probs = np.tile(np.array([[first / 255, 0.0, 1 - first / 255]]), (n, 1))
                    path = self.tmpdir / f'dose_boundary_{n}.bgen'
                    with BgenWriter(path, n, samples=[f's{i}' for i in range(n)],
                                    compression=None) as bfile:
                        bfile.add_variant('v', 'rs', '01', 10, ['A', 'C'], probs,
                                          ploidy=2, bit_depth=8)
                    with BgenReader(path) as bfile:
                        dose = np.asarray(next(iter(bfile)).alt_dosage)

                    values = np.unique(dose)
                    self.assertEqual(len(values), 1,
                                     f'{n} identical samples read back as {len(values)} '
                                     f'distinct dosages, {values}, so the vectorised and '
                                     f'scalar dosage reads disagree')
                    # and the dosage must still be the one stored, 2 - 2 * first / 255
                    np.testing.assert_allclose(dose, 2.0 - 2.0 * first / 255,
                                               atol=1e-6, rtol=0)

    def test_reading_8bit_above_maximum_across_vector_boundary(self):
        ''' the vectorised 8-bit read must clamp an over-large sum like the scalar one

        The third probability is inferred as 255 - first - second, which goes negative for a
        bgen storing values that sum above the maximum. Older versions of this package wrote
        such files. The scalar read clamps the index to zero; the vectorised read has its own
        clamp, and without it the gather would index before the lookup table and return
        whatever lies in front of it in memory.

        test_reading_probabilities_above_maximum covers the same ground for three samples,
        which only ever reach the scalar loop. Vary the count so the vectorised half handles
        the affected samples too.
        '''
        for n in [3, 8, 9, 16, 17, 100]:
            with self.subTest(n_samples=n):
                path = self.tmpdir / f'above_max_{n}.bgen'
                geno = np.tile(np.array([[1.0, 0.0, 0.0]]), (n, 1))
                with BgenWriter(path, n, samples=[f's{i}' for i in range(n)],
                                compression=None) as bfile:
                    bfile.add_variant('var1', 'rs1', 'chr1', 10, ['A', 'C'], geno,
                                      ploidy=2, bit_depth=8)

                # each sample stores (255, 0), so rewrite to (255, 255), summing to 510 and
                # making the inferred final probability 255 - 510 = -255
                raw = bytearray(path.read_bytes())
                patched = 0
                for i in range(len(raw) - 1):
                    if raw[i] == 255 and raw[i + 1] == 0:
                        raw[i + 1] = 255
                        patched += 1
                self.assertEqual(patched, n)
                path.write_bytes(bytes(raw))

                with BgenReader(path, delay_parsing=True) as bfile:
                    probs = next(iter(bfile)).probabilities

                self.assertTrue(np.all((probs >= 0) & (probs <= 1)),
                                f'probabilities out of range: {probs}')
                # clamped to zero, rather than gathered from before the lookup table
                self.assertTrue(np.all(probs[:, 2] == 0),
                                f'expected zeroed final probability: {probs}')
                # every sample is identical, so the vectorised and scalar halves must agree
                for col in range(3):
                    self.assertEqual(len(np.unique(probs[:, col])), 1,
                                     f'column {col} disagrees across the vector boundary: '
                                     f'{np.unique(probs[:, col])}')

    def test_malformed_values_stay_in_range(self):
        ''' a malformed bgen must not yield probabilities or dosages outside their range

        A probability lies in 0-1 and a dosage in 0-ploidy. When the stored values sum above
        what the bit depth allows, the inferred remainder goes negative, so each decoder
        clamps it. There are five: the 8-bit probability and dosage paths each have a
        vectorised half and a scalar tail, and the general path covers other bit depths.
        Every one has to clamp, or the values it returns are impossible.

        Sample counts straddle the vector strides, so some runs are handled wholly by the
        scalar code and others hand over to it partway through. Samples are identical, so
        any variation in the result means the two halves of a decoder disagree.
        '''
        for bit_depth in (8, 16, 32):
            for n in (4, 12, 40, 64):
                with self.subTest(bit_depth=bit_depth, n=n):
                    path = self.tmpdir / f'range_{bit_depth}_{n}.bgen'
                    patched = self._write_above_maximum(path, n, bit_depth, 'all')
                    self.assertEqual(patched, n)

                    with BgenReader(path, delay_parsing=True) as bfile:
                        var = next(iter(bfile))
                        probs = var.probabilities
                        alt = np.asarray(var.alt_dosage, dtype=float)
                    with BgenReader(path, delay_parsing=True) as bfile:
                        minor = np.asarray(next(iter(bfile)).minor_allele_dosage,
                                           dtype=float)

                    self.assertTrue(np.all((probs >= 0) & (probs <= 1)),
                                    f'probabilities outside 0-1: {np.unique(probs)}')
                    for label, dose in (('alt', alt), ('minor', minor)):
                        self.assertTrue(np.all((dose >= 0) & (dose <= 2)),
                                        f'{label} dosage outside 0-2: {np.unique(dose)}')
                    # identical samples, so each decoder's halves must agree
                    for col in range(probs.shape[1]):
                        self.assertEqual(len(np.unique(probs[:, col])), 1,
                                         f'probability column {col} differs across the '
                                         f'vector boundary: {np.unique(probs[:, col])}')
                    self.assertEqual(len(np.unique(alt)), 1,
                                     f'alt dosage differs across the vector boundary: '
                                     f'{np.unique(alt)}')
                    # the clamp puts the ref dosage at the ploidy, so the alt dosage the
                    # swap derives from it is zero. Every decoder has to agree on which end
                    # to clamp to, or the same malformed file reads differently by bit depth
                    self.assertTrue(np.all(alt == 0),
                                    f'expected zeroed alt dosage, got {np.unique(alt)}')
                    self.assertTrue(np.all(probs[:, -1] == 0),
                                    f'expected zeroed final probability, got '
                                    f'{np.unique(probs[:, -1])}')

    def _read_warnings(self, path, prop):
        ''' read every variant via prop and return the malformed-bgen warnings logged
        '''
        messages = []

        class Collect(logging.Handler):
            def emit(self, record):
                messages.append(record.getMessage())

        handler = Collect()
        logger = logging.getLogger()
        logger.addHandler(handler)
        try:
            with BgenReader(path, delay_parsing=True) as bfile:
                for var in bfile:
                    getattr(var, prop)
        finally:
            logger.removeHandler(handler)
        return [x for x in messages if 'malformed' in x]

    def test_reading_probabilities_above_maximum_warns(self):
        ''' a bgen storing probabilities above the maximum has to say so

        Such a file is malformed, and the probability inferred for the affected samples is
        not the one the file describes, so reading it silently hands back numbers that
        cannot be trusted. Older versions of this package wrote these files.

        The detection lives in five separate decoders, and each has to notice on its own.
        Only samples 8 and beyond reach the vectorised 8-bit halves, so a variant with
        fewer than that tests only the scalar code; putting the bad samples at the front or
        the back of a large variant separates the vector loop from its tail. Bit depths
        above 8 take the general decoder instead, and the dosage properties take two more
        decoders again. Each combination below is therefore a different code path, not a
        repeat.
        '''
        # (bit depth, sample count, which samples to corrupt, what it exercises)
        #
        # bit depths are kept to whole bytes so the patching below can find the stored
        # values; 12 bits would pack them across byte boundaries. 16 bits reaches the same
        # general decoder that 12 would.
        cases = [
            (8, 4, 'all', '8-bit scalar only, too few samples to vectorise'),
            (8, 64, 'first', '8-bit vector loop'),
            (8, 12, 'last', '8-bit scalar tail after a vector pass'),
            (16, 4, 'all', 'general decoder, no vectorisation'),
            (16, 64, 'first', 'general decoder over many samples'),
            (16, 40, 'last', 'general decoder, corrupt samples at the end'),
        ]
        for bit_depth, n, which, description in cases:
            for prop in ('probabilities', 'alt_dosage', 'minor_allele_dosage'):
                with self.subTest(bit_depth=bit_depth, n=n, corrupt=which, prop=prop):
                    path = self.tmpdir / f'malformed_{bit_depth}_{n}_{which}.bgen'
                    self._write_above_maximum(path, n, bit_depth, which)
                    messages = self._read_warnings(path, prop)
                    self.assertEqual(len(messages), 1,
                                     f'{description}: expected one warning for {prop}, '
                                     f'got {messages}')
                    self.assertIn('rs0/var0', messages[0])

    def _write_above_maximum(self, path, n, bit_depth, which):
        ''' write a bgen whose stored probabilities sum above what the bit depth allows

        The writer will not produce one, since it rejects such probabilities, so write a
        valid variant and then raise the stored values by patching the bytes. Every sample
        stores its first two probabilities as (max, 0), so setting the second to max as
        well makes the pair sum to twice the maximum.

        Bit depths are whole bytes so the stored values can be located by their byte
        pattern. which selects the samples to corrupt: 'first' and 'last' put them at one
        end of the variant, so the vectorised and scalar halves of the 8-bit decoders can be
        told apart, since only one of them then sees a bad sample.
        '''
        geno = np.tile(np.array([[1.0, 0.0, 0.0]]), (n, 1))
        with BgenWriter(path, n, samples=[f's{i}' for i in range(n)],
                        compression=None) as bfile:
            bfile.add_variant('var0', 'rs0', '01', 10, ['A', 'C'], geno, ploidy=2,
                              bit_depth=bit_depth)

        width = (bit_depth + 7) // 8      # bytes per stored probability
        maxval = (1 << bit_depth) - 1
        top = maxval.to_bytes(width, 'little')
        zero = (0).to_bytes(width, 'little')
        raw = bytearray(path.read_bytes())

        # find each sample's pair of stored probabilities, which the writer laid down in
        # sample order, and raise the second of the pair to the maximum
        starts = []
        pos = raw.find(top + zero)
        while pos != -1:
            starts.append(pos)
            pos = raw.find(top + zero, pos + 2 * width)
        self.assertEqual(len(starts), n,
                         f'found {len(starts)} samples to patch, expected {n}')

        if which == 'all':
            targets = starts
        elif which == 'first':
            targets = starts[:4]
        else:
            targets = starts[-4:]
        for start in targets:
            raw[start + width:start + 2 * width] = top
        path.write_bytes(bytes(raw))
        return len(targets)

    def test_valid_bgen_does_not_warn_as_malformed(self):
        ''' a well formed bgen must never be reported as malformed

        The check for an over-large sum is a comparison against the remainder the decoders
        already compute. In the general path that remainder is a float subtraction, so it
        can land marginally below zero on a valid variant through rounding alone, and a
        naive test would warn about ordinary files at some bit depths. A false warning is
        worse than none, since it would teach users to ignore the real one.
        '''
        rng = np.random.default_rng(17)
        for bit_depth in range(1, 33):
            with self.subTest(bit_depth=bit_depth):
                geno = rng.dirichlet([1, 1, 1], size=25)
                path = self.tmpdir / f'valid_{bit_depth}.bgen'
                with BgenWriter(path, 25, samples=[f's{i}' for i in range(25)]) as bfile:
                    bfile.add_variant('v', 'rs', '01', 10, ['A', 'C'], geno,
                                      bit_depth=bit_depth)
                for prop in ('probabilities', 'alt_dosage'):
                    self.assertEqual(self._read_warnings(path, prop), [],
                                     f'valid bgen warned at bit_depth={bit_depth} '
                                     f'via {prop}')

    def test_malformed_warning_does_not_carry_between_variants(self):
        ''' one malformed variant must not make the valid ones look malformed too

        Each variant owns its own genotype state, so the flag recording an over-large sum is
        scoped to a single variant. This locks that in: were the flag ever moved to the file
        or reader, one bad variant would condemn every variant read after it, and the warning
        would stop being useful for locating the problem.
        '''
        n = 24
        good = np.tile(np.array([[0.2, 0.3, 0.5]]), (n, 1))
        # the malformed variant sits in the middle, so there are valid variants read both
        # before and after it
        path = self.tmpdir / 'carryover.bgen'
        with BgenWriter(path, n, samples=[f's{i}' for i in range(n)],
                        compression=None) as bfile:
            for i in range(5):
                geno = np.tile(np.array([[1.0, 0.0, 0.0]]), (n, 1)) if i == 2 else good
                bfile.add_variant(f'var{i}', f'rs{i}', '01', i + 1, ['A', 'C'], geno,
                                  ploidy=2, bit_depth=8)

        # only var2 stores (255, 0) per sample, the valid variants store other values
        raw = bytearray(path.read_bytes())
        patched = 0
        for i in range(len(raw) - 1):
            if raw[i] == 255 and raw[i + 1] == 0:
                raw[i + 1] = 255
                patched += 1
        self.assertEqual(patched, n)
        path.write_bytes(bytes(raw))

        for prop in ('probabilities', 'alt_dosage', 'minor_allele_dosage'):
            with self.subTest(prop=prop):
                messages = self._read_warnings(path, prop)
                self.assertEqual(len(messages), 1,
                                 f'expected only var2 to warn, got {messages}')
                self.assertIn('rs2/var2', messages[0])

        # and reading one valid variant repeatedly must stay quiet, even after the
        # malformed one has been read through the same reader
        with BgenReader(path, delay_parsing=True) as bfile:
            variants = list(bfile)
            variants[2].probabilities          # trips the flag
            messages = []

            class Collect(logging.Handler):
                def emit(self, record):
                    messages.append(record.getMessage())

            handler = Collect()
            logging.getLogger().addHandler(handler)
            try:
                for _ in range(3):
                    variants[0].probabilities
                    variants[4].probabilities
            finally:
                logging.getLogger().removeHandler(handler)
            self.assertEqual([x for x in messages if 'malformed' in x], [],
                             'a valid variant warned after a malformed one was read')

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

    def test_empty_alleles_are_written_as_given(self):
        ''' an empty allele is stored and read back unchanged

        The bgen spec doesn't require alleles to be non-empty, so these are written as
        given, as sample IDs are. The minor allele is worked out when it is asked for,
        so an empty string coming back means the empty allele is the minor one, rather
        than the allele being unavailable.
        '''
        geno = np.array([[0.1, 0.8, 0.1], [0.5, 0.25, 0.25]])
        for i, alleles in enumerate([['A', ''], ['', 'A'], ['', '']]):
            with self.subTest(alleles=alleles):
                path = self.tmpdir / f'empty_{i}.bgen'
                with BgenWriter(path, 2, samples=['a', 'b']) as bfile:
                    bfile.add_variant('var1', 'rs1', 'chr1', 10, alleles, geno)
                with BgenReader(path) as bfile:
                    var = next(iter(bfile))
                    self.assertEqual(list(var.alleles), alleles)
                    # whichever allele is the minor one, it is one of the two given
                    self.assertIn(var.minor_allele, alleles)

    def test_empty_allele_reported_only_when_it_is_the_minor_one(self):
        ''' an empty minor allele tracks the frequencies, so it is a real answer

        The empty string has to be reported when the empty allele is the less common
        one, and the other allele when it is not. That is what separates a truthful
        answer from a value that just means nothing was worked out.
        '''
        n_samples = 100
        # 4 hom ref and 32 het gives the first allele a frequency of 0.2
        geno = np.zeros((n_samples, 3))
        geno[:4] = [1, 0, 0]
        geno[4:36] = [0, 1, 0]
        geno[36:] = [0, 0, 1]

        for alleles, expected in [(['', 'A'], ''), (['A', ''], 'A')]:
            with self.subTest(alleles=alleles):
                path = self.tmpdir / f'emptyminor_{alleles[0]}.bgen'
                with BgenWriter(path, n_samples,
                                samples=[f's{i}' for i in range(n_samples)]) as bfile:
                    bfile.add_variant('var1', 'rs1', 'chr1', 10, alleles, geno)
                with BgenReader(path) as bfile:
                    var = next(iter(bfile))
                    # the first allele is the rarer one, so it is the minor allele
                    self.assertEqual(var.minor_allele, expected)
                    self.assertAlmostEqual(np.nanmean(var.minor_allele_dosage),
                                           0.4, delta=0.02)

    def test_duplicate_alleles_warn_but_are_written(self):
        ''' duplicate alleles are allowed, with a warning

        The spec permits them and the genotype data is still coherent, so the variant
        is written. The minor allele just cannot say which of the two it means.
        '''
        geno = np.array([[0.1, 0.8, 0.1], [0.5, 0.25, 0.25]])
        path = self.tmpdir / 'dupalleles.bgen'
        with BgenWriter(path, 2, samples=['a', 'b']) as bfile:
            with self.assertLogs(level='WARNING') as logs:
                bfile.add_variant('var1', 'rs1', 'chr1', 10, ['A', 'A'], geno)
        self.assertTrue(any('duplicate alleles' in x for x in logs.output))

        with BgenReader(path) as bfile:
            var = next(iter(bfile))
            self.assertEqual(list(var.alleles), ['A', 'A'])
            # the variant still reads back, and names one of its alleles
            self.assertEqual(var.minor_allele, 'A')

    def test_distinct_alleles_do_not_warn(self):
        ''' the warning must not fire for ordinary variants

        Alleles differing only by case are distinct, so they are not duplicates.
        '''
        geno = np.array([[0.1, 0.8, 0.1], [0.5, 0.25, 0.25]])
        for alleles in [['A', 'C'], ['a', 'A'], ['GT', 'GA']]:
            with self.subTest(alleles=alleles):
                path = self.tmpdir / 'distinct.bgen'
                # assertNoLogs needs python 3.10, so collect the records directly
                messages = []

                class Collect(logging.Handler):
                    def emit(self, record):
                        messages.append(record.getMessage())

                handler = Collect()
                logger = logging.getLogger()
                logger.addHandler(handler)
                try:
                    with BgenWriter(path, 2, samples=['a', 'b']) as bfile:
                        bfile.add_variant('var1', 'rs1', 'chr1', 10, alleles, geno)
                finally:
                    logger.removeHandler(handler)
                self.assertEqual([x for x in messages if 'duplicate' in x], [])
                with BgenReader(path) as bfile:
                    self.assertEqual(list(next(iter(bfile)).alleles), alleles)

    def test_long_alleles_are_still_accepted(self):
        ''' indels can be long, so length is not what makes an allele invalid
        '''
        geno = np.array([[0.1, 0.8, 0.1], [0.5, 0.25, 0.25]])
        alleles = ['A', 'ACGT' * 50]
        path = self.tmpdir / 'longallele.bgen'
        with BgenWriter(path, 2, samples=['a', 'b']) as bfile:
            bfile.add_variant('var1', 'rs1', 'chr1', 10, alleles, geno)
        with BgenReader(path) as bfile:
            self.assertEqual(list(next(iter(bfile)).alleles), alleles)

    def _break_bgen_fd(self, path):
        ''' close the descriptor the writer holds for its bgen

        Writes to the bgen then fail, while the .bgi index keeps working, so the failure
        is confined to the data file the way an I/O error on it would be. Limiting the
        file size instead would break the sqlite index too.
        '''
        target = os.path.realpath(str(path))
        for name in os.listdir('/proc/self/fd'):
            try:
                if os.readlink(f'/proc/self/fd/{name}') == target:
                    os.close(int(name))
                    return True
            except (OSError, ValueError):
                pass
        return False

    @unittest.skipUnless(sys.platform.startswith('linux'),
                         'needs /proc to find the bgen descriptor')
    def test_write_error_is_not_followed_by_a_bogus_offset(self):
        ''' a writer that has already failed must not report made up file offsets

        The offsets come from tellp(), which returns -1 once the stream has failed. That
        used to be cast to an unsigned 64 bit offset, so a caller that caught a write
        error and carried on recorded 18446744073709551615 as a variant's position in the
        index, with no error at all. Anything reading that index would then seek to
        nonsense.
        '''
        geno = np.array([[0.1, 0.8, 0.1]] * 3)
        src = self.tmpdir / 'offsets_src.bgen'
        with BgenWriter(src, 3, samples=['a', 'b', 'c']) as bfile:
            bfile.add_variant('var1', 'rs1', 'chr1', 10, ['A', 'C'], geno)
        reader = BgenReader(src)
        self.addCleanup(reader.close)
        variant = next(iter(reader))

        path = self.tmpdir / 'offsets.bgen'
        bfile = BgenWriter(path, 3, samples=['a', 'b', 'c'])
        bfile.add_variant('var1', 'rs1', 'chr1', 10, ['A', 'C'], geno)
        self.assertTrue(self._break_bgen_fd(path), 'could not find the bgen descriptor')

        # each call has to keep refusing, rather than inventing an offset. A caller that
        # catches the first error and carries on is exactly how a bogus offset used to
        # reach the index
        for i in range(3):
            with self.assertRaises(OSError) as ctx:
                bfile.add_variant(f'v{i}', f'rs{i}', 'chr1', i + 20, ['A', 'C'], geno)
            self.assertIn('position', str(ctx.exception))

        # the direct copy is the site that used to return (uint64)-1 without any error
        with self.assertRaises(OSError) as ctx:
            bfile.add_variant_direct(variant)
        self.assertIn('position', str(ctx.exception))

        # nothing bogus should have reached the index
        bgi = Path(str(path) + '.bgi')
        starts = []
        if bgi.exists():
            con = sqlite3.connect(bgi)
            starts = [row[0] for row in
                      con.execute('select file_start_position from Variant').fetchall()]
            con.close()
        self.assertTrue(all(0 < x < 2 ** 32 for x in starts), starts)

    @unittest.skipUnless(hasattr(os, 'mkfifo'), 'needs fifos, which windows lacks')
    def test_writing_to_a_pipe_reports_a_useful_error(self):
        ''' a bgen needs a seekable output, and saying so beats an iostream error

        The header records where the variants start and how many there are, and both are
        only known at the end, so the writer seeks back to fill them in. On a pipe the
        position cannot even be read, which used to surface as
        "basic_ios::clear: iostream error" from somewhere in the middle of a write.
        '''
        path = self.tmpdir / 'pipe.fifo'
        os.mkfifo(path)

        def drain():
            with open(path, 'rb') as handle:
                handle.read()

        # a reader has to be attached, or opening the fifo for writing blocks
        reader = threading.Thread(target=drain, daemon=True)
        reader.start()
        try:
            with self.assertRaises(OSError) as ctx:
                bfile = BgenWriter(path, 2, samples=['a', 'b'])
                geno = np.array([[0.1, 0.8, 0.1], [0.5, 0.25, 0.25]])
                bfile.add_variant('var1', 'rs1', 'chr1', 10, ['A', 'C'], geno)
                bfile.close()
            self.assertIn('position', str(ctx.exception))
        finally:
            reader.join(timeout=30)
