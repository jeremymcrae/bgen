''' check a bgen decodes to the same floats whichever code path handles a sample

8-bit data has vectorised decoders, plus a scalar loop for the samples left over, and
which runs depends on the CPU. These tests compare the paths through values any CPU
gives, and check exact dosages, which also depend on the -ffp-contract=off build flag.
'''

from pathlib import Path
import tempfile
import unittest

import numpy as np

from bgen import BgenReader, BgenWriter

class TestDecodeConsistency(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.tmpdir = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def write(self, name, geno, ploidy=2, phased=False):
        path = self.tmpdir / name
        with BgenWriter(path, len(geno), compression=None) as bfile:
            bfile.add_variant('v', 'rs', '1', 1, ['A', 'C'], geno, ploidy=ploidy,
                              phased=phased, bit_depth=8)
        return path

    def read(self, path, attr):
        with BgenReader(path) as bfile:
            return getattr(next(iter(bfile)), attr)

    def lookup_table(self):
        ''' the float each stored 8-bit value decodes to, via the unphased decoder
        '''
        k = np.arange(300) % 256
        geno = np.stack([k / 255, np.zeros(len(k)), 1 - k / 255], axis=1)
        path = self.write('unphased.bgen', geno)
        probs = self.read(path, 'probabilities')
        table = np.empty(256, dtype=np.float32)
        table[k] = probs[:, 0]
        # every occurrence of a value must have decoded the same way
        self.assertTrue(np.array_equal(table[k], probs[:, 0]))
        return table

    def test_phased_8bit_matches_the_lookup_table(self):
        ''' phased 8-bit probabilities decode exactly as the unphased ones do
        '''
        table = self.lookup_table()
        # odd sizes leave haplotypes for the scalar tail after any vectorised blocks
        for n in [1, 7, 8, 9, 33, 300, 1001]:
            k = (np.arange(n * 2) * 7) % 256
            hap = k / 255
            geno = np.empty((n, 4))
            geno[:, 0] = hap[0::2]
            geno[:, 1] = 1 - hap[0::2]
            geno[:, 2] = hap[1::2]
            geno[:, 3] = 1 - hap[1::2]
            path = self.write(f'phased_{n}.bgen', geno, phased=True)
            probs = self.read(path, 'probabilities')
            first = probs[:, [0, 2]].ravel()
            second = probs[:, [1, 3]].ravel()
            self.assertTrue(np.array_equal(first, table[k]), f'n={n}')
            self.assertTrue(np.array_equal(second, table[255 - k]), f'n={n}')

    def test_phased_dosage_is_not_fused(self):
        ''' a phased dosage is each haplotype scaled, then summed, without fusing
        '''
        n = 1001
        rng = np.random.default_rng(3)
        k = rng.integers(0, 256, (n, 2))
        geno = np.empty((n, 4))
        geno[:, 0] = k[:, 0] / 255
        geno[:, 1] = 1 - k[:, 0] / 255
        geno[:, 2] = k[:, 1] / 255
        geno[:, 3] = 1 - k[:, 1] / 255
        path = self.write('dosage.bgen', geno, phased=True)
        alt = self.read(path, 'alt_dosage')

        # float32 throughout, and numpy rounds after every operation, so this is exact
        factor = np.float32(1) / np.float32(255)
        scaled = k.astype(np.float32) * factor
        ref = scaled[:, 0] + scaled[:, 1]
        self.assertTrue(np.array_equal(alt, np.float32(2) - ref))

        # a mix of ploidies uses the general loop, which has to give the same values
        ploidy = np.where(np.arange(n) % 3 == 0, 1, 2).astype(np.uint8)
        geno[ploidy == 1, 2:] = np.nan
        path = self.write('dosage_mixed.bgen', geno, ploidy=ploidy, phased=True)
        mixed = self.read(path, 'alt_dosage')
        ref[ploidy == 1] = scaled[ploidy == 1, 0]
        self.assertTrue(np.array_equal(mixed, ploidy.astype(np.float32) - ref))

    def test_inferred_probability_is_never_negative(self):
        ''' the probability inferred as the remainder is never negative
        '''
        n = 2000
        rng = np.random.default_rng(11)
        for bit_depth in [9, 12, 16, 20, 31]:
            first = rng.random(n)
            geno = np.stack([first, 1 - first, np.zeros(n)], axis=1)
            path = self.tmpdir / f'remainder_{bit_depth}.bgen'
            with BgenWriter(path, n, compression=None) as bfile:
                bfile.add_variant('v', 'rs', '1', 1, ['A', 'C'], geno,
                                  bit_depth=bit_depth)
            probs = self.read(path, 'probabilities')
            self.assertTrue((probs >= 0).all(), f'bit_depth={bit_depth}: {probs.min()}')
            # rounding can still leave it a float step above zero, but no further
            self.assertTrue((probs[:, 2] < 1e-6).all(), f'bit_depth={bit_depth}')

        # and the reference files, which are not written by this package
        data = Path(__file__).parent / 'data'
        for name in ['example.12bits.bgen', 'example.16bits.bgen']:
            with BgenReader(data / name) as bfile:
                for var in bfile:
                    probs = var.probabilities
                    self.assertFalse((probs < 0).any(), f'{name} {var.rsid}')

    def test_phased_varying_ploidy_layout(self):
        ''' each sample's haplotypes fill its row in order, padded with nan after them
        '''
        rng = np.random.default_rng(5)
        n = 500
        ploidy = rng.integers(0, 4, n).astype(np.uint8)
        max_ploidy = int(ploidy.max())
        geno = np.full((n, max_ploidy * 2), np.nan)
        for i, p in enumerate(ploidy):
            k = rng.integers(0, 256, p)
            geno[i, 0:2 * p:2] = k / 255
            geno[i, 1:2 * p:2] = 1 - k / 255
        # mark some samples as missing
        missing = np.nonzero(ploidy > 0)[0][:5]
        geno[missing] = np.nan

        path = self.write('ragged.bgen', geno, ploidy=ploidy, phased=True)
        probs = self.read(path, 'probabilities')
        self.assertEqual(probs.shape, geno.shape)
        self.assertEqual(probs.dtype, np.float32)
        self.assertTrue(np.array_equal(np.isnan(probs), np.isnan(geno)))
        self.assertTrue(np.allclose(probs, geno, atol=0.5 / 255, equal_nan=True))
