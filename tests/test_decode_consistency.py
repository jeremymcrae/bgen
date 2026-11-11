''' check a bgen decodes to the same floats whichever code path handles a sample

8-bit data has vectorised decoders, plus a scalar loop for the samples left over, and
which runs depends on the CPU. These tests compare the paths through values any CPU
gives.
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
