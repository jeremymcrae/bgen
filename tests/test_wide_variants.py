''' test variants whose per sample probability counts are large

Two counts have to be wide enough for variants with many alleles, or with a high ploidy:

  * the writer counts a sample's missing probabilities. A sample's probability count is
    C(ploidy + n_alleles - 1, n_alleles - 1), which passes 65535 at 6 alleles and a ploidy
    of 21.
  * the reader's probability_bytes() sizes a variant from the ploidies it declares. For
    12 or more alleles the count at a ploidy of 63 exceeds 32 bits, which must not stop a
    variant with a ploidy of 1-2 being read.

The arrays here stay small: it is the count per sample that has to be large, not the
cohort, so a handful of samples is enough.
'''

import math
from pathlib import Path
import tempfile
import unittest

import numpy as np

from bgen import BgenReader, BgenWriter


class TestWideVariants(unittest.TestCase):
    ''' check variants with large per sample probability counts round trip '''

    def setUp(self):
        # TemporaryDirectory rather than mkdtemp, which never cleans up after itself
        self.tmp = tempfile.TemporaryDirectory()
        self.tmpdir = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def _write_and_read(self, path, n_alleles, ploidy, n_samples=4, missing=()):
        ''' write a variant with the given alleles and ploidy, then read it back

        ploidy may be a scalar or a per sample array. Rows for missing samples are
        all nan, which is how the format stores a missing sample.
        '''
        ploidy_arr = np.asarray(ploidy, dtype=np.uint8) if not np.isscalar(ploidy) \
            else None
        widest_ploidy = int(max(ploidy_arr)) if ploidy_arr is not None else int(ploidy)
        widest = math.comb(widest_ploidy + n_alleles - 1, n_alleles - 1)

        probs = np.full((n_samples, widest), np.nan)
        for i in range(n_samples):
            this_ploidy = int(ploidy_arr[i]) if ploidy_arr is not None else int(ploidy)
            k = math.comb(this_ploidy + n_alleles - 1, n_alleles - 1)
            probs[i, :k] = 1.0 / k
        for i in missing:
            probs[i, :] = np.nan

        alleles = [f'A{i}' for i in range(n_alleles)]
        with BgenWriter(path, n_samples,
                        samples=[f's{i}' for i in range(n_samples)],
                        compression=None) as bfile:
            bfile.add_variant(varid='v', rsid='rs', chrom='01', pos=1,
                              alleles=alleles, genotypes=probs,
                              ploidy=ploidy if ploidy_arr is None else ploidy_arr)

        with BgenReader(path) as bfile:
            var = next(iter(bfile))
            return probs, var.probabilities

    def test_fully_missing_sample_beyond_16_bits(self):
        ''' a missing sample with more than 65535 probabilities is stored as missing
        '''
        # each case exceeds 65535 probabilities per sample, and the one below it does not
        for n_alleles, ploidy in [(6, 21), (5, 33), (7, 16)]:
            with self.subTest(n_alleles=n_alleles, ploidy=ploidy):
                n_probs = math.comb(ploidy + n_alleles - 1, n_alleles - 1)
                self.assertGreater(n_probs, 65535)

                path = self.tmpdir / f'miss{n_alleles}_{ploidy}.bgen'
                _, got = self._write_and_read(path, n_alleles, ploidy, n_samples=2,
                                              missing=(0, ))
                self.assertEqual(got.shape, (2, n_probs))
                # the missing sample stays missing, and the other one does not
                self.assertTrue(np.all(np.isnan(got[0])))
                self.assertFalse(np.any(np.isnan(got[1])))

    def test_missing_sample_below_the_boundary(self):
        ''' the case just under 65535 probabilities has to keep working

        6 alleles at a ploidy of 20 gives 53130 probabilities, which fits a uint16.
        '''
        path = self.tmpdir / 'below.bgen'
        n_probs = math.comb(20 + 6 - 1, 6 - 1)
        self.assertLess(n_probs, 65535)
        _, got = self._write_and_read(path, 6, 20, n_samples=2, missing=(0, ))
        self.assertEqual(got.shape, (2, n_probs))
        self.assertTrue(np.all(np.isnan(got[0])))

    def test_many_alleles_with_varying_ploidy(self):
        ''' a multiallelic variant with a low varying ploidy reads back

        For 12 or more alleles, the probability count at a ploidy of 63 does not fit
        32 bits, which must not matter for a variant that never reaches that ploidy.
        '''
        ploidy = [1, 2, 1, 2]
        for n_alleles in [12, 16, 20, 24, 30]:
            with self.subTest(n_alleles=n_alleles):
                # the count for a ploidy of 63, which these variants never reach
                unreachable = math.comb(63 + n_alleles - 1, n_alleles - 1)
                self.assertGreater(unreachable, 0xFFFFFFFF)

                path = self.tmpdir / f'ma{n_alleles}.bgen'
                _, got = self._write_and_read(path, n_alleles, ploidy)
                widest = math.comb(2 + n_alleles - 1, n_alleles - 1)
                self.assertEqual(got.shape, (4, widest))
                # the ploidy 1 samples fill fewer columns than the ploidy 2 samples
                for i, ploid in enumerate(ploidy):
                    k = math.comb(ploid + n_alleles - 1, n_alleles - 1)
                    self.assertFalse(np.any(np.isnan(got[i, :k])),
                                     f'sample {i} should have {k} probabilities')

    def test_many_alleles_constant_ploidy_still_reads(self):
        ''' many alleles at a constant ploidy can be read '''
        for n_alleles in [12, 24]:
            with self.subTest(n_alleles=n_alleles):
                path = self.tmpdir / f'const{n_alleles}.bgen'
                _, got = self._write_and_read(path, n_alleles, 2)
                widest = math.comb(2 + n_alleles - 1, n_alleles - 1)
                self.assertEqual(got.shape, (4, widest))
                self.assertFalse(np.any(np.isnan(got)))

    def test_dosage_of_wide_variant_is_refused(self):
        ''' dosage needs a biallelic variant, so these must still be refused clearly '''
        path = self.tmpdir / 'dose.bgen'
        self._write_and_read(path, 12, [1, 2, 1, 2])
        with BgenReader(path) as bfile:
            var = next(iter(bfile))
            with self.assertRaises(ValueError):
                var.minor_allele_dosage


if __name__ == '__main__':
    unittest.main()
