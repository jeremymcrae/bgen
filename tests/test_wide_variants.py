''' test variants whose per sample probability counts are large

Two counters used to be too narrow for variants with many alleles, or with a high ploidy.
Both were found while profiling the writer, and both reject or mangle data that is
perfectly legal:

  * missing_genotypes counted nans in a uint16. A sample's probability count is
    C(ploidy + n_alleles - 1, n_alleles - 1), which passes 65535 at 6 alleles and a ploidy
    of 21, so a fully missing sample counted 244 instead of 65780 and the writer rejected
    it as only partially missing.
  * the reader's probability_bytes() filled a 64-entry table by calling get_max_probs for
    every ploidy from 0 to 63, regardless of the ploidy the variant actually declares. For
    12 or more alleles the ploidy 63 entry exceeds a 32-bit count, so a variant with a
    varying ploidy of 1-2 could be written but not read back.

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
        self.tmpdir = Path(tempfile.mkdtemp())

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

        The nan count wrapped a uint16 here, so the writer rejected the variant with
        'samples with any missing genotype must encode all as missing' even though every
        value for that sample was nan.
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

        6 alleles at a ploidy of 20 gives 53130 probabilities, which fits a uint16, so
        this passed before the fix and must still pass.
        '''
        path = self.tmpdir / 'below.bgen'
        n_probs = math.comb(20 + 6 - 1, 6 - 1)
        self.assertLess(n_probs, 65535)
        _, got = self._write_and_read(path, 6, 20, n_samples=2, missing=(0, ))
        self.assertEqual(got.shape, (2, n_probs))
        self.assertTrue(np.all(np.isnan(got[0])))

    def test_many_alleles_with_varying_ploidy(self):
        ''' a multiallelic variant with a low varying ploidy reads back

        The reader sized a per ploidy table by computing every ploidy up to 63, and for
        12 or more alleles that count does not fit 32 bits - so these variants wrote
        successfully but raised 'more than 2^32 probabilities' when read.
        '''
        ploidy = [1, 2, 1, 2]
        for n_alleles in [12, 16, 20, 24, 30]:
            with self.subTest(n_alleles=n_alleles):
                # the ploidy 63 entry is what used to be computed unnecessarily
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
        ''' the constant ploidy path was never affected, so check it is unchanged '''
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
