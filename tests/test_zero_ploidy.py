''' tests for dosage of samples which have no alleles

A ploidy zero sample has just one possible genotype, the empty one, and since the
final probability of each group is never stored it occupies no bytes at all.

The dosage reader used to read one probability per sample regardless, so for these
it consumed the next sample's data, shifted every later sample's bit offset along,
and ran off the end of the block. The probability reader handled them correctly, so
the two disagreed.
'''

from pathlib import Path
import tempfile
import unittest

import numpy as np

from bgen import BgenReader, BgenWriter

nan = float('nan')


def genotypes(ploidy):
    ''' genotype probabilities matching a per sample ploidy array

    Every sample of the same ploidy gets the same probabilities, so any sample read
    at a shifted offset stands out against its peers. A sample stores one
    probability fewer than it has genotypes, so later columns are left as nan.

    A ploidy zero sample holds the empty genotype with certainty. That has to be
    spelled out, since an all nan row means the sample is missing instead.
    '''
    geno = np.full((len(ploidy), 3), nan)
    geno[ploidy == 2, :3] = [0.2, 0.5, 0.3]
    geno[ploidy == 1, :2] = [0.4, 0.6]
    geno[ploidy == 0, 0] = 1.0
    return geno


class TestZeroPloidy(unittest.TestCase):
    ''' check dosages for samples without any alleles '''

    def setUp(self):
        self.tmpdir = tempfile.TemporaryDirectory()
        self.path = Path(self.tmpdir.name) / 'zero_ploidy.bgen'

    def tearDown(self):
        self.tmpdir.cleanup()

    def write(self, ploidy, geno=None, **kwargs):
        ''' write a single variant with the given per sample ploidy '''
        n = len(ploidy)
        if geno is None:
            geno = genotypes(ploidy)
        with BgenWriter(self.path, n, samples=[f's{i}' for i in range(n)]) as bfile:
            bfile.add_variant('var', 'rs1', '01', 1, ['A', 'C'], geno,
                              ploidy=ploidy, **kwargs)

    def read(self, attr):
        ''' read one attribute of the first variant '''
        with BgenReader(self.path) as bfile:
            return np.asarray(getattr(next(iter(bfile)), attr), dtype=np.float64)

    def test_zero_ploidy_dosage_is_zero(self):
        ''' a sample with no alleles carries no copies of either allele '''
        ploidy = np.array([0, 2, 0, 2], dtype=np.uint8)
        self.write(ploidy)
        for attr in ['alt_dosage', 'minor_allele_dosage']:
            dose = self.read(attr)
            self.assertTrue((dose[ploidy == 0] == 0).all(), attr)

    def test_other_samples_keep_their_dosage(self):
        ''' zero ploidy samples must not shift the samples after them

        Every diploid here has identical probabilities, so all must get the same
        dosage no matter where the zero ploidy samples sit.
        '''
        for pattern in [[0, 2], [2, 0], [0, 0, 2], [2, 2, 0], [0, 1, 2]]:
            for reps in [1, 3, 50]:
                ploidy = np.array(pattern * reps, dtype=np.uint8)
                with self.subTest(pattern=pattern, reps=reps):
                    self.write(ploidy)
                    dose = self.read('alt_dosage')
                    diploid = dose[ploidy == 2]
                    if diploid.size:
                        self.assertTrue(np.allclose(diploid, 1.1, atol=0.01),
                                        f'diploid dosages differ: {diploid}')
                    haploid = dose[ploidy == 1]
                    if haploid.size:
                        self.assertTrue(np.allclose(haploid, 0.6, atol=0.01),
                                        f'haploid dosages differ: {haploid}')

    def test_dosage_agrees_with_probabilities(self):
        ''' the two readers must not disagree about the same variant '''
        ploidy = np.array([(i % 3) for i in range(30)], dtype=np.uint8)
        self.write(ploidy)
        probs = self.read('probabilities')
        dose = self.read('alt_dosage')
        expected = np.where(ploidy == 2, probs[:, 1] + 2 * np.nan_to_num(probs[:, 2]),
                            np.where(ploidy == 1, probs[:, 1], 0.0))
        np.testing.assert_allclose(dose, expected, atol=0.01)

    def test_zero_ploidy_at_the_end_of_the_cohort(self):
        ''' the final sample having no alleles used to read past the block

        Nothing follows it to misalign, so the dosages come out right either way,
        but the read still ran off the end of the buffer.
        '''
        ploidy = np.array([2] * 19 + [0], dtype=np.uint8)
        self.write(ploidy)
        dose = self.read('alt_dosage')
        self.assertTrue(np.allclose(dose[:19], 1.1, atol=0.01))
        self.assertEqual(dose[19], 0)

    def test_every_sample_has_zero_ploidy(self):
        ''' a variant nobody has any alleles for is all zero dosage '''
        ploidy = np.zeros(10, dtype=np.uint8)
        geno = np.full((10, 1), 1.0)
        self.write(ploidy, geno=geno)
        dose = self.read('alt_dosage')
        self.assertTrue((dose == 0).all(), f'{dose}')

    def test_zero_ploidy_across_bit_depths(self):
        ''' the offsets are counted in bits, so check several widths '''
        ploidy = np.array([0, 2, 1, 0, 2], dtype=np.uint8)
        for bit_depth in [1, 2, 4, 8, 12, 16, 24]:
            with self.subTest(bit_depth=bit_depth):
                self.write(ploidy, bit_depth=bit_depth)
                dose = self.read('alt_dosage')
                self.assertTrue((dose[ploidy == 0] == 0).all())
                # a low bit depth rounds the probabilities heavily, but the
                # samples sharing a ploidy still have to agree with each other
                diploid = dose[ploidy == 2]
                self.assertAlmostEqual(diploid[0], diploid[1], places=5)

    def test_missing_zero_ploidy_sample(self):
        ''' a zero ploidy sample can also be missing rather than empty

        The writer takes an all nan row to mean missing, and a missing sample stores
        as many probabilities as a called one - so this stores none either, and used
        to shift the samples after it in the same way.
        '''
        ploidy = np.array([0, 2, 0, 2], dtype=np.uint8)
        geno = np.full((4, 3), nan)
        geno[ploidy == 2, :3] = [0.2, 0.5, 0.3]
        self.write(ploidy, geno=geno)
        dose = self.read('alt_dosage')
        self.assertTrue(np.isnan(dose[ploidy == 0]).all(), f'{dose}')
        self.assertTrue(np.allclose(dose[ploidy == 2], 1.1, atol=0.01), f'{dose}')

    def test_phased_zero_ploidy(self):
        ''' phased data reads one probability per haplotype, so zero means none

        A phased zero ploidy sample has no cells at all, so an all nan row is the
        only way to write it - which makes it missing.
        '''
        ploidy = np.array([0, 2, 0, 2], dtype=np.uint8)
        geno = np.full((4, 4), nan)
        geno[ploidy == 2] = [0.3, 0.7, 0.3, 0.7]
        self.write(ploidy, geno=geno, phased=True)
        dose = self.read('alt_dosage')
        self.assertTrue(np.isnan(dose[ploidy == 0]).all(), f'{dose}')
        self.assertTrue(np.allclose(dose[ploidy == 2], 1.4, atol=0.01), f'{dose}')

    def test_ploidy_above_two_still_raises(self):
        ''' the ploidy check moved earlier, so make sure it still fires '''
        ploidy = np.array([0, 3, 2], dtype=np.uint8)
        geno = np.full((3, 4), nan)
        geno[0, :1] = [1.0]
        geno[1] = [0.1, 0.2, 0.3, 0.4]
        geno[2, :3] = [0.2, 0.5, 0.3]
        self.write(ploidy, geno=geno)
        with BgenReader(self.path) as bfile:
            var = next(iter(bfile))
            with self.assertRaises(ValueError):
                var.alt_dosage


if __name__ == '__main__':
    unittest.main()
