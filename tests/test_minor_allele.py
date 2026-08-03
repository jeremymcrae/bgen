''' tests for identifying the minor allele and its dosage

The minor allele is picked by estimating the frequency of the first allele from a
sample of the cohort, then comparing it against 0.5. That estimate used to divide
the summed dosage by a fixed batch size of 100 samples times a ploidy of two,
rather than by the alleles it had actually summed, so it was wrong whenever:

  - the cohort was smaller than the batch size, or did not divide evenly by it
  - the samples were not diploid
  - some samples were missing, and so contributed no alleles

Each understates the frequency, which flips the comparison and hands back the major
allele's dosage from minor_allele_dosage.
'''

from pathlib import Path
import tempfile
import unittest

import numpy as np

from bgen import BgenReader, BgenWriter

nan = float('nan')


def diploid_genotypes(n_samples, af0, missing=()):
    ''' genotype probabilities for a biallelic diploid variant

    Args:
        n_samples: number of samples
        af0: frequency of the first allele
        missing: indices of samples to mark as missing

    Returns:
        array of probabilities, ordered as [hom ref, het, hom alt] per sample
    '''
    geno = np.zeros((n_samples, 3))
    n_hom = int(round(n_samples * af0 * af0))
    n_het = int(round(n_samples * 2 * af0 * (1 - af0)))
    geno[:n_hom] = [1, 0, 0]
    geno[n_hom:n_hom + n_het] = [0, 1, 0]
    geno[n_hom + n_het:] = [0, 0, 1]
    if len(missing) > 0:
        geno[list(missing)] = nan
    return geno


def haploid_genotypes(n_samples, af0):
    ''' genotype probabilities for a biallelic haploid variant
    '''
    geno = np.zeros((n_samples, 2))
    n_ref = int(round(n_samples * af0))
    geno[:n_ref] = [1, 0]
    geno[n_ref:] = [0, 1]
    return geno


def expected(geno, ploidy):
    ''' work out the minor allele and its mean dosage from the genotypes

    Checks every sample rather than the sampled estimate the library uses, so it is
    independent of the code under test. Returns (None, None) when both alleles are
    exactly equally frequent, since either answer is then defensible.
    '''
    present = ~np.isnan(geno[:, 0])
    if ploidy == 2:
        dose0 = geno[present, 0] * 2 + geno[present, 1]
    else:
        dose0 = geno[present, 0]
    af0 = dose0.sum() / (ploidy * present.sum())
    if af0 == 0.5:
        return None, None
    minor_is_first = af0 < 0.5
    dose_minor = dose0 if minor_is_first else (ploidy - dose0)
    return ('A' if minor_is_first else 'C'), dose_minor.mean()


class TestMinorAllele(unittest.TestCase):
    ''' check the minor allele and its dosage are identified correctly
    '''

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.tmpdir = Path(self.tmp.name)

    def tearDown(self):
        try:
            self.tmp.cleanup()
        except:
            pass

    def write(self, path, geno, ploidy=2):
        ''' write a single variant bgen and hand back the parsed variant
        '''
        n_samples = len(geno)
        with BgenWriter(path, n_samples,
                        samples=[f's{i}' for i in range(n_samples)]) as bfile:
            bfile.add_variant('var1', 'rs1', 'chr1', 10, ['A', 'C'], geno,
                              ploidy=ploidy)

    def check(self, path, geno, ploidy=2, tolerance=0.01):
        ''' check the reported minor allele and dosage against the genotypes

        Returns False without checking when the alleles are exactly equally
        frequent, since either one is then a defensible answer.
        '''
        want_allele, want_dose = expected(geno, ploidy)
        if want_allele is None:
            return False
        with BgenReader(path) as bfile:
            var = next(iter(bfile))
            # the dosage has to be read before minor_allele is populated
            dose = var.minor_allele_dosage
            self.assertEqual(var.minor_allele, want_allele)
            self.assertAlmostEqual(np.nanmean(dose), want_dose,
                                   delta=tolerance)
        return True

    def test_cohort_smaller_than_the_batch_size(self):
        ''' check the minor allele is right for cohorts under 100 samples

        The estimate divided by a fixed 100 samples however many it summed, so a
        small cohort's frequency was understated tenfold at ten samples.
        '''
        checked = 0
        for n_samples in [3, 10, 37, 50, 99]:
            for af0 in [0.2, 0.6, 0.8, 0.95]:
                geno = diploid_genotypes(n_samples, af0)
                path = self.tmpdir / f'small_{n_samples}_{af0}.bgen'
                self.write(path, geno)
                with self.subTest(n_samples=n_samples, af0=af0):
                    checked += self.check(path, geno)
        # guard against the ties above quietly emptying this test out
        self.assertGreater(checked, 15)

    def test_cohort_not_a_multiple_of_the_batch_size(self):
        ''' check cohorts which don't divide evenly by the batch size

        The stride only covers exactly a batch of samples when the cohort divides
        evenly by it, so these were scaled by the wrong divisor too.
        '''
        checked = 0
        for n_samples in [101, 150, 199, 250, 299]:
            for af0 in [0.3, 0.35, 0.4, 0.45, 0.6]:
                geno = diploid_genotypes(n_samples, af0)
                path = self.tmpdir / f'odd_{n_samples}_{af0}.bgen'
                self.write(path, geno)
                with self.subTest(n_samples=n_samples, af0=af0):
                    checked += self.check(path, geno)
        self.assertGreater(checked, 20)

    def test_haploid_samples(self):
        ''' check the minor allele is right for haploid data

        The frequency assumed every sample was diploid, so haploid data read as
        half its true frequency.
        '''
        checked = 0
        for n_samples in [10, 100, 199, 500]:
            for af0 in [0.2, 0.8, 0.9]:
                geno = haploid_genotypes(n_samples, af0)
                path = self.tmpdir / f'hap_{n_samples}_{af0}.bgen'
                self.write(path, geno, ploidy=1)
                with self.subTest(n_samples=n_samples, af0=af0):
                    checked += self.check(path, geno, ploidy=1)
        self.assertGreater(checked, 10)

    def test_missing_samples(self):
        ''' check missing samples are left out of the frequency

        A missing sample has no called alleles, but its dosage was still zero when
        the frequency was taken, so it pulled the first allele's frequency down.
        '''
        checked = 0
        for n_samples in [100, 199, 500]:
            for af0 in [0.6, 0.9]:
                for fraction in [0.3, 0.5, 0.8]:
                    n_missing = max(int(n_samples * fraction), 1)
                    step = n_samples / n_missing
                    missing = np.unique((np.arange(n_missing) * step).astype(int))
                    geno = diploid_genotypes(n_samples, af0, missing=missing)
                    path = self.tmpdir / f'miss_{n_samples}_{af0}_{fraction}.bgen'
                    self.write(path, geno)
                    with self.subTest(n_samples=n_samples, af0=af0,
                                      fraction=fraction):
                        checked += self.check(path, geno)
        self.assertGreater(checked, 15)

    def test_missing_samples_are_still_nan(self):
        ''' check missing samples come back as nan, and only those samples

        They are marked before the minor allele is chosen now, so check that
        reordering did not lose the nans, or spread them any further.
        '''
        n_samples = 200
        missing = [0, 5, 17, 100, 199]
        geno = diploid_genotypes(n_samples, 0.7, missing=missing)
        path = self.tmpdir / 'nan.bgen'
        self.write(path, geno)

        for prop in ['minor_allele_dosage', 'alt_dosage']:
            with BgenReader(path) as bfile:
                var = next(iter(bfile))
                dose = getattr(var, prop)
                with self.subTest(prop=prop):
                    self.assertEqual(sorted(np.where(np.isnan(dose))[0]),
                                     missing)

    def test_every_sample_missing(self):
        ''' check a variant with no called samples at all doesn't raise

        There are no alleles to take a frequency from, so this only has to come back
        as all nan rather than divide by zero.
        '''
        n_samples = 150
        geno = np.full((n_samples, 3), nan)
        path = self.tmpdir / 'allmiss.bgen'
        self.write(path, geno)
        with BgenReader(path) as bfile:
            var = next(iter(bfile))
            dose = var.minor_allele_dosage
            self.assertTrue(np.isnan(dose).all())

    def test_frequency_either_side_of_one_half(self):
        ''' check the choice is decisive within a single allele of the boundary
        '''
        checked = 0
        for n_samples in [100, 500, 1000]:
            for delta in [-2, -1, 1, 2]:
                # shift the count of the first allele a whole allele either way
                af0 = (n_samples + delta) / (2 * n_samples)
                n_ref = int(round(2 * n_samples * af0))
                geno = np.zeros((n_samples, 3))
                n_hom = n_ref // 2
                n_het = n_ref % 2
                geno[:n_hom] = [1, 0, 0]
                if n_het:
                    geno[n_hom] = [0, 1, 0]
                geno[n_hom + n_het:] = [0, 0, 1]
                path = self.tmpdir / f'edge_{n_samples}_{delta}.bgen'
                self.write(path, geno)
                with self.subTest(n_samples=n_samples, delta=delta):
                    checked += self.check(path, geno)
        # none of these sit on the boundary, so every one must have been checked
        self.assertEqual(checked, 12)

    def test_alt_dosage_is_unaffected(self):
        ''' check alt_dosage still reports the second allele, not the minor one

        alt_dosage doesn't depend on which allele is minor, so it must give the same
        answer whichever way the frequency falls.
        '''
        for n_samples in [10, 199, 500]:
            for af0 in [0.2, 0.8]:
                geno = diploid_genotypes(n_samples, af0)
                path = self.tmpdir / f'alt_{n_samples}_{af0}.bgen'
                self.write(path, geno)
                with BgenReader(path) as bfile:
                    var = next(iter(bfile))
                    dose = var.alt_dosage
                want = (geno[:, 2] * 2 + geno[:, 1]).mean()
                with self.subTest(n_samples=n_samples, af0=af0):
                    self.assertAlmostEqual(np.nanmean(dose), want, delta=0.01)


if __name__ == '__main__':
    unittest.main()
