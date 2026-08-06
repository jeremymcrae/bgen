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

    def test_frequency_the_sampling_cannot_resolve(self):
        ''' check the allele is right when the sampled estimate never settles

        The estimate stops early once a confidence interval clears 0.5, which never
        happens for a frequency near enough to it. The whole cohort then decides, on a
        path the small cohorts above do not reach, so cover it directly. A frequency of
        one allele either side of the boundary is unambiguous, so there is a right
        answer to check against however the frequency is summed.
        '''
        checked = 0
        for n_samples in [500, 999, 1000, 2000]:
            for delta in [-2, -1, 1, 2]:
                n_ref = n_samples + delta
                geno = np.zeros((n_samples, 3))
                n_hom, n_het = n_ref // 2, n_ref % 2
                geno[:n_hom] = [1, 0, 0]
                if n_het:
                    geno[n_hom] = [0, 1, 0]
                geno[n_hom + n_het:] = [0, 0, 1]
                path = self.tmpdir / f'unsettled_{n_samples}_{delta}.bgen'
                self.write(path, geno)
                with self.subTest(n_samples=n_samples, delta=delta):
                    checked += self.check(path, geno)
        self.assertEqual(checked, 16)

    def test_missing_samples_with_a_frequency_near_one_half(self):
        ''' missing samples must be excluded from the whole cohort sum too

        The sampled estimate already skipped them, but a frequency near 0.5 falls
        through to summing every sample, which is a separate piece of code. A missing
        sample there would either poison the total, since its dosage is nan, or be
        counted as a called sample and understate the frequency.

        The frequency is aimed just *below* 0.5, so the first allele is the minor one.
        A poisoned total gives a nan frequency, and nan <= 0.5 is false, so the answer
        comes back as the second allele and the check fails. Aiming above 0.5 would
        hide that, since the broken and the correct answer would agree.
        '''
        checked = 0
        for n_samples in [500, 1000]:
            for n_missing in [1, 7, 32, 100]:
                # spread them out, so they land in different vector lanes
                step = max(n_samples // n_missing, 1)
                missing = list(range(0, n_samples, step))[:n_missing]
                n_called = n_samples - len(missing)
                # two alleles short of half, among the samples actually called
                n_ref = n_called - 2
                geno = np.zeros((n_samples, 3))
                skip = set(missing)
                n_hom, n_het = n_ref // 2, n_ref % 2
                filled = 0
                for i in range(n_samples):
                    if i in skip:
                        continue
                    if filled < n_hom:
                        geno[i] = [1, 0, 0]
                    elif filled == n_hom and n_het:
                        geno[i] = [0, 1, 0]
                    else:
                        geno[i] = [0, 0, 1]
                    filled += 1
                geno[missing] = nan
                path = self.tmpdir / f'nearmiss_{n_samples}_{n_missing}.bgen'
                self.write(path, geno)
                with self.subTest(n_samples=n_samples, n_missing=n_missing):
                    # the first allele must be the minor one for this to have teeth
                    self.assertEqual(expected(geno, 2)[0], 'A')
                    checked += self.check(path, geno)
                    # the missing samples must still be the only nan dosages
                    with BgenReader(path) as bfile:
                        var = next(iter(bfile))
                        dose = var.minor_allele_dosage
                    self.assertEqual(sorted(np.where(np.isnan(dose))[0]),
                                     sorted(missing))
        self.assertEqual(checked, 8)

    def test_samples_with_differing_ploidy(self):
        ''' a variant whose samples differ in ploidy needs its alleles counted per sample

        With one ploidy for the whole cohort the alleles come from the count of called
        samples times that ploidy. That shortcut is wrong as soon as the samples differ,
        and taking it inflates the divisor for haploid samples, so the frequency reads
        low and the wrong allele comes back as minor.

        Half the cohort is haploid here, and the frequency sits far enough above 0.5 that
        the strided estimate settles on it, while treating the haploid samples as diploid
        puts the estimate below a half and flips the answer.

        The cohort has to be large for this to reach the strided sum at all. That estimate
        only returns an answer once a 5 sigma interval clears 0.5, which at a frequency of
        0.57 takes about 1150 samples, and it visits a hundredth of the cohort per stride.
        A few hundred samples would fall through to the whole cohort pass instead and
        never exercise the code this covers.
        '''
        checked = 0
        for n_samples in [5000, 20000]:
            n_hap = n_samples // 2
            n_dip = n_samples - n_hap
            ploidy = np.full(n_samples, 2, dtype=np.uint8)
            ploidy[:n_hap] = 1
            # max_probs for a diploid biallelic sample, haploid rows use the first two
            geno = np.zeros((n_samples, 3))
            # every haploid sample carries the reference allele, so they contribute
            # n_hap alleles rather than the 2 * n_hap the shortcut would assume
            geno[:n_hap, 0] = 1.0
            geno[:n_hap, 2] = nan
            # split the diploid half so the true frequency clears 0.5 while the count
            # the shortcut would use falls below it
            n_ref_hom = int(round(0.36 * n_dip))
            geno[n_hap:n_hap + n_ref_hom] = [1, 0, 0]
            geno[n_hap + n_ref_hom:] = [0, 0, 1]

            n_ref = n_hap + 2 * n_ref_hom
            n_alleles = n_hap + 2 * n_dip
            path = self.tmpdir / f'mixploidy_{n_samples}.bgen'
            with BgenWriter(path, n_samples,
                            samples=[f's{i}' for i in range(n_samples)]) as bfile:
                bfile.add_variant('var1', 'rs1', 'chr1', 10, ['A', 'C'], geno,
                                  ploidy=ploidy, bit_depth=8)

            with self.subTest(n_samples=n_samples):
                # the reference allele is the major one, so 'C' is minor. Counting the
                # haploid samples as diploid would put the frequency below a half and
                # report 'A' instead
                self.assertGreater(n_ref / n_alleles, 0.5)
                self.assertLess(n_ref / (2 * n_samples), 0.5)
                with BgenReader(path) as bfile:
                    var = next(iter(bfile))
                    self.assertEqual(var.minor_allele, 'C')
                checked += 1
        self.assertEqual(checked, 2)


class TestMinorAlleleWithoutDosage(unittest.TestCase):
    ''' the minor allele must not depend on a dosage having been read first

    Which allele is minor is worked out while computing dosages, and used to be left
    in an attribute by that call. Asking for it without reading a dosage returned an
    empty string, so the answer depended on the order the properties were touched.
    '''

    def setUp(self):
        self.tempdir = tempfile.TemporaryDirectory()
        self.tmpdir = Path(self.tempdir.name)

    def tearDown(self):
        self.tempdir.cleanup()

    def write(self, path, geno, ploidy=2, alleles=None):
        n_samples = len(geno)
        with BgenWriter(path, n_samples,
                        samples=[f's{i}' for i in range(n_samples)]) as bfile:
            bfile.add_variant('var1', 'rs1', 'chr1', 10,
                              alleles if alleles else ['A', 'C'], geno,
                              ploidy=ploidy)

    def test_minor_allele_without_reading_a_dosage(self):
        ''' the allele is correct on a freshly opened variant
        '''
        for af0, want in [(0.2, 'A'), (0.8, 'C')]:
            geno = diploid_genotypes(200, af0)
            path = self.tmpdir / f'fresh_{af0}.bgen'
            self.write(path, geno)
            with self.subTest(af0=af0):
                with BgenReader(path) as bfile:
                    var = next(iter(bfile))
                    # nothing else is read from the variant first
                    self.assertEqual(var.minor_allele, want)

    def test_minor_allele_does_not_depend_on_access_order(self):
        ''' every order of reading the properties gives the same allele
        '''
        geno = diploid_genotypes(200, 0.2)
        path = self.tmpdir / 'order.bgen'
        self.write(path, geno)

        def fresh():
            bfile = BgenReader(path)
            return bfile, next(iter(bfile))

        seen = []
        # allele first, then again after each kind of dosage read
        bfile, var = fresh()
        seen.append(var.minor_allele)
        seen.append(var.minor_allele)
        var.alt_dosage
        seen.append(var.minor_allele)
        var.minor_allele_dosage
        seen.append(var.minor_allele)
        bfile.close()
        # and the other way round, dosage before the allele is ever touched
        bfile, var = fresh()
        var.minor_allele_dosage
        seen.append(var.minor_allele)
        bfile.close()
        bfile, var = fresh()
        var.alt_dosage
        seen.append(var.minor_allele)
        bfile.close()
        # reading the probabilities first must not change it either
        bfile, var = fresh()
        var.probabilities
        seen.append(var.minor_allele)
        bfile.close()

        self.assertEqual(len(set(seen)), 1, f'answers varied by order: {seen}')
        self.assertEqual(seen[0], 'A')

    def test_minor_allele_is_never_empty(self):
        ''' a biallelic variant always names one of its alleles
        '''
        for af0 in [0.05, 0.2, 0.5, 0.8, 0.95]:
            geno = diploid_genotypes(150, af0)
            path = self.tmpdir / f'named_{af0}.bgen'
            self.write(path, geno)
            with self.subTest(af0=af0):
                with BgenReader(path) as bfile:
                    var = next(iter(bfile))
                    allele = var.minor_allele
                    self.assertNotEqual(allele, '')
                    self.assertIn(allele, var.alleles)

    def test_minor_allele_across_every_variant_of_a_file(self):
        ''' the allele is right for later variants too, not just the first

        The dosage state is per variant, so a file with several variants checks that
        one variant's answer does not stand in for another's.
        '''
        n_samples = 120
        fracs = [0.1, 0.9, 0.3, 0.7]
        path = self.tmpdir / 'many.bgen'
        with BgenWriter(path, n_samples,
                        samples=[f's{i}' for i in range(n_samples)]) as bfile:
            for i, af0 in enumerate(fracs):
                bfile.add_variant(f'var{i}', f'rs{i}', 'chr1', 10 + i, ['A', 'C'],
                                  diploid_genotypes(n_samples, af0))
        want = ['A' if af0 <= 0.5 else 'C' for af0 in fracs]
        with BgenReader(path) as bfile:
            got = [var.minor_allele for var in bfile]
        self.assertEqual(got, want)

    def test_minor_allele_matches_the_dosage_it_reports(self):
        ''' the named allele is the one minor_allele_dosage counts

        The mean dosage of the minor allele has to be at or below half the ploidy,
        since by definition it is the less common one.
        '''
        for af0 in [0.1, 0.25, 0.75, 0.9]:
            geno = diploid_genotypes(200, af0)
            path = self.tmpdir / f'pair_{af0}.bgen'
            self.write(path, geno)
            with self.subTest(af0=af0):
                with BgenReader(path) as bfile:
                    var = next(iter(bfile))
                    allele = var.minor_allele
                    dose = var.minor_allele_dosage
                self.assertLessEqual(np.nanmean(dose), 1.0)
                self.assertEqual(allele, var.alleles[0] if af0 <= 0.5
                                 else var.alleles[1])

    def test_minor_allele_with_missing_samples(self):
        ''' missing samples must not be counted towards the frequency
        '''
        geno = diploid_genotypes(200, 0.2, missing=[0, 5, 17, 100, 199])
        path = self.tmpdir / 'miss.bgen'
        self.write(path, geno)
        with BgenReader(path) as bfile:
            var = next(iter(bfile))
            self.assertEqual(var.minor_allele, 'A')

    def test_minor_allele_needs_two_alleles(self):
        ''' a multiallelic variant has no single minor allele to report

        There is no dosage for more than two alleles, so this has to say so, rather
        than hand back an empty string as though the variant had no minor allele.
        '''
        n_samples = 10
        # three alleles, so six unphased diploid genotypes per sample
        geno = np.tile([0.5, 0.1, 0.1, 0.1, 0.1, 0.1], (n_samples, 1))
        path = self.tmpdir / 'multi.bgen'
        with BgenWriter(path, n_samples,
                        samples=[f's{i}' for i in range(n_samples)]) as bfile:
            bfile.add_variant('var1', 'rs1', 'chr1', 10, ['A', 'C', 'G'], geno)
        with BgenReader(path) as bfile:
            var = next(iter(bfile))
            with self.assertRaises(ValueError):
                var.minor_allele

    def test_minor_allele_on_a_closed_bgen(self):
        ''' reading the allele after closing must not touch the freed handle
        '''
        geno = diploid_genotypes(50, 0.2)
        path = self.tmpdir / 'closed.bgen'
        self.write(path, geno)
        bfile = BgenReader(path)
        var = next(iter(bfile))
        bfile.close()
        with self.assertRaises(ValueError):
            var.minor_allele

    def test_repeated_reads_agree(self):
        ''' the cached answer matches what a fresh read gives

        The allele is worked out once and kept, so a stale or half updated cache
        would show up as the second read disagreeing with the first.
        '''
        geno = diploid_genotypes(300, 0.35, missing=[1, 2, 3])
        path = self.tmpdir / 'cache.bgen'
        self.write(path, geno)
        with BgenReader(path) as bfile:
            var = next(iter(bfile))
            answers = [var.minor_allele for _ in range(5)]
            var.minor_allele_dosage
            answers += [var.minor_allele for _ in range(5)]
        self.assertEqual(len(set(answers)), 1)
        # and a fresh variant agrees with the cached answer
        with BgenReader(path) as bfile:
            var = next(iter(bfile))
            self.assertEqual(var.minor_allele, answers[0])

    def test_dosage_is_unchanged_by_reading_the_allele_first(self):
        ''' asking for the allele must not disturb the dosages that follow

        Finding the allele computes dosages into a buffer of its own, and for layout 1
        it also spots the missing samples, so a later dosage read has to come back
        with exactly the values it would have given on its own.
        '''
        geno = diploid_genotypes(200, 0.3, missing=[4, 9, 150])
        path = self.tmpdir / 'undisturbed.bgen'
        self.write(path, geno)
        for prop in ['minor_allele_dosage', 'alt_dosage']:
            with self.subTest(prop=prop):
                with BgenReader(path) as bfile:
                    var = next(iter(bfile))
                    plain = np.asarray(getattr(var, prop)).copy()
                with BgenReader(path) as bfile:
                    var = next(iter(bfile))
                    var.minor_allele
                    after = np.asarray(getattr(var, prop)).copy()
                np.testing.assert_array_equal(np.isnan(plain), np.isnan(after))
                both = ~np.isnan(plain)
                np.testing.assert_allclose(plain[both], after[both])


if __name__ == '__main__':
    unittest.main()
