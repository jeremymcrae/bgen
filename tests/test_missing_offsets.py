''' tests for probability offsets of missing samples in phased variable ploidy data

Phased samples occupy one row per haplotype, so the probabilities for a sample start
after the haplotypes of every sample before it. That offset used to be found by
summing the ploidy of all preceding samples once per missing sample, which costs the
whole cohort per missing sample and so grows with the square of the sample count. The
offset is now carried forward between missing samples, since they are visited in
ascending order.

The values written were always correct, so these tests cover the scaling as well as
the offsets themselves.
'''

from pathlib import Path
import tempfile
import time
import unittest

import numpy as np

from bgen import BgenReader, BgenWriter

nan = float('nan')


def phased_genotypes(ploidy, n_alleles=2):
    ''' phased probabilities for a per sample ploidy array

    Each haplotype gets probabilities that identify both the sample and the haplotype,
    so a row read at the wrong offset holds recognisably wrong values rather than
    something that could pass for its neighbour. Haplotypes a sample does not have are
    left as nan, since the array is sized for the largest ploidy present.
    '''
    max_ploidy = int(ploidy.max())
    geno = np.full((len(ploidy), max_ploidy * n_alleles), nan)
    for sample in range(len(ploidy)):
        for hap in range(int(ploidy[sample])):
            # a value that varies by sample and haplotype, kept clear of 0 and 1 so
            # that low bit depths do not collapse neighbouring rows together
            first = 0.1 + 0.8 * (((sample * max_ploidy + hap) % 7) / 7)
            base = hap * n_alleles
            geno[sample, base] = first
            rest = (1 - first) / (n_alleles - 1)
            geno[sample, base + 1:base + n_alleles] = rest
    return geno


def write_phased(path, ploidy, missing, n_alleles=2, bit_depth=8, n_variants=1):
    ''' write a phased bgen with variable ploidy and the given samples missing
    '''
    geno = phased_genotypes(ploidy, n_alleles)
    expected = geno.copy()
    for idx in missing:
        geno[idx, :] = nan
        expected[idx, :] = nan
    n_samples = len(ploidy)
    with BgenWriter(path, n_samples,
                    samples=[f's{i}' for i in range(n_samples)]) as bfile:
        for i in range(n_variants):
            bfile.add_variant(f'var{i}', f'rs{i}', 'chr1', 10 + i,
                              [c for c in 'ACGT'[:n_alleles]], geno,
                              ploidy=ploidy, phased=1, bit_depth=bit_depth)
    return expected


class TestMissingSampleOffsets(unittest.TestCase):
    ''' check missing samples are marked at the right offsets, and scale linearly
    '''

    def setUp(self):
        self.tempdir = tempfile.TemporaryDirectory()
        self.tmpdir = Path(self.tempdir.name)

    def tearDown(self):
        self.tempdir.cleanup()

    def assert_probs_match(self, path, expected, bit_depth=8):
        ''' every sample matches, and exactly the missing ones are nan
        '''
        # probabilities come back as float32, so at high bit depths the encoding is
        # finer than the type can represent and that sets the floor on the tolerance
        tol = max(2.0 / (2 ** bit_depth - 1), 1e-6)
        bfile = BgenReader(path, delay_parsing=True)
        n_variants = 0
        for var in bfile:
            probs = var.probabilities
            n_variants += 1
            self.assertEqual(probs.shape, expected.shape)
            # compare nan placement separately, since nan != nan
            np.testing.assert_array_equal(np.isnan(probs), np.isnan(expected))
            both = ~np.isnan(expected)
            self.assertTrue(np.all(np.abs(probs[both] - expected[both]) <= tol))
        bfile.close()
        self.assertGreater(n_variants, 0)

    def test_missing_samples_at_edges(self):
        ''' the first and last sample are the offsets most easily got wrong
        '''
        ploidy = np.array([1, 2, 1, 2, 1, 2, 1, 2], dtype=np.uint8)
        for missing in [[0], [7], [0, 7], [0, 1], [6, 7]]:
            with self.subTest(missing=missing):
                path = self.tmpdir / 'edges.bgen'
                expected = write_phased(path, ploidy, missing)
                self.assert_probs_match(path, expected)

    def test_consecutive_missing_samples(self):
        ''' neighbouring missing samples share no work, so the offsets must still move
        '''
        ploidy = np.array([2, 1, 3, 1, 2, 3, 1, 2, 3, 1], dtype=np.uint8)
        path = self.tmpdir / 'consecutive.bgen'
        expected = write_phased(path, ploidy, [2, 3, 4, 5])
        self.assert_probs_match(path, expected)

    def test_all_and_no_samples_missing(self):
        ''' the two extremes of the missing list
        '''
        ploidy = np.array([1, 2, 3, 2, 1], dtype=np.uint8)
        for missing in [[], list(range(5))]:
            with self.subTest(n_missing=len(missing)):
                path = self.tmpdir / 'extremes.bgen'
                expected = write_phased(path, ploidy, missing)
                self.assert_probs_match(path, expected)

    def test_missing_offsets_across_ploidies_and_depths(self):
        ''' the offset scales by ploidy and the values by bit depth, so vary both
        '''
        for ploidy_vals in [[1, 2], [2, 3], [1, 2, 3, 4], [3, 3]]:
            for bit_depth in [1, 2, 8, 16, 32]:
                for n_alleles in [2, 3]:
                    with self.subTest(ploidy=ploidy_vals, bit_depth=bit_depth,
                                      n_alleles=n_alleles):
                        # repeat the pattern so missing samples fall on varying ploidy
                        ploidy = np.array((ploidy_vals * 4)[:8], dtype=np.uint8)
                        path = self.tmpdir / 'varied.bgen'
                        expected = write_phased(path, ploidy, [1, 4, 7],
                                                n_alleles=n_alleles,
                                                bit_depth=bit_depth)
                        self.assert_probs_match(path, expected, bit_depth=bit_depth)

    def test_multiple_variants_reuse_offsets(self):
        ''' the running total is per variant, so a second variant must not inherit it
        '''
        ploidy = np.array([1, 2, 3, 1, 2, 3], dtype=np.uint8)
        path = self.tmpdir / 'multi.bgen'
        expected = write_phased(path, ploidy, [0, 3, 5], n_variants=4)
        self.assert_probs_match(path, expected)

    def _time_read(self, n_samples):
        ''' best time to read every probability of an all missing cohort

        All samples missing is the worst case, since the offset is recomputed for
        each one. Ploidy alternates so the variable ploidy branch is taken.
        '''
        ploidy = np.resize(np.array([1, 2], dtype=np.uint8), n_samples)
        path = self.tmpdir / f'scale{n_samples}.bgen'
        write_phased(path, ploidy, range(n_samples), n_variants=2)
        best = None
        for _ in range(3):
            bfile = BgenReader(path, delay_parsing=True)
            start = time.perf_counter()
            for var in bfile:
                var.probabilities
            elapsed = time.perf_counter() - start
            bfile.close()
            best = elapsed if best is None else min(best, elapsed)
        path.unlink()
        return best / n_samples

    def test_missing_offsets_scale_linearly(self):
        ''' cost per sample must not grow with the size of the cohort

        Re-summing the preceding ploidy values for each missing sample costs the
        whole cohort per missing sample, so the time per sample grows in step with
        the sample count, while carrying the total forward leaves it flat.

        The cohorts have to be large enough for that growth to outweigh the fixed
        cost of parsing each sample, which does not change with cohort size. At
        these sizes the re-summing version measured a growth of 2.4 to 2.9, and the
        carried version 1.0, so the threshold sits between the two with roughly a
        1.5x margin either way.
        '''
        small = self._time_read(20000)
        large = self._time_read(80000)
        growth = large / small
        self.assertLess(growth, 1.6,
                        f'per sample cost grew {growth:.1f}x when the cohort grew 4x, '
                        'which suggests the missing sample offsets are being '
                        're-summed rather than carried forward')


if __name__ == '__main__':
    unittest.main()
