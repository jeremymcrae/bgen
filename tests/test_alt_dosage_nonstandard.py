

from pathlib import Path
import tempfile
import unittest

import numpy as np

from bgen import BgenReader, BgenWriter

def get_alt_dose(variant):
    ''' compute alt dosage from variant probabilities
    
    This accounts for ploidy and phased status.
    '''
    probs = variant.probabilities
    haploid = variant.ploidy == 1
    
    alt_dose = np.empty(len(probs))
    if variant.is_phased:
        alt_dose[~haploid] = probs[~haploid, 1] + probs[~haploid, 3]
        alt_dose[haploid] = probs[haploid, 1]
    else:
        alt_dose[~haploid] = 2 * probs[~haploid, 2] + probs[~haploid, 1]
        alt_dose[haploid] = probs[haploid, 1]
    
    return alt_dose

class TestAltDosage(unittest.TestCase):
    ''' class to make sure alt dosage is computed correctly
    '''
    
    def setUp(self):
        ''' set path to folder with test data
        '''
        self.folder = Path(__file__).parent /  "data"
    
    def test_alt_dosage_nonstandard(self):
        ''' variant.alt_dosage is correct with variable ploidy and with phased data
        '''
        path = self.folder / 'alt_dosage_check.bgen'
        with BgenReader(path) as bfile:
            for variant in bfile:
                dose = variant.alt_dosage
                self.assertTrue((dose >= 0).all())
                self.assertTrue((dose == get_alt_dose(variant)).all())

    def test_alt_dosage_at_every_bit_depth(self):
        ''' dosage must match the probabilities at all bit depths

        A probability of one fills a 32 bit integer at the maximum bit depth, so
        doubling it for a diploid genotype used to wrap, turning a dosage of two
        into one. Dosages below one were unaffected, hence the range checked here.
        '''
        geno = np.array([[1.0, 0.0, 0.0],
                         [0.0, 1.0, 0.0],
                         [0.0, 0.0, 1.0],
                         [0.5, 0.5, 0.0],
                         [0.0, 0.5, 0.5],
                         [0.25, 0.5, 0.25],
                         ])
        n = len(geno)
        with tempfile.TemporaryDirectory() as folder:
            for bit_depth in range(1, 33):
                path = Path(folder) / f'depth_{bit_depth}.bgen'
                with BgenWriter(path, n, samples=[f's{i}' for i in range(n)]) as bfile:
                    bfile.add_variant('var1', 'rs1', 'chr1', 10, ['A', 'C'], geno,
                                      bit_depth=bit_depth)

                # compare against the probabilities read back from the same file,
                # so the rounding the bit depth imposes is already accounted for
                with BgenReader(path) as bfile:
                    expected = get_alt_dose(next(iter(bfile)))
                with BgenReader(path) as bfile:
                    dose = next(iter(bfile)).alt_dosage

                with self.subTest(bit_depth=bit_depth):
                    np.testing.assert_allclose(
                        dose, expected, atol=1e-6,
                        err_msg=f'alt dosage at bit_depth={bit_depth}')
                    self.assertTrue((dose >= 0).all())
                    self.assertTrue((dose <= 2).all())

    def test_alt_dosage_in_32_bit_example_file(self):
        ''' every variant of the shipped 32 bit example file had wrong dosages '''
        path = self.folder / 'example.32bits.bgen'
        with BgenReader(path) as bfile:
            for variant in bfile:
                dose = variant.alt_dosage
                expected = get_alt_dose(variant)
                np.testing.assert_allclose(dose, expected, atol=1e-6,
                                           err_msg=f'{variant.rsid}')

    def test_minor_allele_at_the_maximum_bit_depth(self):
        ''' the wrapped dosage also picked the wrong minor allele

        The first allele's frequency is estimated from the dosages, so wrapping
        them flipped the comparison against one half.
        '''
        n = 40
        # make the first allele clearly the major one, so the second is minor
        geno = np.zeros((n, 3))
        geno[:, 0] = 1.0
        geno[:4] = [0.0, 1.0, 0.0]
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / 'major_first.bgen'
            with BgenWriter(path, n, samples=[f's{i}' for i in range(n)]) as bfile:
                bfile.add_variant('var1', 'rs1', 'chr1', 10, ['A', 'C'], geno,
                                  bit_depth=32)

            with BgenReader(path) as bfile:
                var = next(iter(bfile))
                minor_dose = var.minor_allele_dosage
                self.assertEqual(var.minor_allele, 'C')
                np.testing.assert_allclose(minor_dose, get_alt_dose(var), atol=1e-6)
