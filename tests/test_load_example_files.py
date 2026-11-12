
from pathlib import Path
import tempfile
import unittest


from bgen import BgenReader, BgenWriter

from tests.utils import load_gen_data, load_vcf_data, load_haps_data, arrays_equal

class TestExampleBgens(unittest.TestCase):
    ''' class to make sure we can load bgen files
    '''
    @classmethod
    def setUpClass(cls):
        cls.gen_data = load_gen_data()
        cls.vcf_data = load_vcf_data()
        cls.haps_data = load_haps_data()
    
    def setUp(self):
        ''' set path to folder with test data
        '''
        self.folder = Path(__file__).parent /  "data"
        # files written by the tests go here, rather than among the test data
        self.tmp = tempfile.TemporaryDirectory()
        self.tmpdir = Path(self.tmp.name)
        self.addCleanup(self.tmp.cleanup)
    
    def test_load_example_genotypes_bit_depths(self):
        ''' check parsing genotypes from the example files with different bit depths
        '''
        for path in self.folder.glob('example.*bits.bgen'):
            bit_depth = int(path.stem.split('.')[1].strip('bits'))
            with BgenReader(str(path)) as bfile:
                for var, g in zip(bfile, self.gen_data):
                    self.assertEqual(g, var)
                    self.assertTrue(arrays_equal(g.probabilities, var.probabilities, bit_depth))
    
    def test_zstd_compressed(self):
        ''' check we can parse genotypes from zstd compressed geno probabilities
        '''
        path = self.folder / 'example.16bits.zstd.bgen'
        with BgenReader(str(path)) as bfile:
            for var, g in zip(bfile, self.gen_data):
                self.assertEqual(g, var)
                self.assertTrue(arrays_equal(g.probabilities, var.probabilities, 16))
    
    def test_v11(self):
        ''' check we can open a bgen in v1.1 format, and parse genotypes correctly
        '''
        path = self.folder / 'example.v11.bgen'
        with BgenReader(str(path)) as bfile:
            bit_depth = 16
            for var, g in zip(bfile, self.gen_data):
                self.assertEqual(g, var)
                self.assertTrue(arrays_equal(g.probabilities, var.probabilities, bit_depth))
    
    def test_Path(self):
        ''' check we can open bgen files from Path objects
        '''
        path = self.folder / 'example.v11.bgen'
        with BgenReader(path):
            pass
    
    def test_load_haplotypes_bgen(self):
        ''' check we can open a bgen with haplotypes, and parse genotypes correctly
        '''
        path = self.folder / 'haplotypes.bgen'
        with BgenReader(str(path)) as bfile:
            bit_depth = 16
            for var, g in zip(bfile, self.haps_data):
                self.assertEqual(g, var)
                self.assertTrue(arrays_equal(g.probabilities, var.probabilities, bit_depth))
    
    def test_load_complex_file(self):
        ''' make sure we can open a complex bgen file
        '''
        path = self.folder / 'complex.bgen'
        with BgenReader(path) as bfile:
            bit_depth = 16
            for var, g in zip(bfile, self.vcf_data):
                self.assertEqual(g, var)
                self.assertTrue(arrays_equal(g.probabilities, var.probabilities, bit_depth))
                self.assertTrue(all(x == y for x, y in zip(g.ploidy, var.ploidy)))
    
    def test_load_complex_files(self):
        ''' make sure we can open the complex bgen files
        '''
        
        for path in self.folder.glob('complex.*.bgen'):
            bit_depth = int(path.stem.split('.')[1].strip('bits'))
            with BgenReader(path) as bfile:
                for var, g in zip(bfile, self.vcf_data):
                    self.assertEqual(g, var)
                    self.assertTrue(arrays_equal(g.probabilities, var.probabilities, bit_depth))
    
    def test_load_sample_file(self):
        ''' make sure we can open a sample file for a bgen
        '''
        bgen_path = self.folder / 'complex.bgen'
        sample_path = self.folder / 'complex.sample'
        orig = BgenReader(bgen_path)
        orig_samples = orig.samples
        
        # construct a bgen file without any sample IDs inside (otherwise 
        # the internal IDs take priority). The variant data is copied across as
        # raw bytes, so the layout and compression have to match the source
        bgen_path = self.tmpdir / 'temp.bgen'
        with BgenWriter(bgen_path, n_samples=len(orig.samples),
                        layout=orig.header.layout,
                        compression=orig.header.compression) as bfile:
            for var in orig:
                bfile.add_variant_direct(var)
        orig.close()
        
        # check bgen files without internal IDs or an external file instead
        # use numeric IDs (converted to strings)
        with BgenReader(bgen_path) as bfile:
            numeric_ids = [f'{x}' for x in range(len(orig_samples))]
            self.assertEqual(numeric_ids, bfile.samples)
        
        # reading sample IDs from the corresponding sample file should give
        # identical IDs
        with BgenReader(bgen_path, sample_path) as bfile:
            self.assertEqual(orig_samples, bfile.samples)
        
        # check we raise an error with too few sample IDs
        missing_path = self.tmpdir / 'empty.sample'
        missing = open(missing_path, 'wt')
        missing.close()
        
        with self.assertRaises(ValueError):
            BgenReader(bgen_path, missing_path)
        
        # check we raise an error with too many sample IDs
        extra_path = self.tmpdir / 'extra.sample'
        with open(extra_path, 'wt') as extra:
            extra.write('id\n0\nsample_0\nsample_1\nsample_2\nsample_3\nsample_4\n')
        
        with self.assertRaises(ValueError):
            BgenReader(bgen_path, extra_path)
    
    def test_load_sample_file_delimiters(self):
        ''' sample files can be tab delimited, or have windows line endings
        '''
        bgen_path = self.tmpdir / 'temp.bgen'
        expected = ['sample_0', 'sample_1', 'sample_2', 'sample_3']
        
        # build a bgen without internal sample IDs, so the sample file is used
        orig = BgenReader(self.folder / 'complex.bgen')
        with BgenWriter(bgen_path, n_samples=len(orig.samples),
                        layout=orig.header.layout,
                        compression=orig.header.compression) as bfile:
            for var in orig:
                bfile.add_variant_direct(var)
        orig.close()
        
        rows = ['id', '0'] + expected
        for label, text in [
                ('space delimited', ''.join(f'{x} 0\n' for x in rows)),
                ('tab delimited', ''.join(f'{x}\t0\n' for x in rows)),
                ('windows line endings', ''.join(f'{x}\r\n' for x in rows)),
                ('no trailing newline', '\n'.join(rows)),
                # read in binary, a blank windows line is a lone carriage return, so
                # it has to be recognised as empty rather than taken as an ID
                ('windows blank last line', ''.join(f'{x}\r\n' for x in rows) + '\r\n'),
                ('blank last line', ''.join(f'{x}\n' for x in rows) + '\n'),
                ]:
            sample_path = self.tmpdir / 'delim.sample'
            with open(sample_path, 'w', newline='') as handle:
                handle.write(text)
            with BgenReader(bgen_path, sample_path) as bfile:
                self.assertEqual(bfile.samples, expected, label)
            sample_path.unlink()
    
    def test_load_missing_file(self):
        ''' check passing in a path to a missing file fails gracefully
        '''
        with self.assertRaises(ValueError):
            BgenReader('/zzz/jjj/qqq.bgen')
    
    def test_load_missing_sample_file(self):
        ''' a missing sample file fails gracefully, for a bgen without its own IDs
        '''
        path = self.folder / 'example.v11.bgen'
        with self.assertRaises(ValueError):
            BgenReader(path, '/zzz/jjj/qqq.sample')
