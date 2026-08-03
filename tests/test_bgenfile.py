
from pathlib import Path
import struct
import subprocess
import sys
import tempfile
import textwrap
import unittest
import warnings

import numpy as np

from bgen import BgenReader

from tests.utils import load_gen_data

try:
    import resource
    HAS_RLIMIT = hasattr(resource, 'RLIMIT_AS')
except ImportError:
    # windows has no resource module, so the memory capped tests are skipped there
    HAS_RLIMIT = False

# ceiling for the memory capped tests. Comfortably fits the interpreter, numpy and a
# sub-megabyte bgen, but not one string per sample for a corrupt sample count
MEMORY_CAP = 2 * 1024 ** 3

def run_capped(code):
    ''' run code in a subprocess under a memory cap
    
    The cap has to apply to the whole process, so this cannot run in-process. The
    child inherits this process' sys.path, so it imports the same build under test
    rather than any installed copy.
    '''
    preamble = textwrap.dedent(f'''
        import resource, sys
        sys.path[:] = {sys.path!r}
        resource.setrlimit(resource.RLIMIT_AS, ({MEMORY_CAP}, {MEMORY_CAP}))
        ''')
    return subprocess.run([sys.executable, '-c', preamble + textwrap.dedent(code)],
                          capture_output=True, text=True)

class TestBgenReader(unittest.TestCase):
    ''' class to make sure BgenReader works correctly
    '''
    
    @classmethod
    def setUpClass(cls):
        cls.gen_data = load_gen_data()
    
    def setUp(self):
        ''' set path to folder with test data
        '''
        self.folder = Path(__file__).parent /  "data"
        self._tempdir = tempfile.TemporaryDirectory()
        self.tmpdir = self._tempdir.name
        self.addCleanup(self._tempdir.cleanup)
    
    def test_context_handler_closed_bgen_samples(self):
        ''' no samples available from exited BgenReader
        '''
        path = self.folder / 'example.16bits.zstd.bgen'
        with BgenReader(path) as bfile:
            self.assertTrue(len(bfile.samples) > 0)
        
        with self.assertRaises(ValueError):
            bfile.samples
    
    def test_context_handler_closed_bgen_varids(self):
        ''' no varids available from exited BgenReader
        '''
        path = self.folder / 'example.16bits.zstd.bgen'
        with BgenReader(path) as bfile:
            self.assertTrue(len(bfile.varids()) > 0)
        
        with self.assertRaises(ValueError):
            bfile.varids()
    
    def test_context_handler_closed_bgen_rsids(self):
        ''' no rsids available from exited BgenReader
        '''
        path = self.folder / 'example.16bits.zstd.bgen'
        with BgenReader(path) as bfile:
            self.assertTrue(len(bfile.rsids()) > 0)
        
        with self.assertRaises(ValueError):
            bfile.rsids()
    
    def test_context_handler_closed_bgen_positions(self):
        ''' no positions available from exited BgenReader
        '''
        path = self.folder / 'example.16bits.zstd.bgen'
        with BgenReader(path) as bfile:
            self.assertTrue(len(bfile.positions()) > 0)
        
        with self.assertRaises(ValueError):
            bfile.positions()
    
    def test_context_handler_closed_bgen_length(self):
        ''' error raised if accessing length of exited BgenReader
        '''
        path = self.folder / 'example.16bits.zstd.bgen'
        with BgenReader(path) as bfile:
            self.assertTrue(len(bfile) > 0)
        
        with self.assertRaises(ValueError):
             len(bfile)
    
    def test_context_handler_closed_bgen_slice(self):
        ''' error raised if slicing variant from exited BgenReader
        '''
        path = self.folder / 'example.16bits.zstd.bgen'
        with BgenReader(path) as bfile:
            self.assertTrue(len(bfile) > 0)
        
        with self.assertRaises(ValueError):
             var = bfile[0]
    
    def test_context_handler_closed_bgen_at_position(self):
        ''' error raised if getting variant at position from exited BgenReader
        '''
        path = self.folder / 'example.16bits.zstd.bgen'
        with BgenReader(path) as bfile:
            self.assertTrue(len(bfile) > 0)
        
        with self.assertRaises(ValueError):
             var = bfile.at_position(100)
    
    def test_context_handler_closed_bgen_with_rsid(self):
        ''' error raised if getting variant with rsid from exited BgenReader
        '''
        path = self.folder / 'example.16bits.zstd.bgen'
        with BgenReader(path) as bfile:
            self.assertTrue(len(bfile) > 0)
        
        with self.assertRaises(ValueError):
             var = bfile.with_rsid('rs111')
    
    def test_context_handler_variant_data_not_loaded(self):
        ''' error raised if we try to access variant data after closing BgenFile
        '''
        path = self.folder / 'example.16bits.zstd.bgen'
        with BgenReader(path) as bfile:
            var = next(bfile)
        
        with self.assertRaises(ValueError):
            # cannot access data after file closed
            var.minor_allele_dosage
    
    def test_context_handler_variant_data_loaded(self):
        '''no error raised for variant from closed BgenReader, IF data is already loaded
        '''
        path = self.folder / 'example.16bits.zstd.bgen'
        with BgenReader(path) as bfile:
            var = next(bfile)
            var.minor_allele_dosage # load data while file still open
        
        # can access data after file closed, but only because the file was read 
        # previously while still open
        dose = var.minor_allele_dosage
        self.assertTrue(isinstance(dose, np.ndarray))
    
    def test_fetch(self):
        ''' can fetch variants within a genomic region
        '''
        chrom, start, stop = '01', 5000, 50000
        bfile = BgenReader(self.folder / 'example.16bits.bgen')
        self.assertTrue(bfile._check_for_index(str(self.folder / 'example.16bits.bgen')))
        
        self.assertTrue(list(bfile.fetch('02')) == [])
    
    def test_fetch_whole_chrom(self):
        ''' fetching just with chrom gives all variants on chromosome
        '''
        chrom, start, stop = '01', 5000, 50000
        bfile = BgenReader(self.folder / 'example.16bits.bgen')
        
        # test fetching a whole chromosome
        sortkey = lambda x: (x.chrom, x.pos)
        for x, y in zip(sorted(bfile.fetch(chrom), key=sortkey), sorted(self.gen_data, key=sortkey)):
            self.assertEqual(x.rsid, y.rsid)
            self.assertEqual(x.chrom, y.chrom)
            self.assertEqual(x.pos, y.pos)
    
    def test_fetch_after_position(self):
        ''' fetching variants with chrom and start gives all variants after pos
        '''
        chrom, start, stop = '01', 5000, 50000
        bfile = BgenReader(self.folder / 'example.16bits.bgen')
        
        sortkey = lambda x: (x.chrom, x.pos)
        gen_vars = [x for x in sorted(self.gen_data, key=sortkey) if start <= x.pos]
        for x, y in zip(sorted(bfile.fetch(chrom, start), key=sortkey), gen_vars):
            self.assertEqual(x.rsid, y.rsid)
            self.assertEqual(x.chrom, y.chrom)
            self.assertEqual(x.pos, y.pos)
    
    def test_fetch_in_region(self):
        ''' fetching variants with chrom, start, stop gives variants in region
        '''
        chrom, start, stop = '01', 5000, 50000
        bfile = BgenReader(self.folder / 'example.16bits.bgen')
        
        sortkey = lambda x: (x.chrom, x.pos)
        gen_vars = [x for x in sorted(self.gen_data, key=sortkey) if start <= x.pos <= stop]
        for x, y in zip(sorted(bfile.fetch(chrom, start, stop), key=sortkey), gen_vars):
            self.assertEqual(x.rsid, y.rsid)
            self.assertEqual(x.chrom, y.chrom)
            self.assertEqual(x.pos, y.pos)
        
        # check that we don't get any variants in a region without any
        self.assertEqual(list(bfile.fetch(chrom, start * 1000, stop * 1000)), [])
    
    def test_corrupt_sample_count(self):
        ''' an impossible sample count is rejected rather than sizing an allocation
        
        nsamples comes straight from the header, so a corrupt count used to be
        handed to resize(), which asked for over a hundred GB from a tiny file.
        '''
        data = bytearray((self.folder / 'example.16bits.bgen').read_bytes())
        # nsamples is the 4 byte field at offset 12
        struct.pack_into('<I', data, 12, 0xFFFFFFFF)
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / 'corrupt.bgen'
            path.write_bytes(bytes(data))
            with self.assertRaises(ValueError):
                BgenReader(path, delay_parsing=True)
    
    def test_corrupt_sample_count_without_sample_ids(self):
        ''' a corrupt sample count does not allocate when the IDs are not in the bgen
        
        The file size check only applies to the sample block, so with the sample ID
        flag cleared there is nothing in the bgen to bound the count. Opening used
        to build one string per claimed sample, so a tiny file cost gigabytes.
        Reading genotypes is what detects the bad count.
        '''
        path = self._corrupt_count_bgen(self.tmpdir, 100_000_000)
        bfile = BgenReader(path, delay_parsing=True)
        self.assertEqual(bfile.header.nsamples, 100_000_000)
        var = next(iter(bfile))
        with self.assertRaises(ValueError):
            var.probabilities
        bfile.close()
    
    @unittest.skipUnless(HAS_RLIMIT, 'needs RLIMIT_AS to cap memory')
    def test_opening_bgen_without_sample_ids_does_not_allocate(self):
        ''' opening must not allocate memory proportional to the sample count
        
        Runs under a memory cap in a subprocess, since the point is the size of the
        allocation rather than the result. Building 100 million IDs needs several GB,
        so the cap fails the old behaviour while leaving ample room for the file
        itself, which is under a megabyte.
        '''
        path = self._corrupt_count_bgen(self.tmpdir, 100_000_000)
        out = run_capped(f'''
            from bgen import BgenReader
            bfile = BgenReader({str(path)!r}, delay_parsing=True)
            assert bfile.header.nsamples == 100_000_000
            bfile.close()
            print('opened')
            ''')
        self.assertEqual(out.returncode, 0, msg=out.stderr[-2000:])
        self.assertIn('opened', out.stdout)
    
    def _corrupt_count_bgen(self, folder, nsamples):
        ''' copy a bgen, overstate its sample count and clear its sample ID flag
        '''
        data = bytearray((self.folder / 'example.16bits.bgen').read_bytes())
        # nsamples is the 4 byte field at offset 12
        struct.pack_into('<I', data, 12, nsamples)
        # clear the has_sample_ids flag, the top bit of the flags field that ends
        # the header
        header_length = struct.unpack_from('<I', data, 4)[0]
        flags = struct.unpack_from('<I', data, header_length)[0]
        struct.pack_into('<I', data, header_length, flags & ~(1 << 31))
        path = Path(folder) / f'corrupt.{nsamples}.bgen'
        path.write_bytes(bytes(data))
        return path
    
    def test_corrupt_sample_count_with_sample_file(self):
        ''' a corrupt sample count does not allocate for an external sample file
        
        The IDs come from the .sample file, so the bgen cannot bound the count. The
        list used to be sized from the count before reading, so it allocated for
        every claimed sample before finding the file only held a few.
        '''
        path = self._corrupt_count_bgen(self.tmpdir, 100_000_000)
        sample_path = Path(self.tmpdir) / 'few.sample'
        sample_path.write_text('ID_1 ID_2 missing\n0 0 0\nsample_1 sample_1 0\n')
        with self.assertRaises(ValueError):
            BgenReader(path, sample_path=str(sample_path), delay_parsing=True)
    
    @unittest.skipUnless(HAS_RLIMIT, 'needs RLIMIT_AS to cap memory')
    def test_sample_file_mismatch_does_not_allocate(self):
        ''' a sample file shorter than the claimed count must not allocate first
        '''
        path = self._corrupt_count_bgen(self.tmpdir, 100_000_000)
        sample_path = Path(self.tmpdir) / 'few2.sample'
        sample_path.write_text('ID_1 ID_2 missing\n0 0 0\nsample_1 sample_1 0\n')
        code = f'''
            from bgen import BgenReader
            try:
                BgenReader({str(path)!r}, sample_path={str(sample_path)!r},
                           delay_parsing=True)
            except ValueError:
                print('rejected')
            '''
        out = run_capped(code)
        self.assertEqual(out.returncode, 0, msg=out.stderr[-2000:])
        self.assertIn('rejected', out.stdout)
    
    def test_generated_sample_ids(self):
        ''' a bgen without sample IDs still gets sequential placeholder IDs
        
        These are built on demand now, so check they are unchanged, and that asking
        twice gives the same answer.
        '''
        path = self.folder / 'example.v11.bgen'
        bfile = BgenReader(path, delay_parsing=True)
        samples = bfile.samples
        self.assertEqual(len(samples), 500)
        self.assertEqual(samples[:3], ['0', '1', '2'])
        self.assertEqual(samples[-1], '499')
        self.assertEqual(bfile.samples, samples)
        bfile.close()
    
    def test_sample_ids_from_sample_file(self):
        ''' IDs from an external sample file are read in order
        '''
        data = bytearray((self.folder / 'example.16bits.bgen').read_bytes())
        header_length = struct.unpack_from('<I', data, 4)[0]
        flags = struct.unpack_from('<I', data, header_length)[0]
        struct.pack_into('<I', data, header_length, flags & ~(1 << 31))
        struct.pack_into('<I', data, 12, 3)
        path = Path(self.tmpdir) / 'nosamples.bgen'
        path.write_bytes(bytes(data))
        sample_path = Path(self.tmpdir) / 'three.sample'
        sample_path.write_text('ID_1 ID_2 missing\n0 0 0\n'
                               'first first 0\nsecond second 0\nthird third 0\n')
        bfile = BgenReader(path, sample_path=str(sample_path), delay_parsing=True)
        self.assertEqual(bfile.samples, ['first', 'second', 'third'])
        bfile.close()
    
    def test_sample_block_length_mismatch(self):
        ''' a sample block whose length disagrees with its IDs is rejected
        
        The variants are found through the header offset, not by walking the sample
        block, so a wrong block length used to pass unnoticed even though the file is
        internally inconsistent.
        '''
        for delta in [100, -100, 2, -2]:
            with self.subTest(delta=delta):
                data = bytearray((self.folder / 'example.16bits.bgen').read_bytes())
                # the sample block starts after the header, and opens with its own
                # length
                header_length = struct.unpack_from('<I', data, 4)[0]
                pos = 4 + header_length
                block_length = struct.unpack_from('<I', data, pos)[0]
                struct.pack_into('<I', data, pos, block_length + delta)
                path = Path(self.tmpdir) / f'block{delta}.bgen'
                path.write_bytes(bytes(data))
                with self.assertRaises(ValueError):
                    BgenReader(path, delay_parsing=True)
    
    def test_sample_count_beyond_block_capacity(self):
        ''' a count larger than the sample block could hold is rejected up front
        
        The block length caps how many IDs can follow it, since each carries a two
        byte length prefix. A count above that cannot be satisfied whatever the block
        contains, so it is rejected without reading the IDs.
        '''
        data = bytearray((self.folder / 'example.16bits.bgen').read_bytes())
        header_length = struct.unpack_from('<I', data, 4)[0]
        block_at = 4 + header_length
        block_length = struct.unpack_from('<I', data, block_at)[0]
        capacity = (block_length - 8) // 2
        # keep the count in the block matching the header, so the mismatch under test
        # is the capacity rather than the two counts disagreeing
        for claimed in [capacity + 1, capacity * 2, 0xFFFFFFFF]:
            with self.subTest(claimed=claimed):
                struct.pack_into('<I', data, 12, claimed)
                struct.pack_into('<I', data, block_at + 4, claimed)
                path = Path(self.tmpdir) / f'cap{claimed}.bgen'
                path.write_bytes(bytes(data))
                with self.assertRaises(ValueError):
                    BgenReader(path, delay_parsing=True)
    
    def test_sample_block_length_counts_ids_exactly(self):
        ''' the block length must match the IDs to the byte
        
        The check adds two bytes of length prefix per ID to a running total, so it has
        to cope with IDs that are empty, long, or contain spaces, and be off by one in
        neither direction.
        '''
        cases = {'plain': [b'a', b'bb', b'ccc'],
                 'empty id': [b'a', b'', b'c'],
                 'spaces': [b'first last', b'x'],
                 'long id': [b'z' * 1000, b'y'],
                 'single': [b'only']}
        for label, ids in cases.items():
            for delta in [0, 1, -1]:
                with self.subTest(ids=label, delta=delta):
                    path = self._bgen_with_ids(ids, block_delta=delta)
                    if delta == 0:
                        bfile = BgenReader(path, delay_parsing=True)
                        self.assertEqual(bfile.samples,
                                         [x.decode('utf8') for x in ids])
                        bfile.close()
                    else:
                        with self.assertRaises(ValueError):
                            BgenReader(path, delay_parsing=True)
    
    def _bgen_with_ids(self, ids, block_delta=0):
        ''' write a bgen holding exactly these sample IDs and no variants
        '''
        blob = b''.join(struct.pack('<H', len(x)) + x for x in ids)
        block = struct.pack('<II', 8 + len(blob) + block_delta, len(ids)) + blob
        header_length = 20
        flags = (2 << 2) | (1 << 31)
        data = struct.pack('<IIII', header_length + len(block), header_length, 0,
                           len(ids))
        data += b'bgen' + struct.pack('<I', flags) + block
        path = Path(self.tmpdir) / f'ids{len(ids)}_{block_delta}_{len(blob)}.bgen'
        path.write_bytes(data)
        return path
    
    def test_valid_sample_block_still_accepted(self):
        ''' the block length check must not reject an untouched file
        '''
        path = self.folder / 'example.16bits.bgen'
        bfile = BgenReader(path, delay_parsing=True)
        self.assertEqual(len(bfile.samples), 500)
        self.assertEqual(bfile.samples[0], 'sample_001')
        bfile.close()
    
    def test_corrupt_variant_count(self):
        ''' an impossible variant count does not exhaust memory
        
        nvariants only bounds how much is reserved up front, so this has to fail
        by running out of file, not by trying to allocate for every variant.
        '''
        data = bytearray((self.folder / 'example.16bits.bgen').read_bytes())
        # nvariants is the 4 byte field at offset 8
        struct.pack_into('<I', data, 8, 0xFFFFFFFF)
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / 'corrupt.bgen'
            path.write_bytes(bytes(data))
            bfile = BgenReader(path, delay_parsing=True)
            with self.assertRaises(ValueError):
                bfile.rsids()
            bfile.close()
    
    def test_truncated_bgen_iteration(self):
        ''' iterating a truncated bgen raises rather than stopping early
        
        Iteration used to stop at whatever the file held, so a truncated bgen was
        indistinguishable from a complete, shorter one.
        '''
        data = (self.folder / 'example.16bits.bgen').read_bytes()
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / 'truncated.bgen'
            path.write_bytes(data[:len(data) * 90 // 100])
            
            with BgenReader(path, delay_parsing=True) as bfile:
                with self.assertRaises(ValueError):
                    for var in bfile:
                        pass
        
        # but a complete file still iterates to the end without complaint
        with BgenReader(self.folder / 'example.16bits.bgen') as bfile:
            self.assertEqual(len([x for x in bfile]), len(bfile))
    
    def test_truncated_bgen_reparse(self):
        ''' a failed parse must not leave half-parsed variants behind
        
        parse_all_variants() used to resize the variants up front, so a failure
        partway left a full sized vector of empty variants. Its 'already parsed'
        check then passed, and later calls silently returned empty rsids and
        uninitialised positions instead of raising.
        '''
        data = (self.folder / 'example.16bits.bgen').read_bytes()
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / 'truncated.bgen'
            # keep the header, but cut nearly all of the variant data
            path.write_bytes(data[:len(data) // 50])
            
            bfile = BgenReader(path, delay_parsing=True)
            # every attempt has to raise, not just the first
            for _ in range(3):
                with self.assertRaises(ValueError):
                    bfile.rsids()
            with self.assertRaises(ValueError):
                bfile.positions()
            bfile.close()
    
    def test_drop_variants_out_of_range(self):
        ''' dropping variants by an out of range index raises an error
        
        Out of range indices used to build an invalid iterator in the c++ code,
        which segfaulted rather than raising.
        '''
        path = self.folder / 'example.16bits.bgen'
        with warnings.catch_warnings():
            # drop_variants is deprecated, but still needs to not segfault
            warnings.simplefilter('ignore', DeprecationWarning)
            for idx in [-1, -100, 10000000]:
                with BgenReader(path, delay_parsing=True) as bfile:
                    with self.assertRaises(IndexError):
                        bfile.drop_variants([idx])
            
            # a bad index must not drop any of the valid ones alongside it
            with BgenReader(path, delay_parsing=True) as bfile:
                before = len(bfile)
                with self.assertRaises(IndexError):
                    bfile.drop_variants([0, 1, -1])
                self.assertEqual(len(bfile), before)
    
    def test_drop_variants_duplicate_indices(self):
        ''' dropping variants with duplicated indices raises an error
        '''
        path = self.folder / 'example.16bits.bgen'
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', DeprecationWarning)
            with BgenReader(path, delay_parsing=True) as bfile:
                with self.assertRaises(ValueError):
                    bfile.drop_variants([1, 1])
    
    def test_drop_variants(self):
        ''' dropping variants by index drops the expected number
        '''
        path = self.folder / 'example.16bits.bgen'
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', DeprecationWarning)
            with BgenReader(path, delay_parsing=True) as bfile:
                before = len(bfile)
                bfile.drop_variants([0, 2, 5])
                self.assertEqual(len(bfile), before - 3)
