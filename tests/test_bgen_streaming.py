
from pathlib import Path
import os
import struct
import subprocess
import tempfile
import unittest
import sys

import numpy as np

from bgen import BgenReader, BgenWriter

from tests.utils import load_gen_data, arrays_equal

try:
    import resource
    HAS_RLIMIT = hasattr(resource, 'RLIMIT_AS')
except ImportError:
    # windows has no resource module, so the memory capped test is skipped there
    HAS_RLIMIT = False

# how long to give a child before treating it as stuck. These children read a bgen from
# a pipe and exit, which takes well under a second, so this only trips on a real hang.
# Without it a child that never exits leaves the whole test run waiting forever.
CHILD_TIMEOUT = 120

# a child that caps its own address space has to import before capping. OpenBLAS sizes
# its thread arenas by core count, and if RLIMIT_AS is too tight to mmap them its init
# retries forever rather than failing, so capping first hangs on machines with enough
# cores. Capping after the imports leaves the allocation under test just as constrained.
CAP_AFTER_IMPORT = '''
import resource
# measure what the imports already took, so the headroom below does not depend on how
# much address space this machine's numpy needed. Only linux publishes this, and
# elsewhere the cap just falls back to being an absolute one, as it used to be
used = 0
try:
    with open('/proc/self/status') as handle:
        for line in handle:
            if line.startswith('VmSize'):
                used = int(line.split()[1]) * 1024
except OSError:
    pass
resource.setrlimit(resource.RLIMIT_AS, (used + 1024 ** 3, used + 1024 ** 3))
'''

def run_piped(code, data):
    ''' run code in a subprocess with data fed to its stdin

    The bgen is piped in so the reader sees a stream it cannot seek, which is the point
    of these tests. It runs in a subprocess so the parent's own stdin is left alone,
    and so a crash shows up as a signal rather than taking the test run down.

    The child gets this process' sys.path, so it imports the build under test rather
    than any installed copy.
    '''
    env = dict(os.environ, PYTHONPATH=os.pathsep.join(p for p in sys.path if p))
    try:
        return subprocess.run([sys.executable, '-c', code], input=data,
                              capture_output=True, env=env, timeout=CHILD_TIMEOUT)
    except subprocess.TimeoutExpired as e:
        # say what was being run, since a bare timeout gives no clue which child stuck
        out = (e.stdout or b'').decode('utf8', 'replace')
        err = (e.stderr or b'').decode('utf8', 'replace')
        raise AssertionError(
            f'child did not exit within {CHILD_TIMEOUT}s, so it is stuck rather than '
            f'slow.\ncode:\n{code}\nstdout so far:\n{out}\nstderr so far:\n{err}')

class TestBgenStream(unittest.TestCase):
    ''' class to make sure BgenReader works correctly
    '''
    
    @classmethod
    def setUpClass(cls):
        cls.gen_data = load_gen_data()
    
    def setUp(self):
        ''' set path to folder with test data
        '''
        # can't write > 2 ** 16 bytes to linux write_fd pipe before something stalls
        self.max_buff = 65536
        self.folder = Path(__file__).parent /  "data"
    
    @unittest.skipIf(sys.platform == "win32", "windows lacks /dev/stdin")
    def test_bgen_streaming_parse_all_variants(self):
        ''' check we can parse all variants of a bgen from stdin
        
        Anything which needs every variant at once (rsids, varids, chroms,
        positions, indexing, drop_variants) parses them all up front, which
        collects the variants into a list. A streamed bgen loads its genotypes
        as each variant is constructed, so this used to hand the same genotype
        buffers to two variants, and freeing them twice aborted the process.
        '''
        path = self.folder / 'example.16bits.zstd.bgen'
        code = ('from bgen import BgenReader\n'
                'b = BgenReader("/dev/stdin")\n'
                'assert len(b.rsids()) == 199\n'
                'assert len(b.positions()) == 199\n'
                'print("ok")\n')
        proc = run_piped(code, path.read_bytes())
        # a double free aborts, which shows up as a negative (signal) returncode
        self.assertEqual(proc.returncode, 0, msg=proc.stderr.decode('utf8', 'replace'))
        self.assertEqual(proc.stdout.decode('utf8').strip(), 'ok')
    
    @unittest.skipIf(sys.platform == "win32", "windows lacks /dev/stdin")
    @unittest.skipUnless(HAS_RLIMIT, 'needs RLIMIT_AS to cap memory')
    def test_streamed_bgen_corrupt_sample_count(self):
        ''' a corrupt sample count on stdin must not size an allocation
        
        The file size check in the reader needs to seek to the end, so it is skipped
        for a stream. That left the sample block, whose length and count both come
        from the file, as the only thing sizing the sample list, so a stream that
        claims millions of samples used to allocate for them all before running out
        of data.
        
        The cap is applied after the imports, since capping first hangs OpenBLAS's
        thread init.
        '''
        n = 50_000_000
        # a header claiming n samples, and a sample block sized to match, but with
        # no ID data at all behind it
        block = struct.pack('<II', 8 + n * 2, n)
        header_length = 20
        flags = (2 << 2) | (1 << 31)
        data = struct.pack('<IIII', header_length + len(block), header_length, 0, n)
        data += b'bgen' + struct.pack('<I', flags) + block
        code = ('from bgen import BgenReader\n'
                + CAP_AFTER_IMPORT +
                'try:\n'
                '    BgenReader("/dev/stdin")\n'
                'except ValueError:\n'
                '    print("ok")\n')
        proc = run_piped(code, data)
        self.assertEqual(proc.returncode, 0,
                         msg=proc.stderr.decode('utf8', 'replace'))
        self.assertEqual(proc.stdout.decode('utf8').strip(), 'ok')
    
    @unittest.skipIf(sys.platform == "win32", "windows lacks /dev/stdin")
    @unittest.skipUnless(HAS_RLIMIT, 'needs RLIMIT_AS to cap memory')
    def test_streamed_bgen_corrupt_variant_count(self):
        ''' a corrupt variant count on stdin must not size an allocation
        
        The reserve for the variant list is bounded by what the file could hold, but a
        stream cannot be seeked, so its size is unknown and that bound was skipped. A
        Variant is a few hundred bytes, so a stream claiming millions of them used to
        ask for several GB up front and fail with a memory error, rather than reporting
        the file as truncated.
        
        The cap is applied after the imports, since capping first hangs OpenBLAS's
        thread init.
        '''
        # a header with no variant data behind it, so the parse should stop immediately
        for nvariants in [10_000_000, 4_294_967_295]:
            with self.subTest(nvariants=nvariants):
                header_length = 20
                flags = 2 << 2  # layout 2, and no sample IDs in the bgen
                data = struct.pack('<IIII', header_length, header_length,
                                   nvariants, 1)
                data += b'bgen' + struct.pack('<I', flags)
                code = ('from bgen import BgenReader\n'
                        + CAP_AFTER_IMPORT +
                        'try:\n'
                        # the reserve happens when the variants are parsed, which is
                        # what asking for the rsids does
                        '    BgenReader("/dev/stdin").rsids()\n'
                        'except ValueError:\n'
                        '    print("ok")\n')
                proc = run_piped(code, data)
                self.assertEqual(proc.returncode, 0,
                                 msg=proc.stderr.decode('utf8', 'replace'))
                self.assertEqual(proc.stdout.decode('utf8').strip(), 'ok')
    
    @unittest.skipIf(sys.platform == "win32", "windows lacks /dev/stdin")
    def test_streamed_bgen_with_many_variants(self):
        ''' a streamed bgen with more variants than the reserve cap still reads
        
        The reserve on a stream is capped, so a bgen holding more variants than the cap
        grows the list as they arrive. This makes sure that growth is not mistaken for
        corruption, and does not drop or reorder anything, by checking every variant is
        present and in order.
        '''
        # more variants than MAX_VARIANT_RESERVE in reader.cpp, so the list has to grow
        nvariants = (1 << 16) + 5
        genotypes = np.array([[0.7, 0.2, 0.1]] * 2)
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / 'many_variants.bgen'
            with BgenWriter(path, n_samples=2, samples=['a', 'b']) as bfile:
                for idx in range(nvariants):
                    bfile.add_variant(varid=f'v{idx}', rsid=f'rs{idx}', chrom='1',
                                      pos=idx + 1, alleles=['A', 'C'],
                                      genotypes=genotypes, bit_depth=8)
            code = ('from bgen import BgenReader\n'
                    'b = BgenReader("/dev/stdin")\n'
                    'ids = b.rsids()\n'
                    f'assert len(ids) == {nvariants}, len(ids)\n'
                    f'assert ids == [f"rs{{i}}" for i in range({nvariants})], "ids differ"\n'
                    'print("ok")\n')
            proc = run_piped(code, path.read_bytes())
        self.assertEqual(proc.returncode, 0,
                         msg=proc.stderr.decode('utf8', 'replace'))
        self.assertEqual(proc.stdout.decode('utf8').strip(), 'ok')
    
    @unittest.skipIf(sys.platform == "win32", "windows lacks /dev/stdin")
    def test_streamed_bgen_sample_block_length_mismatch(self):
        ''' a bad sample block length is caught on a stream as well as a file
        
        The check counts the bytes the IDs occupy rather than asking the stream for its
        position, since tellg() returns -1 on a pipe. Reading the same bytes from a
        pipe must therefore reach the same verdict as reading them from a file.
        '''
        path = self.folder / 'example.16bits.bgen'
        for delta in [10, -10]:
            with self.subTest(delta=delta):
                data = bytearray(path.read_bytes())
                header_length = struct.unpack_from('<I', data, 4)[0]
                at = 4 + header_length
                block_length = struct.unpack_from('<I', data, at)[0]
                struct.pack_into('<I', data, at, block_length + delta)
                code = ('from bgen import BgenReader\n'
                        'try:\n'
                        '    BgenReader("/dev/stdin", delay_parsing=True)\n'
                        '    print("accepted")\n'
                        'except ValueError:\n'
                        '    print("rejected")\n')
                proc = run_piped(code, bytes(data))
                self.assertEqual(proc.returncode, 0,
                                 msg=proc.stderr.decode('utf8', 'replace'))
                self.assertEqual(proc.stdout.decode('utf8').strip(), 'rejected')
    
    @unittest.skipIf(sys.platform == "win32", "windows lacks /dev/stdin")
    def test_streamed_bgen_with_many_samples(self):
        ''' a streamed bgen with more samples than the reserve cap still reads
        
        A stream cannot be seeked, so its size is unknown and the sample count cannot
        be checked against it. The IDs are therefore read into a list that grows past
        a fixed reserve, and this makes sure that growth is not mistaken for
        corruption.
        '''
        n = (1 << 20) + 5
        ids = b''.join(struct.pack('<H', len(s)) + s
                       for s in (b'sample_%d' % i for i in range(n)))
        block = struct.pack('<II', 8 + len(ids), n) + ids
        header_length = 20
        flags = (2 << 2) | (1 << 31)
        data = struct.pack('<IIII', header_length + len(block), header_length, 0, n)
        data += b'bgen' + struct.pack('<I', flags) + block
        code = ('from bgen import BgenReader\n'
                'b = BgenReader("/dev/stdin")\n'
                's = b.samples\n'
                f'assert len(s) == {n}, len(s)\n'
                'assert s[0] == "sample_0"\n'
                f'assert s[-1] == "sample_{n - 1}"\n'
                'print("ok")\n')
        proc = run_piped(code, data)
        self.assertEqual(proc.returncode, 0,
                         msg=proc.stderr.decode('utf8', 'replace'))
        self.assertEqual(proc.stdout.decode('utf8').strip(), 'ok')
    
    @unittest.skipIf(sys.platform == "win32", "haven't figured out file handle " \
                                              "duplication and writing on windows " \
                                              "at the required buffer size")
    def test_bgen_streaming(self):
        ''' check that was can open a bgen file from stdin
        '''
        # Save original stdin and stdout file descriptors
        original_stdin_fd = os.dup(sys.stdin.fileno())
        original_stdout_fd = os.dup(sys.stdout.fileno())
        read_fd, write_fd = os.pipe()
        
        path = self.folder / 'example.16bits.zstd.bgen'
        with open(path, 'rb') as handle:
            data = handle.read()
        
        try:
            os.write(write_fd, data[:self.max_buff])
            os.close(write_fd) # Close the write end once done writing
            
            # Redirect stdin file descriptor to the read end of the pipe
            os.dup2(read_fd, sys.stdin.fileno())
            os.close(read_fd) # Close the original read end after duplicating
            
            with BgenReader(sys.stdin) as bfile:
                header = bfile.header
                self.assertEqual(header.nvariants, 199)
                self.assertEqual(header.nsamples, 500)
                
                var1 = next(bfile)
                var1.alt_dosage # load dosage for the first variant
                g = self.gen_data[0]
                self.assertTrue(arrays_equal(g.probabilities, var1.probabilities, bit_depth=16))
                
                var2 = next(bfile)  # don't load any genotype probabilities for this variant
        finally:
            # Restore original stdin and stdout file descriptors
            os.dup2(original_stdin_fd, sys.stdin.fileno())
            os.close(original_stdin_fd)
            sys.stdout = sys.__stdout__ # Use sys.__stdout__ to get the true original stdout
            # Close the duplicated original file descriptors
            os.close(original_stdout_fd)
        
        # can't get header for closed bgen file
        with self.assertRaises(ValueError):
            bfile.header
        self.assertEqual(var1.pos, 2000)
        self.assertEqual(var2.pos, 3000)
        
        # we should still be able to access genotype probabilities for the first variant, as
        # these were loaded while the file was still open
        g = self.gen_data[0]
        self.assertTrue(arrays_equal(g.probabilities, var1.probabilities, bit_depth=16))
        
        # # in contrast, we shouldn't be able to load the genotype probabilities for the 
        # # second variant, as these were not loaded while the bgen was open. It should raise
        # # and error, but currently does not. I haven't figured out why not yet.
        # with self.assertRaises(ValueError):
        #     var2.probabilities
        