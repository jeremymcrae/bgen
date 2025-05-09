
from pathlib import Path
import os
import unittest
import sys

from bgen import BgenReader

from tests.utils import load_gen_data, arrays_equal

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
        