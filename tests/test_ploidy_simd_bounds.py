''' test that the vectorised ploidy helpers stay inside their arrays

fast_ploidy_sum and fast_range load a block of bytes per iteration, so their loop guards
have to account for the width of the load, not just the step. The SSE4 sum stepped four
bytes at a time but loaded sixteen, and its guard allowed a final iteration whose load ran
up to three bytes past the end of the ploidy array. Only the low four bytes are converted,
so the values were always right, but the read itself is out of bounds and faults if the
array happens to end within three bytes of a mapping boundary.

None of this is reachable from Python on a machine with AVX2, since fast_ploidy_sum only
calls the SSE4 version when AVX2 is missing, and there is no way to ask for the SSE4 path
at runtime. So these tests compile the helpers directly and call them with the ploidy
array butted up against an unreadable page, which turns any overread into a crash rather
than something that depends on allocator slack.

The same helpers also used to widen the bytes to 32 bits and accumulate there, which
overflowed partway through a large cohort, so that is checked here too.
'''

import os
from pathlib import Path
import shutil
import subprocess
import sys
import sysconfig
import tempfile
import unittest

SRC = Path(__file__).parent.parent / 'src'

# drive the helpers on an array that ends exactly at an unreadable page, so a load which
# runs past the end faults instead of quietly reading whatever follows
PROBE = r'''
#include <cstdint>
#include <cstdio>
#include <unistd.h>
#include <sys/mman.h>

#include "utils.cpp"

using namespace bgen;

int main() {
  long page = sysconf(_SC_PAGESIZE);
  int failures = 0;
  // sizes that are not multiples of four are what push the final iteration over the end
  for (std::uint32_t size = 1; size <= 300; size++) {
    std::uint8_t * region = (std::uint8_t *) mmap(nullptr, page * 2,
                                                  PROT_READ | PROT_WRITE,
                                                  MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
    if (region == MAP_FAILED) { return 2; }
    if (mprotect(region + page, page, PROT_NONE) != 0) { return 2; }

    std::uint8_t * x = region + page - size;
    std::uint64_t want = 0;
    std::uint8_t want_min = 255, want_max = 0;
    for (std::uint32_t j = 0; j < size; j++) {
      x[j] = (std::uint8_t) ((j * 31 + 5) % 256);
      want += x[j];
      if (x[j] < want_min) { want_min = x[j]; }
      if (x[j] > want_max) { want_max = x[j]; }
    }

    // the SSE4 helpers are static and skipped on this cpu, so call them directly, then
    // finish the remainder the way fast_ploidy_sum does
    std::uint32_t i = 0;
    std::uint32_t sz = size;
    std::uint64_t got = ploidy_sum_sse4(x, sz, i);
    for (; i < size; i++) { got += x[i]; }
    if (got != want) { failures++; }

    std::uint8_t got_min = 255, got_max = 0;
    size_t k = 0;
    sz = size;
    range_sse4(x, sz, k, got_min, got_max);
    for (; k < size; k++) {
      if (x[k] < got_min) { got_min = x[k]; }
      if (x[k] > got_max) { got_max = x[k]; }
    }
    if (got_min != want_min || got_max != want_max) { failures++; }

    // and the same through the public entry points, which take the avx2 path here
    sz = size;
    if (fast_ploidy_sum(x, sz) != want) { failures++; }
    sz = size;
    Range r = fast_range(x, sz);
    if (r._min != want_min || r._max != want_max) { failures++; }

    munmap(region, page * 2);
  }
  printf("%d\n", failures);
  return 0;
}
'''


# sum a large array of the largest possible ploidy, which is where a 32-bit accumulator
# wraps. This allocates one byte per sample, so it stays affordable
OVERFLOW_PROBE = r'''
#include <cstdint>
#include <cstdio>
#include <new>
#include <vector>

#include "utils.cpp"

using namespace bgen;

int main() {
  int failures = 0;
  // max_ploidy is read as a whole byte from the file, and the constant ploidy branch of
  // parse_ploidy memsets the array to it without masking, so 255 is reachable
  const std::uint8_t worst = 255;
  // well past where 32-bit lanes wrapped, which was about 33.7M for avx2 and 16.8M for sse4
  for (std::uint32_t size : {16843024u, 33686048u, 40000000u}) {
    std::vector<std::uint8_t> values;
    try {
      values.assign(size, worst);
    } catch (const std::bad_alloc &) {
      continue;
    }
    std::uint64_t want = (std::uint64_t) size * worst;
    std::uint32_t n = size;
    if (fast_ploidy_sum(values.data(), n) != want) { failures++; }

    // the sse4 helper is skipped on a cpu with avx2, so drive it directly as well
    std::uint32_t i = 0;
    n = size;
    std::uint64_t got = ploidy_sum_sse4(values.data(), n, i);
    for (; i < size; i++) { got += values[i]; }
    if (got != want) { failures++; }
  }
  printf("%d\n", failures);
  return 0;
}
'''


def compiler():
    ''' the compiler python itself was built with, so the probe matches the extension
    '''
    cxx = sysconfig.get_config_var('CXX') or os.environ.get('CXX', 'c++')
    return cxx.split()


@unittest.skipUnless(sys.platform.startswith('linux'),
                     'needs mmap and mprotect to place a guard page')
@unittest.skipUnless(shutil.which(compiler()[0]) is not None,
                     'needs a c++ compiler to build the probe')
class TestPloidySimdBounds(unittest.TestCase):
    ''' check the vectorised ploidy helpers do not read past their arrays
    '''

    @classmethod
    def setUpClass(cls):
        ''' build the probes once, since compiling is the slow part
        '''
        cls.tmpdir = tempfile.mkdtemp()
        cls.binaries = {}
        cls.build_errors = {}
        for name, code in [('bounds', PROBE), ('overflow', OVERFLOW_PROBE)]:
            source = Path(cls.tmpdir) / f'{name}.cpp'
            source.write_text(code)
            binary = Path(cls.tmpdir) / name
            build = subprocess.run(compiler() + ['-std=c++11', '-O2', f'-I{SRC}',
                                                 '-o', str(binary), str(source)],
                                   capture_output=True, text=True, timeout=300)
            cls.binaries[name] = binary
            cls.build_errors[name] = (build.stderr if build.returncode != 0 else None)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmpdir, ignore_errors=True)

    def run_probe(self, name):
        ''' run one of the compiled probes and return its output
        '''
        if self.build_errors[name] is not None:
            self.skipTest(f'could not build the {name} probe: '
                          f'{self.build_errors[name][-500:]}')
        proc = subprocess.run([str(self.binaries[name])], capture_output=True,
                              text=True, timeout=300)
        # a read past the end of the array lands on the guard page, which is a signal
        self.assertGreaterEqual(proc.returncode, 0,
                                msg=f'crashed with signal {-proc.returncode}, so a '
                                    'helper read past the end of the ploidy array')
        self.assertEqual(proc.returncode, 0, msg=proc.stderr[-500:])
        return proc.stdout.strip()

    def test_ploidy_helpers_stay_in_bounds(self):
        ''' the helpers must not read past a ploidy array that ends at a guard page
        '''
        self.assertEqual(self.run_probe('bounds'), '0',
                         msg='a helper returned the wrong sum or range')

    def test_ploidy_sum_does_not_overflow(self):
        ''' the sums must stay correct for a big cohort at the largest ploidy

        The bytes used to be widened to 32 bits and accumulated there, so both an
        individual lane and the sum of the lanes wrapped partway through a large cohort.
        The size at which that happened depended on the ploidy values, so it could not be
        guarded by a fixed sample count. They are accumulated into 64-bit counters now,
        which cannot overflow for any array a 32-bit size can describe.
        '''
        self.assertEqual(self.run_probe('overflow'), '0',
                         msg='a ploidy sum overflowed')


if __name__ == '__main__':
    unittest.main()
