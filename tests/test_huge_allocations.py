''' test that a small bgen cannot drive a huge allocation

Several sizes in a bgen come from 32-bit fields, and a corrupt or hostile file can set
them far beyond what the file itself holds. Where those fields size an allocation before
anything checks them, a file of a hundred bytes can ask for gigabytes of memory:

  * an allele's length prefix. read_prefixed_string used to resize a std::string to the
    claimed length before reading any bytes, and resize() zero fills, so the memory was
    genuinely committed rather than merely reserved. A 105 byte bgen claiming a 3 GB
    allele grew the process by 3071 MB, an amplification of about 30 million to one.
  * the sample count, when the bgen carries no sample IDs. Building the placeholder IDs
    is deferred so that the count cannot size an allocation at open time, but a caller
    can ask for the IDs without ever reading a variant, and layout 2 does not repeat the
    sample count per variant, so nothing else validates it.

Both are checked here by running a child with a capped address space. Without the fix the
child dies with MemoryError or std::bad_alloc; with it, the file is rejected as malformed.
The cap is what makes these safe to run: the point is a claim too large to serve, so
serving it would page out the machine.

Allele lengths cannot simply be capped, since indels and structural variants carry real
sequence, so a legitimately long allele is checked here too.
'''

import os
from pathlib import Path
import resource
import struct
import subprocess
import sys
import tempfile
import unittest

import numpy as np

from bgen import BgenReader, BgenWriter

# enough address space for the interpreter and numpy, but far less than the claims below
CHILD_MEMORY_CAP = 1200 * 1024 * 1024

# read a bgen in a child with a capped address space, reporting what came back. The cap is
# set in the child rather than the parent so the test process keeps its own memory.
CHILD = r'''
import resource, sys
cap = int(sys.argv[3])
resource.setrlimit(resource.RLIMIT_AS, (cap, cap))
from bgen import BgenReader
path, action = sys.argv[1], sys.argv[2]
try:
    bfile = BgenReader(path, delay_parsing=True)
    if action == 'samples':
        print(f'SAMPLES {len(bfile.samples)}')
    else:
        for var in bfile:
            var.probabilities
        print('PARSED')
except MemoryError as err:
    print(f'MEMORY {err}')
except Exception as err:
    print(f'{type(err).__name__} {err}')
'''


def run_capped(path, action):
    ''' read a bgen in a child whose address space is capped

    A pre-check failure shows up as the reader's own error, whereas an allocation sized
    straight from the file shows up as MemoryError, bad_alloc or a fatal signal.
    '''
    env = dict(os.environ, PYTHONPATH=os.pathsep.join(p for p in sys.path if p))
    done = subprocess.run([sys.executable, '-c', CHILD, str(path), action,
                           str(CHILD_MEMORY_CAP)], capture_output=True, text=True,
                          env=env, timeout=600)
    out = done.stdout.strip()
    if done.returncode != 0 and not out:
        # died without reporting, e.g. killed by a signal
        return f'DIED rc={done.returncode} {done.stderr.strip()[-200:]}'
    return out


class TestHugeAllocations(unittest.TestCase):
    ''' check overstated lengths are rejected rather than allocated '''

    def setUp(self):
        self.tmpdir = Path(tempfile.mkdtemp())

    def _write(self, path, n_samples=4, alleles=('A', 'C'), sample_ids=True):
        ''' write a small valid bgen to corrupt a single field of '''
        n_probs = 3
        probs = np.full((n_samples, n_probs), 1 / n_probs)
        samples = [f's{i}' for i in range(n_samples)] if sample_ids else None
        with BgenWriter(path, n_samples, samples=samples, compression=None) as bfile:
            bfile.add_variant(varid='v', rsid='rs', chrom='01', pos=1,
                              alleles=list(alleles), genotypes=probs)
        return path

    def _allele_length_offset(self, path):
        ''' byte offset of the first allele's 4 byte length prefix

        The variant header is varid, rsid and chrom as 2-byte length prefixed strings,
        then a 4-byte position and a 2-byte allele count.
        '''
        data = path.read_bytes()
        with BgenReader(path) as bfile:
            offset = next(iter(bfile)).fileoffset
        for _ in range(3):
            (size, ) = struct.unpack_from('<H', data, offset)
            offset += 2 + size
        return offset + 4 + 2

    def test_overstated_allele_length_is_not_allocated(self):
        ''' an allele claiming more bytes than the file holds is rejected

        The length is a 32-bit field, so the claim can be 4 GB from a file of a hundred
        bytes. Reading in bounded chunks means only what the file can supply is ever
        allocated.
        '''
        path = self._write(self.tmpdir / 'allele.bgen')
        offset = self._allele_length_offset(path)
        for claim in [0xFFFFFFF0, 3 * 1024 * 1024 * 1024, 700 * 1024 * 1024]:
            with self.subTest(claim=claim):
                bad = self.tmpdir / f'allele{claim}.bgen'
                bad.write_bytes(path.read_bytes())
                data = bytearray(bad.read_bytes())
                struct.pack_into('<I', data, offset, claim)
                bad.write_bytes(bytes(data))

                self.assertLess(bad.stat().st_size, 1000)
                outcome = run_capped(bad, 'iterate')
                self.assertNotIn('MEMORY', outcome)
                self.assertNotIn('DIED', outcome)
                self.assertNotIn('bad_alloc', outcome)
                # the allele cannot be read, so the file is short of the variants it
                # claims, which is what the reader reports
                self.assertIn('ValueError', outcome)

    def test_long_alleles_still_read(self):
        ''' a genuinely long allele must survive, so the fix cannot be a cap

        Indels and structural variants store sequence as the allele, so multi-megabyte
        alleles are legitimate and have to round trip intact.
        '''
        for length in [1000, 100000, 2000000]:
            with self.subTest(length=length):
                alleles = ('A', 'C' * length)
                path = self._write(self.tmpdir / f'long{length}.bgen', alleles=alleles)
                with BgenReader(path) as bfile:
                    self.assertEqual(next(iter(bfile)).alleles, list(alleles))

    def _drop_sample_ids(self, path):
        ''' clear the has_sample_ids flag, which is bit 31 of the header's flags word '''
        data = bytearray(path.read_bytes())
        (header_length, ) = struct.unpack_from('<I', data, 4)
        flags_offset = 4 + header_length - 4
        (flags, ) = struct.unpack_from('<I', data, flags_offset)
        struct.pack_into('<I', data, flags_offset, flags & ~(1 << 31))
        path.write_bytes(bytes(data))
        return path

    def test_placeholder_sample_ids_are_bounded(self):
        ''' a bgen without sample IDs cannot claim more samples than it could hold

        With no IDs in the file the reader generates placeholders, and it defers that so
        the count cannot allocate at open time. But asking for the samples reaches it
        without reading a variant, and layout 2 stores no per variant sample count, so
        the deferral alone left the count unchecked.
        '''
        path = self._drop_sample_ids(self._write(self.tmpdir / 'noids.bgen'))
        for claim in [500_000_000, 2_000_000_000, 0xFFFFFFF0]:
            with self.subTest(claim=claim):
                bad = self.tmpdir / f'noids{claim}.bgen'
                data = bytearray(path.read_bytes())
                struct.pack_into('<I', data, 12, claim)
                bad.write_bytes(bytes(data))

                self.assertLess(bad.stat().st_size, 1000)
                outcome = run_capped(bad, 'samples')
                self.assertNotIn('MEMORY', outcome)
                self.assertNotIn('DIED', outcome)
                self.assertIn('ValueError', outcome)
                self.assertIn('cannot describe', outcome)

    def test_placeholder_sample_ids_still_generated(self):
        ''' the real count must still produce placeholder IDs

        The bound has to reject only counts the file cannot back, or ordinary bgens
        without sample IDs would stop working.
        '''
        path = self._drop_sample_ids(self._write(self.tmpdir / 'plain.bgen',
                                                 n_samples=4))
        with BgenReader(path) as bfile:
            self.assertEqual(bfile.samples, ['0', '1', '2', '3'])


if __name__ == '__main__':
    unittest.main()
