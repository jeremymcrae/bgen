

import math
from pathlib import Path
import struct
import tempfile
import unittest

import numpy as np

from bgen import BgenReader, BgenWriter


def write_bgen(path, n_samples, compression=None, layout=2, bit_depth=8,
               ploidy=2, phased=False, n_alleles=2, n_variants=1):
    ''' write a small bgen to patch the genotype data of
    '''
    if phased:
        n_probs = n_alleles * ploidy
    else:
        n_probs = math.comb(ploidy + n_alleles - 1, n_alleles - 1)
    probs = np.full((n_samples, n_probs), 1 / n_probs)
    samples = [f's{i}' for i in range(n_samples)]
    alleles = [f'A{i}' for i in range(n_alleles)]
    with BgenWriter(path, n_samples, samples=samples, layout=layout,
                    compression=compression) as bfile:
        for i in range(n_variants):
            bfile.add_variant(varid=f'v{i}', rsid=f'rs{i}', chrom='01',
                              pos=i + 1, alleles=alleles, genotypes=probs,
                              ploidy=ploidy, phased=phased, bit_depth=bit_depth)
    return path


def block_offset(path, index=0):
    ''' byte offset of a variant's genotype block length field

    The length field sits after the variant header, which is varid, rsid and
    chrom as 2-byte length prefixed strings, a 4-byte position, a 2-byte allele
    count, then each allele as a 4-byte length prefixed string.
    '''
    with BgenReader(path) as bfile:
        offset = [x for x in bfile][index].fileoffset
    data = path.read_bytes()
    for _ in range(3):
        (size, ) = struct.unpack_from('<H', data, offset)
        offset += 2 + size
    offset += 4
    (n_alleles, ) = struct.unpack_from('<H', data, offset)
    offset += 2
    for _ in range(n_alleles):
        (size, ) = struct.unpack_from('<I', data, offset)
        offset += 4 + size
    return offset


def patched(path, edits):
    ''' rewrite a bgen with the given (offset, struct format, value) edits
    '''
    data = bytearray(path.read_bytes())
    for offset, fmt, value in edits:
        struct.pack_into(fmt, data, offset, value)
    path.write_bytes(bytes(data))
    return path


class TestMalformedGenotypes(unittest.TestCase):
    ''' check that malformed layout 2 genotype blocks are rejected

    The block length is stored in the bgen and sizes the decompressed buffer,
    but the fields which say how much data the variant holds (sample count,
    ploidy, allele count and bit depth) are stored separately. Nothing checked
    the two agreed, so the parsing indexed on those fields and read past the end
    of the buffer. Depending on how far past, that either returned fabricated
    probabilities or segfaulted.
    '''

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.tmpdir = Path(self.tmp.name)

    def tearDown(self):
        try:
            self.tmp.cleanup()
        except OSError:
            pass

    def test_block_without_probabilities(self):
        ''' a block holding only its preamble and ploidy is rejected

        This used to report success, with every sample given the fabricated
        probabilities that came of reading past the buffer.
        '''
        n = 500
        path = write_bgen(self.tmpdir / 'a.bgen', n)
        offset = block_offset(path)
        patched(path, [(offset, '<I', 10 + n)])

        with BgenReader(path) as bfile:
            var = next(iter(bfile))
            with self.assertRaises(ValueError):
                var.probabilities

    def test_block_with_half_the_probabilities(self):
        ''' a block cut partway through the probability data is rejected
        '''
        n = 500
        path = write_bgen(self.tmpdir / 'a.bgen', n)
        offset = block_offset(path)
        # a ploidy 2 biallelic sample stores 2 of its 3 probabilities, so at a
        # bit depth of 8 the full block needs 2 bytes per sample
        patched(path, [(offset, '<I', 10 + n + n)])

        with BgenReader(path) as bfile:
            var = next(iter(bfile))
            with self.assertRaises(ValueError):
                var.probabilities

    def test_block_truncated_inside_ploidy(self):
        ''' a block too short for even the ploidy bytes is rejected

        parse_ploidy reads a byte per sample, so this read past the buffer
        before the bit depth was validated. It happened to raise about the bit
        depth afterwards, which hid the out of bounds reads that came first.
        '''
        n = 500
        path = write_bgen(self.tmpdir / 'a.bgen', n)
        offset = block_offset(path)
        patched(path, [(offset, '<I', 20)])

        with BgenReader(path) as bfile:
            var = next(iter(bfile))
            with self.assertRaisesRegex(ValueError, 'truncated'):
                var.probabilities

    def test_ploidy_larger_than_block(self):
        ''' a ploidy which needs more data than the block holds is rejected

        min_ploidy and max_ploidy are single bytes from the file which size the
        probabilities per sample, so raising them made the reader run far past
        the end of the buffer, which segfaulted once the sample count was large
        enough.
        '''
        n = 500
        path = write_bgen(self.tmpdir / 'a.bgen', n)
        offset = block_offset(path)
        # min_ploidy and max_ploidy follow the sample and allele counts
        patched(path, [(offset + 4 + 6, '<B', 255), (offset + 4 + 7, '<B', 255)])

        with BgenReader(path) as bfile:
            var = next(iter(bfile))
            with self.assertRaises(ValueError):
                var.probabilities

    def test_bit_depth_larger_than_block(self):
        ''' a bit depth which needs more data than the block holds is rejected

        A bit depth of 32 needs four times the data of the 8 the file was
        written with, but only the 1 to 32 range was checked.
        '''
        n = 500
        path = write_bgen(self.tmpdir / 'a.bgen', n)
        offset = block_offset(path)
        # the bit depth is the last preamble byte, after the ploidy bytes
        patched(path, [(offset + 4 + 8 + n + 1, '<B', 32)])

        with BgenReader(path) as bfile:
            var = next(iter(bfile))
            with self.assertRaises(ValueError):
                var.probabilities

    def test_truncated_block_raises_every_time(self):
        ''' asking again after a rejection raises again, rather than reading

        The genotype header is only parsed once, so a caller which caught the
        error and asked for the probabilities again skipped straight past the
        check to the reads it exists to prevent.
        '''
        n = 500
        path = write_bgen(self.tmpdir / 'a.bgen', n)
        offset = block_offset(path)
        patched(path, [(offset, '<I', 10 + n)])

        with BgenReader(path) as bfile:
            var = next(iter(bfile))
            for _ in range(3):
                with self.assertRaises(ValueError):
                    var.probabilities

    def test_truncated_block_other_accessors(self):
        ''' the other genotype accessors reject a short block too

        They all parse the same genotype header, so none of them should be able
        to reach the reads.
        '''
        n = 500
        path = write_bgen(self.tmpdir / 'a.bgen', n)
        offset = block_offset(path)
        patched(path, [(offset, '<I', 10 + n)])

        with BgenReader(path) as bfile:
            var = next(iter(bfile))
            for name in ['probabilities', 'minor_allele_dosage', 'alt_dosage']:
                with self.assertRaises(ValueError):
                    getattr(var, name)

    def test_implausible_decompressed_length(self):
        ''' a decompressed length near the 32-bit limit is rejected

        The buffer was sized as the decompressed length plus 8 bytes of read
        padding in 32-bit arithmetic, so lengths this large wrapped to a tiny
        allocation and the padding was then zeroed ~4GB past the end of it.
        This segfaulted for every compression scheme.
        '''
        for compression in [None, 'zlib', 'zstd']:
            for length in [0xFFFFFFF8, 0xFFFFFFFB, 0xFFFFFFFF]:
                path = write_bgen(self.tmpdir / 'a.bgen', 20,
                                  compression=compression)
                offset = block_offset(path)
                # uncompressed blocks have no separate decompressed length, so
                # the block length is what sizes the buffer
                if compression is not None:
                    offset += 4
                patched(path, [(offset, '<I', length)])

                with BgenReader(path) as bfile:
                    var = next(iter(bfile))
                    with self.assertRaises(ValueError):
                        var.probabilities

    def test_block_too_short_for_decompressed_length(self):
        ''' a compressed block shorter than its own length field is rejected

        The 4-byte decompressed length is counted within the block length, so a
        block shorter than that underflowed the compressed length. That asked
        for a ~4GB allocation, which only failed because the read then ran out
        of file, so this checks the length itself is rejected rather than
        accepting the incidental read error the old code gave.
        '''
        for compression in ['zlib', 'zstd']:
            for length in [0, 1, 2, 3]:
                path = write_bgen(self.tmpdir / 'a.bgen', 20,
                                  compression=compression)
                offset = block_offset(path)
                patched(path, [(offset, '<I', length)])

                with BgenReader(path) as bfile:
                    var = next(iter(bfile))
                    with self.assertRaisesRegex(ValueError, 'too short'):
                        var.probabilities

    def test_valid_blocks_still_read(self):
        ''' the checks must not reject any of the blocks the writer produces

        The required size has to match the writer's layout exactly, so this
        covers the combinations where the calculation differs: each compression
        scheme, phased against unphased, and bit depths which do and do not fall
        on byte boundaries.
        '''
        for compression in [None, 'zlib', 'zstd']:
            for bit_depth in [1, 8, 15, 16, 32]:
                for n_alleles in [2, 3]:
                    for phased in [False, True]:
                        path = write_bgen(self.tmpdir / 'a.bgen', 20,
                                          compression=compression,
                                          bit_depth=bit_depth, phased=phased,
                                          n_alleles=n_alleles)
                        with BgenReader(path) as bfile:
                            var = next(iter(bfile))
                            self.assertEqual(len(var.probabilities), 20)

    def test_valid_variable_ploidy_blocks_still_read(self):
        ''' variable ploidy blocks are sized per sample, so check those too
        '''
        n = 30
        for n_alleles in [2, 3]:
            for bit_depth in [1, 8, 16, 32]:
                for phased in [False, True]:
                    path = self.tmpdir / 'a.bgen'
                    ploidy = np.tile([1, 2, 3], n // 3).astype(np.uint8)
                    if phased:
                        probs = np.full((n, n_alleles * 3), np.nan)
                        for i, ploid in enumerate(ploidy):
                            for hap in range(ploid):
                                cols = slice(hap * n_alleles,
                                             (hap + 1) * n_alleles)
                                probs[i, cols] = 1 / n_alleles
                    else:
                        widest = math.comb(3 + n_alleles - 1, n_alleles - 1)
                        probs = np.full((n, widest), np.nan)
                        for i, ploid in enumerate(ploidy):
                            k = math.comb(int(ploid) + n_alleles - 1,
                                          n_alleles - 1)
                            probs[i, :k] = 1 / k
                    with BgenWriter(path, n, samples=[f's{i}' for i in range(n)],
                                    layout=2, compression='zstd') as bfile:
                        bfile.add_variant(
                            varid='v', rsid='rs', chrom='01', pos=1,
                            alleles=[f'A{i}' for i in range(n_alleles)],
                            genotypes=probs, ploidy=ploidy, phased=phased,
                            bit_depth=bit_depth)
                    with BgenReader(path) as bfile:
                        var = next(iter(bfile))
                        self.assertEqual(len(var.probabilities), n)

    def test_malformed_variant_does_not_stop_the_others(self):
        ''' a rejected variant leaves the rest of the file readable

        The ploidy is patched rather than the block length here, since the block
        length is what locates the next variant, and shortening it would make
        the rest of the file genuinely unreadable.
        '''
        n = 100
        path = write_bgen(self.tmpdir / 'a.bgen', n, n_variants=3)
        offset = block_offset(path, index=1)
        patched(path, [(offset + 4 + 6, '<B', 255), (offset + 4 + 7, '<B', 255)])

        with BgenReader(path) as bfile:
            variants = [x for x in bfile]
            self.assertEqual(len(variants[0].probabilities), n)
            with self.assertRaises(ValueError):
                variants[1].probabilities
            self.assertEqual(len(variants[2].probabilities), n)
