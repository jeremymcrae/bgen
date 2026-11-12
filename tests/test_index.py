
from pathlib import Path
import tempfile
import unittest
from unittest import mock

import numpy as np


from bgen import BgenReader, BgenWriter
from bgen.index import Index

from tests.utils import load_gen_data

class TestBgenIndex(unittest.TestCase):
    ''' class to make sure bgen.index.Index works correctly
    '''
    @classmethod
    def setUpClass(cls):
        cls.gen_data = load_gen_data()
    
    def setUp(self):
        ''' set path to folder with test data
        '''
        self.folder = Path(__file__).parent /  "data"
    
    def test_index_opens(self):
        ''' loads index when available
        '''
        with BgenReader(self.folder / 'example.15bits.bgen') as bfile:
            self.assertFalse(bfile._check_for_index(str(self.folder / 'example.15bits.bgen')))
        
        with BgenReader(self.folder / 'example.16bits.bgen') as bfile:
            self.assertTrue(bfile._check_for_index(str(self.folder / 'example.16bits.bgen')))
    
    def test_rechecking_for_an_index_closes_the_previous_one(self):
        ''' an index replaced by a second check is closed, rather than left open '''
        opened = []
        class Recorder(Index):
            def __init__(self, path):
                super().__init__(path)
                opened.append(self)
        
        path = self.folder / 'example.16bits.bgen'
        with mock.patch('bgen.reader.Index', Recorder):
            with BgenReader(path) as bfile:
                bfile._check_for_index(str(path))
                self.assertEqual(len(opened), 2)
                self.assertIsNone(opened[0].conn)
                self.assertIsNotNone(opened[1].conn)
    
    def test_index_fetch(self):
        ''' fetches file offsets
        '''
        chrom = '01'
        start = 5000
        stop = 50000
        
        index = Index(self.folder / 'example.16bits.bgen.bgi')
        self.addCleanup(index.close)
        self.assertTrue(len(list(index.fetch(chrom))) == len(self.gen_data))
        self.assertTrue(len(list(index.fetch('02'))) == 0)
        self.assertTrue(len(list(index.fetch(chrom, start * 100, stop * 100))) == 0)
        
        # check for a whole chromosome
        chrom_offsets = list(index.fetch(chrom))
        self.assertTrue(len(chrom_offsets) > 0)
        self.assertTrue(len(chrom_offsets) == len(self.gen_data))
        
        # check for all variants following a position
        after_pos_offsets = list(index.fetch(chrom, start))
        self.assertTrue(len(after_pos_offsets) > 0)
        self.assertTrue(len(after_pos_offsets) == len([x for x in self.gen_data if start <= x.pos]))
        
        # check for all variants within a region
        in_region_offsets = list(index.fetch(chrom, start, stop))
        self.assertTrue(len(in_region_offsets) > 0)
        self.assertTrue(len(in_region_offsets) == len([x for x in self.gen_data if start <= x.pos <= stop]))
        
        # make sure the queries return different lists
        self.assertTrue(len(chrom_offsets) != len(after_pos_offsets))
        self.assertTrue(len(chrom_offsets) != len(in_region_offsets))
        self.assertTrue(len(after_pos_offsets) != len(in_region_offsets))
    
    def test_index_fetch_stop_without_start(self):
        ''' fetches file offsets for variants up to a position
        
        A stop without a start is the one combination the test above misses.
        '''
        chrom = '01'
        stop = 50000
        
        index = Index(self.folder / 'example.16bits.bgen.bgi')
        self.addCleanup(index.close)
        before_pos_offsets = list(index.fetch(chrom, stop=stop))
        expected = [x for x in self.gen_data if x.pos <= stop]
        self.assertTrue(len(before_pos_offsets) > 0)
        self.assertEqual(len(before_pos_offsets), len(expected))
        
        # the offsets have to be the ones the full scan gives for those variants
        in_region_offsets = list(index.fetch(chrom, 0, stop))
        self.assertEqual(sorted(before_pos_offsets), sorted(in_region_offsets))
        
        # a stop below every variant gives nothing, and one above gives them all
        positions = [x.pos for x in self.gen_data]
        self.assertEqual(list(index.fetch(chrom, stop=min(positions) - 1)), [])
        self.assertEqual(len(list(index.fetch(chrom, stop=max(positions)))),
                         len(self.gen_data))
        
        # and it still has to respect the chromosome
        self.assertEqual(list(index.fetch('02', stop=stop)), [])
    
    def test_reader_fetch_stop_without_start(self):
        ''' BgenReader.fetch yields the variants up to a position
        '''
        chrom = '01'
        stop = 50000
        
        with BgenReader(self.folder / 'example.16bits.bgen') as bfile:
            variants = list(bfile.fetch(chrom, stop=stop))
            expected = [x for x in self.gen_data if x.pos <= stop]
            self.assertEqual(len(variants), len(expected))
            self.assertTrue(all(x.pos <= stop for x in variants))
            self.assertEqual(sorted(x.rsid for x in variants),
                             sorted(x.rsid for x in expected))
    
    def test_index_queries_do_not_disturb_a_live_fetch(self):
        ''' another query must not cut short a fetch which is part way through
        
        fetch is a generator, so it hands out rows as it steps through them, and
        needs a cursor of its own.
        '''
        index = Index(self.folder / 'example.16bits.bgen.bgi')
        self.addCleanup(index.close)
        expected = list(index.fetch('01'))
        self.assertGreater(len(expected), 1)
        
        interleave = {
            'offset_by_rsid': lambda: index.offset_by_rsid('RSID_5'),
            'offset_by_pos': lambda: index.offset_by_pos(5000),
            'offset_by_index': lambda: index.offset_by_index(0),
            'rsids': lambda: index.rsids,
            'chroms': lambda: index.chroms,
            'positions': lambda: index.positions,
            'len': lambda: len(index),
            'nested fetch': lambda: list(index.fetch('01', 5000, 6000)),
        }
        for label, query in interleave.items():
            with self.subTest(query=label):
                got = []
                for offset in index.fetch('01'):
                    got.append(offset)
                    query()
                self.assertEqual(got, expected)
    
    def test_index_supports_concurrent_fetches(self):
        ''' several fetches can be part way through at once '''
        index = Index(self.folder / 'example.16bits.bgen.bgi')
        self.addCleanup(index.close)
        expected = list(index.fetch('01'))
        
        generators = [index.fetch('01') for _ in range(4)]
        # start them all before draining any, so they overlap
        firsts = [next(x) for x in generators]
        self.assertEqual(firsts, [expected[0]] * len(generators))
        for gen in generators:
            self.assertEqual([expected[0]] + list(gen), expected)
    
    def test_reader_fetch_survives_other_reader_calls(self):
        ''' the same holds through the public BgenReader interface '''
        path = self.folder / 'example.16bits.bgen'
        with BgenReader(path) as bfile:
            expected = [x.rsid for x in bfile.fetch('01')]
        self.assertGreater(len(expected), 1)
        
        with BgenReader(path) as bfile:
            rsids = []
            for var in bfile.fetch('01'):
                rsids.append(var.rsid)
                # each of these goes through the index too
                bfile.with_rsid('RSID_5')
                bfile.at_position(5000)
                bfile.rsids()
            self.assertEqual(rsids, expected)
    
    def test_index_rejects_queries_once_closed(self):
        ''' a closed index must not be usable '''
        index = Index(self.folder / 'example.16bits.bgen.bgi')
        list(index.fetch('01'))
        index.close()
        with self.assertRaises(ValueError):
            list(index.fetch('01'))
        with self.assertRaises(ValueError):
            index.offset_by_rsid('RSID_5')
        # closing twice must stay harmless
        index.close()
    
    def test_reader_uses_an_index_without_variants(self):
        ''' an indexed bgen with no variants is still read through its index '''
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'empty.bgen'
            with BgenWriter(path, n_samples=2):
                pass
            with BgenReader(path) as bfile:
                self.assertEqual(len(bfile), 0)
                self.assertEqual(list(bfile.fetch('01')), [])
                self.assertEqual(bfile.with_rsid('rs1'), [])
                self.assertEqual(bfile.at_position(1), [])
                self.assertEqual(bfile.rsids(), [])
                with self.assertRaises(ValueError):
                    bfile.varids()
    
    def test_reader_index_lists_are_fresh_each_time(self):
        ''' changing a list from an indexed reader does not change the next one '''
        with BgenReader(self.folder / 'example.16bits.bgen') as bfile:
            for name in ['rsids', 'chroms', 'positions']:
                with self.subTest(name=name):
                    values = getattr(bfile, name)()
                    expected = list(values)
                    values.append(None)
                    self.assertEqual(getattr(bfile, name)(), expected)
    
    def test_missing_index_is_not_created(self):
        ''' opening an index which does not exist fails, rather than making one '''
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'missing.bgen.bgi'
            with self.assertRaises(FileNotFoundError):
                Index(path)
            self.assertFalse(path.exists())
    
    def _positions_bgen(self, tmp):
        ''' a small bgen with its index, for comparing lookups by position '''
        path = Path(tmp) / 'positions.bgen'
        geno = np.full((2, 3), 1 / 3)
        with BgenWriter(path, 2) as bfile:
            for i, pos in enumerate([20, 10, 30, 10]):
                bfile.add_variant(f'v{i}', f'rs{i}', '01', pos, ['A', 'C'], geno)
        return path
    
    def test_position_lookups_match_without_an_index(self):
        ''' an indexed bgen finds the same variants by position as an unindexed one

        sqlite binds numpy scalars as blobs, which never equal a stored integer, and
        cannot bind a python int beyond 64 bits at all
        '''
        keys = [10, np.int64(10), np.int32(10), np.uint32(10), np.float32(10.0), 10.0,
                10.5, float('nan'), -1, 2 ** 40, 2 ** 63, '10', None]
        with tempfile.TemporaryDirectory() as tmp:
            path = self._positions_bgen(tmp)
            with BgenReader(path) as bfile:
                indexed = [[x.varid for x in bfile.at_position(k)] for k in keys]
            Path(str(path) + '.bgi').unlink()
            with BgenReader(path, delay_parsing=False) as bfile:
                unindexed = [[x.varid for x in bfile.at_position(k)] for k in keys]
        self.assertEqual(indexed, unindexed)
        self.assertEqual(indexed[1], ['v1', 'v3'])
    
    def test_fetch_bounds_from_numpy_and_huge_ints(self):
        ''' region bounds behave the same whatever numeric type gives them '''
        with tempfile.TemporaryDirectory() as tmp:
            path = self._positions_bgen(tmp)
            with BgenReader(path) as bfile:
                fetched = lambda *args: sorted(x.varid for x in bfile.fetch('01', *args))
                for start, stop in [(15, 25), (np.int64(15), np.int64(25)),
                                    (np.float32(15), np.uint32(25)), (15.0, 25.0)]:
                    with self.subTest(start=start, stop=stop):
                        self.assertEqual(fetched(start, stop), ['v0'])
                self.assertEqual(fetched(None, np.int64(0)), [])
                self.assertEqual(fetched(-2 ** 70, 2 ** 70), ['v0', 'v1', 'v2', 'v3'])
                self.assertEqual(fetched(2 ** 63), [])
    
    def test_lookups_come_in_file_order(self):
        ''' variants sharing an rsID or position come in file order, as without an index '''
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'order.bgen'
            geno = np.full((2, 3), 1 / 3)
            with BgenWriter(path, 2) as bfile:
                for i, (chrom, pos) in enumerate([('02', 10), ('01', 30), ('01', 10)]):
                    bfile.add_variant(f'v{i}', 'rsX', chrom, pos, ['A', 'C'], geno)
            for indexed in [True, False]:
                if not indexed:
                    Path(str(path) + '.bgi').unlink()
                with self.subTest(indexed=indexed):
                    with BgenReader(path, delay_parsing=indexed) as bfile:
                        self.assertEqual([x.varid for x in bfile.with_rsid('rsX')],
                                         ['v0', 'v1', 'v2'])
                        self.assertEqual([x.varid for x in bfile.at_position(10)],
                                         ['v0', 'v2'])
    
    def test_truncated_bgen_lookups_say_so(self):
        ''' an index listing variants past the end of a truncated bgen raises ValueError '''
        with tempfile.TemporaryDirectory() as tmp:
            path = self._positions_bgen(tmp)
            with BgenReader(path) as bfile:
                cut = bfile[2].fileoffset
            path.write_bytes(path.read_bytes()[:cut])
            with BgenReader(path) as bfile:
                for name, lookup in [('fetch', lambda: list(bfile.fetch('01'))),
                                     ('with_rsid', lambda: bfile.with_rsid('rs2')),
                                     ('at_position', lambda: bfile.at_position(30))]:
                    with self.subTest(name=name):
                        with self.assertRaisesRegex(ValueError, 'truncated'):
                            lookup()
