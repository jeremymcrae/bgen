
from pathlib import Path
import tempfile
import unittest
from unittest import mock


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
        bfile = BgenReader(self.folder / 'example.15bits.bgen')
        self.assertFalse(bfile._check_for_index(str(self.folder / 'example.15bits.bgen')))
        
        bfile = BgenReader(self.folder / 'example.16bits.bgen')
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
