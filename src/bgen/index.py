
from functools import cached_property
import numbers
import sqlite3
import logging
import os
from typing import Optional, Union

import numpy as np
from numpy.typing import NDArray

# positions are unsigned 32-bit
MAX_POSITION = 0xFFFFFFFF

def _whole_number(value) -> Optional[int]:
    ''' value as a python int, or None if it is not a whole number

    numpy scalars are converted too, since sqlite binds them as blobs, which never
    compare equal to a stored integer
    '''
    if isinstance(value, numbers.Integral):
        return int(value)
    if isinstance(value, numbers.Number):
        try:
            whole = int(value)
        except (TypeError, ValueError, OverflowError):
            # complex, nan or infinite
            return None
        return whole if whole == value else None
    return None

# what a lookup key can match: nothing, one value to look up exactly, or whatever the
# stored values say they equal, for a key whose == only it knows
NO_MATCH, EXACT, COMPARE = range(3)

def lookup_key(key, by_rsid):
    ''' how to look up variants by a key, so that every lookup matches as == would

    rsIDs are stored as strings, and positions as unsigned 32-bit integers. Comparing
    other built in types against those can never succeed, whereas sqlite would convert
    them, so that 1 matched the rsID '1'

    Args:
        key: rsID or position to look up
        by_rsid: whether the key is an rsID, rather than a position

    Returns:
        tuple of NO_MATCH, EXACT or COMPARE, and the value to look up for EXACT
    '''
    if by_rsid:
        if isinstance(key, str):
            try:
                key.encode('utf8')
            except UnicodeEncodeError:
                # a stored rsID is valid utf8, so cannot equal this
                return NO_MATCH, None
            return EXACT, key
        if isinstance(key, (numbers.Number, bytes, type(None))):
            return NO_MATCH, None
        return COMPARE, key
    if isinstance(key, numbers.Number):
        pos = _whole_number(key)
        if pos is None or not 0 <= pos <= MAX_POSITION:
            return NO_MATCH, None
        return EXACT, pos
    if isinstance(key, (str, bytes, type(None))):
        return NO_MATCH, None
    return COMPARE, key

def _bound(value):
    ''' a region bound as a python number which sqlite compares with positions

    numpy scalars would otherwise bind as blobs, and an int beyond 64 bits cannot bind
    at all. Positions are 32-bit, so clamping just outside that range changes nothing
    '''
    if isinstance(value, numbers.Integral):
        value = int(value)
    elif isinstance(value, numbers.Real):
        value = float(value)
    else:
        return value
    if value != value:
        # nan, which no position satisfies
        return value
    return min(max(value, -1), MAX_POSITION + 1)

class Index:
    def __init__(self, path: Union[str, os.PathLike[str]]) -> None:
        logging.debug(f'opening bgen index: {path}')
        self.path = str(path)
        # sqlite would otherwise create an empty database at a missing path
        if not os.path.isfile(self.path):
            raise FileNotFoundError(f'bgen index not found: {self.path}')
        self.conn = sqlite3.connect(self.path)
        
        self._offsets: Optional[NDArray[np.uint64]] = None
    
    def _query(self, query, params=()):
        ''' run a query on a cursor of its own
        
        Each query needs its own cursor, since re-executing one abandons the rows
        it had not yet handed out, and fetch() yields as it steps through its results.
        '''
        if self.conn is None:
            raise ValueError('bgen index is closed')
        cur = self.conn.cursor()
        try:
            yield from cur.execute(query, params)
        finally:
            cur.close()
    
    def fetch(self, chrom, start=None, stop=None):
        ''' get file offsets for variants within a genome region in a bgen file
        
        Args:
            chrom: chromosome that variants must be on
            start: start nucleotide of region. If None, gets offsets for
                variants with positions up to stop, or for all variants on the
                chromosome if stop is also None
            stop: end nucleotide of region. If None, gets offsets for variants
                with positions after start
        
        Yields:
            file offsets for variants within the genome region
        '''
        # only bound the ends which are given, since comparing a position against a
        # NULL bound is never true in sqlite
        query = 'SELECT file_start_position FROM Variant WHERE chromosome=?'
        params = [chrom]
        if start is not None:
            query += ' AND position >= ?'
            params.append(_bound(start))
        if stop is not None:
            query += ' AND position <= ?'
            params.append(_bound(stop))

        for row in self._query(query, params):
            yield row[0]
    
    def _load_offsets(self) -> NDArray[np.uint64]:
        ''' get file offsets of every variant in the bgen, in file order
        '''
        if self._offsets is None:
            query = "SELECT file_start_position FROM Variant ORDER BY file_start_position"
            self._offsets = np.fromiter((x[0] for x in self._query(query)), dtype=np.uint64)
        return self._offsets
    
    def __len__(self) -> int:
        ''' number of variants listed in the index
        '''
        return len(self._load_offsets())
    
    def offset_by_index(self, index) -> int:
        ''' get file offset of bgen variant given a variant index
        '''
        return int(self._load_offsets()[index])
        
    def _offsets_matching(self, key, by_rsid) -> list[int]:
        ''' file offsets of the variants whose rsID, or else position, equals a key
        
        In file order, as a bgen without an index gives them
        '''
        kind, value = lookup_key(key, by_rsid)
        if kind == NO_MATCH:
            return []
        if kind == COMPARE:
            values = self.rsids if by_rsid else self.positions
            offsets = self._load_offsets()
            return [int(offsets[i]) for i, x in enumerate(values) if x == key]
        column = 'rsid' if by_rsid else 'position'
        query = (f"SELECT file_start_position FROM Variant WHERE {column} = ? "
                 "ORDER BY file_start_position")
        return [x[0] for x in self._query(query, (value, ))]
    
    def offset_by_rsid(self, rsid) -> list[int]:
        ''' get file offsets of the bgen variants with a given rsID
        '''
        return self._offsets_matching(rsid, True)
    
    def offset_by_pos(self, pos) -> list[int]:
        ''' get file offsets of the bgen variants at a given position
        '''
        return self._offsets_matching(pos, False)
    
    @cached_property
    def rsids(self) -> list[str]:
        ''' get rsID list for all variants in the bgen file
        '''
        query = "SELECT rsid FROM Variant ORDER BY file_start_position"
        return [x[0] for x in self._query(query)]
    
    @cached_property
    def chroms(self) -> list[str]:
        ''' get chromosome list for all variants in the bgen file
        '''
        query = "SELECT chromosome FROM Variant ORDER BY file_start_position"
        return [x[0] for x in self._query(query)]
    
    @cached_property
    def positions(self) -> list[int]:
        ''' get position list for all variants in the bgen file
        '''
        query = "SELECT position FROM Variant ORDER BY file_start_position"
        return [x[0] for x in self._query(query)]

    def close(self):
        if sqlite3 is None:
            # interpreter shutting down, nothing to clean up
            self.conn = None
            return
        if self.conn is not None:
            self.conn.close()
        self.conn = None
