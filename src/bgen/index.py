
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
        
    def offset_by_rsid(self, rsid) -> list[int]:
        ''' get file offsets of the bgen variants with a given rsID
        '''
        query = "SELECT file_start_position FROM Variant WHERE rsid = ?"
        params = (rsid, )
        return [x[0] for x in self._query(query, params)]
    
    def offset_by_pos(self, pos) -> list[int]:
        ''' get file offsets of the bgen variants at a given position
        '''
        pos = _whole_number(pos)
        if pos is None or not 0 <= pos <= MAX_POSITION:
            # no stored position can equal this
            return []
        query = "SELECT file_start_position FROM Variant WHERE position = ?"
        params = (pos, )
        return [x[0] for x in self._query(query, params)]
    
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
