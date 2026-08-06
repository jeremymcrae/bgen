
import sqlite3
import logging
import os
from typing import Optional, Union

import numpy as np
from numpy.typing import NDArray

class Index:
    def __init__(self, path: Union[str, os.PathLike[str]]) -> None:
        logging.debug(f'opening bgen index: {path}')
        self.path = str(path)
        self.conn = sqlite3.connect(self.path)
        
        self._offsets: Optional[NDArray[np.uint64]] = None
        self._rsids: Optional[list[str]] = None
        self._chroms: Optional[list[str]] = None
        self._positions: Optional[list[int]] = None
    
    def _query(self, query, params=()):
        ''' run a query on a cursor of its own
        
        Each query needs its own cursor, since re-executing one abandons the rows
        it had not yet handed out. fetch() yields as it steps through its results,
        so sharing a cursor let any other query cut a fetch short partway.
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
        if start is None and stop is None:
            query = 'SELECT file_start_position FROM Variant WHERE chromosome=?'
            params = (chrom, )
        elif stop is None:
            query = 'SELECT file_start_position FROM Variant \
                        WHERE chromosome=? AND position >= ?'
            params = (chrom, start)
        elif start is None:
            # a stop without a start needs its own query, since comparing a
            # position against a NULL start is never true in sqlite, so folding
            # this into the query below would silently match nothing
            query = 'SELECT file_start_position FROM Variant \
                        WHERE chromosome=? AND position <= ?'
            params = (chrom, stop)
        else:
            query = 'SELECT file_start_position FROM Variant \
                        WHERE chromosome=? AND position >= ? AND position <= ?'
            params = (chrom, start, stop)

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
        ''' get file offset of bgen variant given a variant index
        '''
        query = "SELECT file_start_position FROM Variant WHERE rsid = ?"
        params = (rsid, )
        return [x[0] for x in self._query(query, params)]
    
    def offset_by_pos(self, pos) -> list[int]:
        ''' get file offset of bgen variant given a variant index
        '''
        query = "SELECT file_start_position FROM Variant WHERE position = ?"
        params = (pos, )
        return [x[0] for x in self._query(query, params)]
    
    @property
    def rsids(self):
        ''' get rsID list for all variants in the bgen file
        '''
        if self._rsids is None:
            query = "SELECT rsid FROM Variant ORDER BY file_start_position"
            self._rsids = [x[0] for x in self._query(query)]
        return self._rsids
    
    @property
    def chroms(self):
        ''' get chromosome list for all variants in the bgen file
        '''
        if self._chroms is None:
            query = "SELECT chromosome FROM Variant ORDER BY file_start_position"
            self._chroms = [x[0] for x in self._query(query)]
        return self._chroms
    
    @property
    def positions(self):
        ''' get position list for all variants in the bgen file
        '''
        if self._positions is None:
            query = "SELECT position FROM Variant ORDER BY file_start_position"
            self._positions = [x[0] for x in self._query(query)]
        return self._positions

    def close(self):
        if sqlite3 is None:
            # interpreter shutting down, nothing to clean up
            self.conn = None
            return
        if self.conn is not None:
            self.conn.close()
        self.conn = None
