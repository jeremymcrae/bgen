
import sqlite3
import logging

import numpy as np

class Index:
    def __init__(self, path):
        logging.debug(f'opening bgen index: {path}')
        self.path = str(path)
        self.conn = sqlite3.connect(self.path)
        self.cur = self.conn.cursor()
        self.dropped_variants = None
        
        self._offsets: list[int] | None = None
        self._rsids: list[str] | None = None
        self._chroms: list[str] | None = None
        self._positions: list[int] | None = None
    
    def fetch(self, chrom, start=None, stop=None):
        ''' get file offsets for variants within a genome region in a bgen file
        
        Args:
            chrom: chromosome that variants must be on
            start: start nucleotide of region. If None, gets offsets for all
                variants on chromosome
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
        else:
            query = 'SELECT file_start_position FROM Variant \
                        WHERE chromosome=? AND position >= ? AND position <= ?'
            params = (chrom, start, stop)

        for row in self.cur.execute(query, params):
            yield row[0]
    
    def offset_by_index(self, index) -> int:
        ''' get file offset of bgen variant given a variant index
        '''
        if self._offsets is None:
            query = "SELECT file_start_position FROM Variant ORDER BY file_start_position"
            self._offsets = np.fromiter((x[0] for x in self.cur.execute(query)), dtype=np.uint64)
        return int(self._offsets[index])
        
    def offset_by_rsid(self, rsid) -> list[int]:
        ''' get file offset of bgen variant given a variant index
        '''
        query = "SELECT file_start_position FROM Variant WHERE rsid = ?"
        params = (rsid, )
        return [x[0] for x in self.cur.execute(query, params)]
    
    def offset_by_pos(self, pos) -> list[int]:
        ''' get file offset of bgen variant given a variant index
        '''
        query = "SELECT file_start_position FROM Variant WHERE position = ?"
        params = (pos, )
        return [x[0] for x in self.cur.execute(query, params)]
    
    @property
    def rsids(self):
        ''' get rsID list for all variants in the bgen file
        '''
        if self._rsids is None:
            query = "SELECT rsid FROM Variant ORDER BY file_start_position"
            self._rsids = [x[0] for x in self.cur.execute(query)]
        return self._rsids
    
    @property
    def chroms(self):
        ''' get chromosome list for all variants in the bgen file
        '''
        if self._chroms is None:
            query = "SELECT chromosome FROM Variant ORDER BY file_start_position"
            self._chroms = [x[0] for x in self.cur.execute(query)]
        return self._chroms
    
    @property
    def positions(self):
        ''' get position list for all variants in the bgen file
        '''
        if self._positions is None:
            query = "SELECT position FROM Variant ORDER BY file_start_position"
            self._positions = [x[0] for x in self.cur.execute(query)]
        return self._positions

    def close(self):
        if sqlite3 is None:
            # interpreter shutting down, nothing to clean up
            self.cur = None
            self.conn = None
            return
        if self.cur is not None:
            self.cur.close()
        self.cur = None
        if self.conn is not None:
            self.conn.close()
        self.conn = None
