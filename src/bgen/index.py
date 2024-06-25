
import bisect
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
        
        self._rsids = None
        self._chroms = None
        self._positions = None
        self._masked = None
    
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

        query = self.cur.execute(query, params)
        while True:
            res = query.fetchone()
            if not res:
                break
            yield res[0]
    
    def offset_by_index(self, index):
        ''' get file offset of bgen variant given a variant index
        '''
        offset = self.cur.execute('''SELECT file_start_position FROM Variant LIMIT 1 OFFSET ?''', (index, )).fetchone()
        return offset[0]
        
    def offset_by_rsid(self, rsid):
        ''' get file offset of bgen variant given a variant index
        '''
        offsets = self.cur.execute("SELECT file_start_position FROM Variant WHERE rsid= ?", (rsid, )).fetchall()
        return [x[0] for x in offsets]
    
    def offset_by_pos(self, pos):
        ''' get file offset of bgen variant given a variant index
        '''
        offsets = self.cur.execute("SELECT file_start_position FROM Variant WHERE position= ?", (pos, )).fetchall()
        return [x[0] for x in offsets]
    
    @property
    def rsids(self):
        ''' get rsID list for all variants in the bgen file
        '''
        if self._rsids is None:
            query = self.cur.execute(
                "SELECT rsid FROM Variant ORDER BY file_start_position")
            self._rsids = [x[0] for x in query.fetchall()]
        return self._rsids
    
    @property
    def chroms(self):
        ''' get chromosome list for all variants in the bgen file
        '''
        if self._chroms is None:
            query = self.cur.execute(
                "SELECT chromosome FROM Variant ORDER BY file_start_position")
            self._chroms = [x[0] for x in query.fetchall()]
        return self._chroms
    
    @property
    def positions(self):
        ''' get position list for all variants in the bgen file
        '''
        if self._positions is None:
            query = self.cur.execute(
                "SELECT position FROM Variant ORDER BY file_start_position")
            self._positions = np.array([x[0] for x in query.fetchall()])
        return self._positions

    def close(self):
        if self.cur is not None:
            self.cur.close()
        self.cur = None
        if self.conn is not None:
            self.conn.close()
        self.conn = None
