
import bisect
import sqlite3
import logging

import numpy as np

class Index:
    def __init__(self, path):
        logging.info(f'opening bgen index: {path}')
        self.path = str(path)
        self.dropped_variants = None
        
        self._rsids = None
        self._chroms = None
        self._positions = None
        self._masked = None
    
    def offset_by_index(self, index):
        ''' get file offset of bgen variant given a variant index
        '''
        with sqlite3.connect(self.path) as conn:
            offset = conn.execute('''SELECT file_start_position FROM Variant LIMIT 1 OFFSET ?''', (index, )).fetchone()
        return offset[0]
        
    def offset_by_rsid(self, rsid):
        ''' get file offset of bgen variant given a variant index
        '''
        with sqlite3.connect(self.path) as conn:
            offsets = conn.execute("SELECT file_start_position FROM Variant WHERE rsid= ?", (rsid, )).fetchall()
        
        if len(offsets) == 0:
            raise ValueError(f'cannot find variant match for {rsid}')
        elif len(offsets) > 1:
            raise ValueError(f'multiple variant matches for {rsid}')
        
        return offsets[0][0]
    
    def offset_by_pos(self, pos):
        ''' get file offset of bgen variant given a variant index
        '''
        with sqlite3.connect(self.path) as conn:
            offsets = conn.execute("SELECT file_start_position FROM Variant WHERE position= ?", (pos, )).fetchall()
        
        if len(offsets) == 0:
            raise ValueError(f'cannot find variant match at pos: {pos}')
        elif len(offsets) > 1:
            raise ValueError(f'multiple variant matches at pos: {pos}')
        
        return offsets[0][0]
    
    @property
    def rsids(self):
        ''' get rsID list for all variants in the bgen file
        '''
        if self._rsids is None:
            with sqlite3.connect(self.path) as conn:
                query = conn.execute("SELECT rsid FROM Variant")
            self._rsids = [x[0] for x in query.fetchall()]
        return self._rsids
    
    @property
    def chroms(self):
        ''' get chromosome list for all variants in the bgen file
        '''
        if self._chroms is None:
            with sqlite3.connect(self.path) as conn:
                query = conn.execute("SELECT chromosome FROM Variant")
            self._chroms = [x[0] for x in query.fetchall()]
        return self._chroms
    
    @property
    def positions(self):
        ''' get position list for all variants in the bgen file
        '''
        if self._positions is None:
            with sqlite3.connect(self.path) as conn:
                query = conn.execute("SELECT position FROM Variant")
            self._positions = np.array([x[0] for x in query.fetchall()])
        return self._positions