
import bisect
import sqlite3
import logging

import numpy as np

class Index:
    def __init__(self, path):
        logging.info(f'opening bgen index: {path}')
        self.conn = sqlite3.connect(path)
        self.cursor = self.conn.cursor()
        self.dropped_variants = None
        
        self._rsids = None
        self._chroms = None
        self._positions = None
        self._offsets = None
        self._masked = None
    
    def offset_by_index(self, index):
        ''' get file offset of bgen variant given a variant index
        '''
        return self.offsets[index]
        
    def offset_by_rsid(self, rsid):
        ''' get file offset of bgen variant given a variant index
        '''
        i = bisect.bisect_left(self._rsids_idx, (rsid, 0))
        indices = []
        while i < len(self._rsids) and self._rsids_idx[i][0] == rsid:
            indices.append(self._rsids_idx[i][1])
            i += 1
        
        if len(indices) == 0:
            raise ValueError(f'cannot find variant match for {rsid}')
        elif len(indices) > 1:
            raise ValueError(f'multiple variant matches for {rsid}')
        
        idx = indices[0]
        return self.offsets[idx]
    
    def offset_by_pos(self, pos):
        ''' get file offset of bgen variant given a variant index
        '''
        i = bisect.bisect_left(self._positions_idx, (pos, 0))
        indices = []
        while i < len(self._positions_idx) and self._positions_idx[i][0] == pos:
            indices.append(self._positions_idx[i][1])
            i += 1
        
        if len(indices) == 0:
            raise ValueError(f'cannot find variant match for {pos}')
        elif len(indices) > 1:
            raise ValueError(f'multiple variant matches for {pos}')
        
        idx = indices[0]
        return self.offsets[idx]
    
    @property
    def rsids(self):
        ''' get rsID list for all variants in the bgen file
        '''
        if self._rsids is None:
            self.cursor.execute("SELECT rsid FROM Variant")
            self._rsids = [x[0] for x in self.cursor.fetchall()]
            self._rsids_idx = sorted(zip(self._rsids, range(len(self._rsids))))
        return self._rsids
    
    @property
    def chroms(self):
        ''' get chromosome list for all variants in the bgen file
        '''
        if self._chroms is None:
            self.cursor.execute("SELECT chromosome FROM Variant")
            self._chroms = [x[0] for x in self.cursor.fetchall()]
        return self._chroms
    
    @property
    def positions(self):
        ''' get position list for all variants in the bgen file
        '''
        if self._positions is None:
            self.cursor.execute("SELECT position FROM Variant")
            self._positions = np.array([x[0] for x in self.cursor.fetchall()])
            self._positions_idx = sorted(zip(self._positions, range(len(self._positions))))
        return self._positions
    
    @property
    def offsets(self):
        ''' get variant file offsets for all variants in the bgen file
        '''
        if self._offsets is None:
            self.cursor.execute("SELECT file_start_position FROM Variant")
            self._offsets = np.array([x[0] for x in self.cursor.fetchall()])
        return self._offsets
