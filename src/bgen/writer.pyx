# cython: language_level=3, boundscheck=False, emit_linenums=True

import logging
import os
from pathlib import Path
import sqlite3
import sys
import time

from libcpp cimport bool
from libcpp.string cimport string
from libcpp.vector cimport vector
from libc.stdint cimport int32_t, uint8_t, uint16_t, uint32_t, uint64_t

import numpy as np

_IS_WIN32 = sys.platform == 'win32'

# bits per stored probability unless the caller asks for something else. This gives a
# step of 1/255, so two decimal places survive a round trip, and it doubles as the
# yardstick for warning about depths that lose more than this
DEFAULT_BIT_DEPTH = 8

cdef extern from 'writer.h' namespace 'bgen':
    cdef cppclass CppBgenWriter:
        # declare class constructor and methods
        CppBgenWriter(const string &path, uint32_t n_samples, const string &free_data,
                   uint32_t compression, uint32_t layout, const vector[string] &samples) except +
        uint64_t write_variant_header(const string &varid, const string &rsid,
                        const string &chrom, uint32_t pos, const vector[string] &alleles) except +
        uint64_t write_variant_direct(const char * data, size_t size) except +
        void encode_genotype_data(uint16_t n_alleles,
                         double *genotypes, uint32_t geno_len, const uint8_t *ploidy,
                         uint8_t min_ploidy, uint8_t max_ploidy,
                         bool phased, uint8_t bit_depth) except +
        uint64_t write_genotype_data() except +
        void close() except +

class Indexer:
    ''' class to automatically index bgen files as they are being constructed
    '''
    def __init__(self, bgen_path):
        self.index_path = Path(str(bgen_path) + '.bgi')
        if self.index_path.exists():
            self.index_path.unlink()
        self.create_time = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())
        self.conn = sqlite3.connect(self.index_path)
        self.cur = self.conn.cursor()
        self.create_tables()
    
    def create_tables(self):
        query = '''CREATE TABLE Metadata (
                    filename TEXT NOT NULL, 
                    file_size INT NOT NULL, 
                    last_write_time INT NOT NULL, 
                    first_1000_bytes BLOB NOT NULL, 
                    index_creation_time INT NOT NULL)'''
        self.cur.execute(query)
        query = '''CREATE TABLE Variant (
                    chromosome TEXT NOT NULL,
                    position INT NOT NULL,
                    rsid TEXT NOT NULL,
                    number_of_alleles INT NOT NULL,
                    allele1 TEXT NOT NULL,
                    allele2 TEXT NULL,
                    file_start_position INT NOT NULL,
                    size_in_bytes INT NOT NULL,
                PRIMARY KEY (chromosome, position, rsid, allele1, allele2, file_start_position))
                WITHOUT ROWID'''
        self.cur.execute(query)
        
        # index the Variant table
        self.cur.execute('CREATE INDEX chrom_index on Variant(chromosome)')
        self.cur.execute('CREATE INDEX pos_index on Variant(position)')
        self.cur.execute('CREATE INDEX rsid_index on Variant(rsid)')
    
    def add_variant(self, chrom, pos, rsid, alleles, offset, size):
        allele_1 = alleles[0]
        # allele2 is part of the primary key, which cannot be NULL, so use '' if absent
        allele_2 = alleles[1] if len(alleles) > 1 else ''
        query = '''INSERT INTO Variant VALUES (?, ?, ?, ?, ?, ?, ?, ?)'''
        params = (chrom, pos, rsid, len(alleles), allele_1, allele_2, offset, size)
        self.cur.execute(query, params)
    
    def add_metadata(self):
        bgen_path = self.index_path.with_suffix('')
        bgen_size = bgen_path.stat().st_size
        bgen_time = int(bgen_path.stat().st_mtime)
        with open(bgen_path, 'rb') as f:
            first_1000_bytes = f.read(1000)
        query = '''INSERT INTO Metadata VALUES (?, ?, ?, ?, ?)'''
        params = (str(bgen_path), bgen_size, bgen_time, first_1000_bytes, self.create_time)
        self.cur.execute(query, params)

    def close(self):
        try:
            self.add_metadata()
        finally:
            self.conn.commit()
            if _IS_WIN32 and time is not None:
                time.sleep(0.01)
            self.cur.close()
            self.conn.close()

cdef class BgenWriter:
    ''' class to write bgen files to disk
    '''
    cdef CppBgenWriter * thisptr
    cdef string path
    cdef bool is_open
    cdef object indexer
    cdef int layout
    cdef object compression
    cdef uint32_t n_samples
    def __cinit__(self, path, uint32_t n_samples, samples=None, compression='zstd',
                  int layout=2, metadata=None):
        # any path-like object can name the file
        path = os.fspath(path)
        
        if compression not in [None, 'zstd', 'zlib']:
            raise ValueError(f'compression type {compression} not one of zlib or zstd')
        
        cdef uint32_t compress_flag=0
        if compression == 'zlib':
            compress_flag = 1
        elif compression == 'zstd':
            compress_flag = 2
        
        if layout not in [1, 2]:
            raise ValueError(f'layout must be 1 or 2: {layout}')
        
        if layout == 1 and compression == 'zstd':
            raise ValueError('layout 1 is not supported with zstd compression')
        
        self.n_samples = n_samples
        self.layout = layout
        self.compression = compression

        # re-define variables into cpp objects
        cdef string _metadata = metadata.encode('utf8') if metadata is not None else b''
        cdef vector[string] _samples = [x.encode('utf8') for x in samples] if samples is not None else []

        self.path = path.encode('utf8')
        
        logging.debug(f'opening CppBgenWriter from {self.path.decode("utf")}')
        self.thisptr = new CppBgenWriter(self.path, n_samples, _metadata, compress_flag, layout, _samples)
        self.is_open = True
        self.indexer = Indexer(path)
    
    def __dealloc__(self):
        # a deallocated writer has nobody left to report errors to, so drop them
        # rather than emitting 'Exception ignored in __dealloc__' noise. Call
        # close() directly (or use a with block) to have errors raised.
        try:
            self.close()
        except Exception:
            pass
    
    def __repr__(self):
        return f'BgenWriter("{self.path.decode("utf8")}")'
    

    def _check_bit_depth_fits(self, genotypes, bit_depth, rsid, varid):
        ''' warn when the bit depth is too coarse for the genotype probabilities

        Stored probabilities step by 1/(2**bit_depth - 1), so a low depth can lose
        information e.g. [0.3, 0.3, 0.4] encoded into a 1 bit file roundtrips to
        [0, 1, 0]. This check permits hard-coded genotypes, and variable probabilities
        by adjusting the tolerance to the bit depth.
        '''
        if bit_depth >= DEFAULT_BIT_DEPTH:
            # only coarser depths than the default need checking
            return
        tolerance = 0.5 / (2 ** DEFAULT_BIT_DEPTH - 1)

        factor = float(2 ** bit_depth - 1)
        with np.errstate(invalid='ignore'):
            error = np.abs(np.round(genotypes * factor) / factor - genotypes)
        if error.size == 0 or np.all(np.isnan(error)):
            # nothing was stored, or every sample is missing, so nothing was lost
            return
        worst = float(np.nanmax(error))
        if worst <= tolerance:
            return

        stores = ('only stores 0 or 1' if bit_depth == 1
                  else f'only stores multiples of 1/{factor:.0f}')
        logging.warning(f'variant {rsid}/{varid} loses up to {worst:.3g} per '
                        f'probability at bit_depth={bit_depth}, which {stores}. Use a '
                        f'higher bit_depth to keep the values you passed in')

    def _validate_genotypes(self, genotypes):
        ''' check the genotype values
        '''
        # ensure the genotypes array is the correct size and type
        genotypes = np.asarray(genotypes, dtype=np.float64)
        if genotypes.ndim != 2:
            raise ValueError('genotypes must be a 2D array')
        
        if genotypes.size == 0:
            raise ValueError('genotypes array must be non-empty')
        
        n_samples = genotypes.shape[0]
        if n_samples != self.n_samples:
            raise ValueError(f'genotypes array must have {self.n_samples} samples, not {n_samples}')
        
        return genotypes

    def _validate_ploidy(self, ploidy, n_samples):
        ''' check the ploidy values
        '''
        if isinstance(ploidy, list):
            ploidy = np.array(ploidy)
        
        # determine ploidy levels
        cdef int32_t min_ploidy, max_ploidy
        cdef const uint8_t[::1] ploidy_arr = None
        if np.isscalar(ploidy) and np.issubdtype(type(ploidy), np.integer):
            min_ploidy, max_ploidy = ploidy, ploidy
        elif isinstance(ploidy, np.ndarray) and np.issubdtype(ploidy.dtype, np.integer):
            if ploidy.ndim != 1 or ploidy.size != n_samples:
                raise ValueError(f'ploidy array must have {n_samples} entries, not {ploidy.size}')

            min_ploidy, max_ploidy = np.min(ploidy), np.max(ploidy)
            if min_ploidy != max_ploidy:
                # contiguous, since the encoder reads through a plain pointer
                ploidy_arr = np.ascontiguousarray(ploidy, dtype=np.uint8)
        else:
            raise ValueError('ploidy must be either integer, or numpy array of integers')
        
        if min_ploidy < 0 or max_ploidy > 63:
            raise ValueError('ploidy values must be between 0 and 63')
        
        return min_ploidy, max_ploidy, ploidy_arr

    def _validate_layout1_data(self, n_alleles, n_genos, phased):
        ''' validate layout 1 data
        '''
        if self.layout != 1:
            return
        if n_alleles != 2:
            raise ValueError('layout 1 requires exactly two alleles')
        elif n_genos != 3:
            raise ValueError('layout 1 requires 3 genotype probabilities per variant')
        elif phased:
            raise ValueError('layout 1 cannot use phased data')

    def add_variant(self, varid, rsid, chrom, uint32_t pos, alleles, 
                    genotypes, ploidy=2, bool phased=False,
                    int bit_depth=DEFAULT_BIT_DEPTH):
        ''' add a variant to the bgen file on disk

        Args:
            varid: variant ID
            rsid: reference SNP ID
            chrom: chromosome the variant is on
            pos: nucleotide position of the variant
            alleles: list of allele strings. Duplicates are allowed, but warn, since
                the minor allele cannot then say which one it means.
            genotypes: numpy array of genotype proabilities, ordered as per the
                bgen samples.
            ploidy: integer for constant ploidy, or numpy array of ploidy values per 
                sample, in same order as genotypes
            phased: whether the genotypes are for phased data or not
            bit_depth: integer from 1-32 (inclusive) for how many bits to store
                each genotype in.
        '''
        if not self.is_open:
            raise ValueError("bgen file is closed")

        # re-define variables into cpp objects
        cdef string _varid = varid.encode('utf8')
        cdef string _rsid = rsid.encode('utf8')
        cdef string _chrom = chrom.encode('utf8')

        alleles = list(alleles)
        if len(alleles) == 0:
            raise ValueError('alleles must be a non-empty list')
        if len(set(alleles)) != len(alleles):
            # bgen spec allows this, and genotypes will be fine, but probably a bad input
            logging.warning(f'variant {rsid}/{varid} has duplicate alleles {alleles}')
        cdef vector[string] _alleles = [x.encode('utf8') for x in alleles]

        if bit_depth < 1 or bit_depth > 32:
            raise ValueError(f'bit_depth must be between 1 and 32: {bit_depth}')
        
        # sanitise the genotypes. They are made C contiguous (e.g. a transposed array is
        # not), since the encoder reads them through a plain pointer
        genotypes = self._validate_genotypes(genotypes)
        cdef uint32_t n_genos = genotypes.shape[1]
        self._check_bit_depth_fits(genotypes, bit_depth, rsid, varid)
        cdef const double[:, ::1] geno_c = np.ascontiguousarray(genotypes)
    
        # validate ploidy levels
        cdef int32_t min_ploidy, max_ploidy
        cdef const uint8_t[::1] ploidy_arr
        min_ploidy, max_ploidy, ploidy_arr = self._validate_ploidy(ploidy, self.n_samples)
        
        self._validate_layout1_data(len(alleles), n_genos, phased)
        
        # the ploidy array is only passed on when the ploidy varies
        cdef const uint8_t * ploidy_ptr = NULL
        if min_ploidy != max_ploidy:
            ploidy_ptr = &ploidy_arr[0]
        
        # encode the genotypes before writing anything, so that a variant with
        # invalid genotypes cannot leave a partial variant in the bgen file
        cdef uint32_t geno_len = self.n_samples * n_genos
        self.thisptr.encode_genotype_data(_alleles.size(), <double *> &geno_c[0, 0],
                                          geno_len, ploidy_ptr, min_ploidy, max_ploidy,
                                          phased, bit_depth)
        
        var_offset = self.thisptr.write_variant_header(_varid, _rsid, _chrom, pos, _alleles)
        end_offset = self.thisptr.write_genotype_data()
        
        self.indexer.add_variant(chrom, pos, rsid, alleles, var_offset,
                                 end_offset - var_offset)

    def add_variant_direct(self, variant):
        ''' insert a BgenVar directly into the bgen file

        This copies the variant's genotype data across as raw bytes, without
        decoding it, so the destination bgen has to be able to hold that data
        unchanged. That means the sample count, layout and compression scheme all
        have to match the bgen the variant came from, otherwise the copied bytes
        would be misinterpreted when read back.

        Args:
            variant: a BgenVar from a bgen file with the same samples, layout and
                compression as this file
        '''
        if not self.is_open:
            raise ValueError("bgen file is closed")

        self._check_compatible(variant)

        chrom = variant.chrom
        pos = int(variant.pos)
        rsid = variant.rsid
        alleles = variant.alleles
        raw = variant.copy_data()
        cdef const uint8_t[::1] data = raw
        if data.shape[0] == 0:
            raise ValueError('cannot copy variant data directly: the variant is empty')
        var_offset = self.thisptr.write_variant_direct(<const char *> &data[0], data.shape[0])

        self.indexer.add_variant(chrom, pos, rsid, alleles, var_offset, len(raw))

    def _check_compatible(self, variant):
        ''' check a variant's data can be copied into this file unchanged
        '''
        n_samples = getattr(variant, 'n_samples', None)
        layout = getattr(variant, 'layout', None)
        compression = getattr(variant, 'compression', None)
        if n_samples is None or layout is None:
            # not a BgenVar, so let copy_data() report the problem
            return

        if n_samples != self.n_samples:
            raise ValueError(f'cannot copy variant data directly: the variant has '
                             f'{n_samples} samples, but this bgen has {self.n_samples}')
        if layout != self.layout:
            raise ValueError(f'cannot copy variant data directly: the variant uses '
                             f'layout {layout}, but this bgen uses layout {self.layout}')
        if compression != self.compression:
            raise ValueError(f'cannot copy variant data directly: the variant uses '
                             f'{compression} compression, but this bgen uses '
                             f'{self.compression}')

    def __enter__(self):
        return self
    
    def __exit__(self, exc_type, exc_value, traceback):
        self.close()
        return False
    
    def close(self):
        try:
            if self.is_open:
                self.is_open = False
                try:
                    # close explicitly, so that errors while finishing the file
                    # off (e.g. a full disk) are raised rather than being
                    # swallowed by the C++ destructor
                    self.thisptr.close()
                finally:
                    del self.thisptr
        finally:
            if self.indexer is not None:
                self.indexer.close()
            self.indexer = None
            if _IS_WIN32 and time is not None:
                time.sleep(0.01)
