# cython: language_level=3, boundscheck=False, emit_linenums=True

import logging
from pathlib import Path
import sqlite3
import sys
import time

from libcpp cimport bool
from libcpp.string cimport string
from libcpp.vector cimport vector
from libc.stdint cimport int8_t, int32_t, uint8_t, uint16_t, uint32_t, uint64_t

from cython.operator cimport dereference as deref

import numpy as np

_IS_WIN32 = sys.platform == 'win32'

cdef extern from "<iostream>" namespace "std":
    cdef cppclass ostream:
        pass

cdef extern from "<iostream>" namespace "std::ios_base":
    cdef cppclass open_mode:
        pass
    cdef open_mode binary

cdef extern from "<fstream>" namespace "std":
    cdef cppclass ofstream(ostream):
        ofstream() except +
        ofstream(const string&) except +
        ofstream(const string&, open_mode) except +

cdef extern from 'writer.h' namespace 'bgen':
    cdef cppclass CppBgenWriter:
        # declare class constructor and methods
        CppBgenWriter(string &path, uint32_t n_samples, string &free_data, 
                   uint32_t compression, uint32_t layout, vector[string] &samples) except +
        uint64_t write_variant_header(string &varid, string &rsid, string &chrom, 
                        uint32_t &pos, vector[string] &alleles, uint32_t _n_samples) except +
        uint64_t write_variant_direct(vector[uint8_t] & data) except +
        void encode_genotype_data(uint16_t n_alleles,
                         double *genotypes, uint32_t geno_len, uint8_t ploidy,
                         bool phased, uint8_t bit_depth) except +
        void encode_genotype_data(uint16_t n_alleles,
                         double *genotypes, uint32_t geno_len, uint8_t *ploidy,
                         uint8_t min_ploidy, uint8_t max_ploidy,
                         bool phased, uint8_t bit_depth) except +
        uint64_t write_genotype_data() except +

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
        allele_2 = alleles[1] if len(alleles) > 1 else None
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
    cdef uint32_t n_samples
    def __cinit__(self, path, uint32_t n_samples, samples=[], compression='zstd',
                  int layout=2, metadata=None):
        if isinstance(path, Path):
            path = str(path)
        
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

        # re-define variables into cpp objects
        cdef string _metadata = metadata.encode('utf8') if metadata is not None else b''
        cdef vector[string] _samples = [x.encode('utf8') for x in samples]

        self.path = path.encode('utf8')
        
        logging.debug(f'opening CppBgenWriter from {self.path.decode("utf")}')
        self.thisptr = new CppBgenWriter(self.path, n_samples, _metadata, compress_flag, layout, _samples)
        self.is_open = True
        self.indexer = Indexer(path)
    
    def __dealloc__(self):
        self.close()
    
    def __repr__(self):
        return f'BgenFile("{self.path.decode("utf8")}")'
    

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
        n_genos = genotypes.shape[1]
        if n_samples != self.n_samples:
            raise ValueError(f'genotypes array must have {self.n_samples} samples, not {n_samples}')
        
        return genotypes, n_samples, n_genos
    
    def _make_contiguous(self, genotypes):
        ''' make genotypes array C contiguous
        '''
        # convert numpy array to C contiguous for storing values on disk. numpy
        # arrays default to C contiguous, so most won't need conversion, but
        # some can be fortran order, e.g. if transposed
        cdef double[:, :] geno_c
        if genotypes.flags['C_CONTIGUOUS']:
            geno_c = genotypes
        else:
            geno_c = np.ascontiguousarray(genotypes)
        
        return geno_c

    def _validate_ploidy(self, ploidy, n_samples):
        ''' check the ploidy values
        '''
        if isinstance(ploidy, list):
            ploidy = np.array(ploidy)
        
        # determine ploidy levels
        cdef int32_t min_ploidy, max_ploidy
        cdef uint8_t[:] ploidy_arr = np.array([], dtype=np.uint8)
        if np.isscalar(ploidy) and np.issubdtype(type(ploidy), np.integer):
            min_ploidy, max_ploidy = ploidy, ploidy
        elif isinstance(ploidy, np.ndarray) and np.issubdtype(ploidy.dtype, np.integer):
            if ploidy.ndim != 1 or ploidy.size != n_samples:
                raise ValueError("ploidy array doesn't must match sample number")

            min_ploidy, max_ploidy = np.min(ploidy), np.max(ploidy)
            if min_ploidy != max_ploidy:
                ploidy_arr = np.asarray(ploidy, dtype=np.uint8)
        else:
            raise ValueError('ploidy must be either integer, or numpy array of integers')
        
        if min_ploidy < 0 or max_ploidy > 63:
            raise ValueError('ploidy values must be between 0 and 63')
        
        return min_ploidy, max_ploidy, ploidy_arr

    def _validate_layout1_data(self, vector[string] _alleles, uint32_t n_samples, uint32_t n_genos, bool phased):
        ''' validate layout 1 data
        '''
        if self.layout == 1 and _alleles.size() != 2:
            raise ValueError('layout 1 requires exactly two alleles')
        elif self.layout == 1 and n_genos != 3:
            raise ValueError('layout 1 requires 3 genotype probabilities per variant')
        elif self.layout == 1 and phased:
            raise ValueError('layout 1 cannot use phased data')


    def add_variant(self, varid, rsid, chrom, uint32_t pos, alleles, 
                    genotypes, ploidy=2, bool phased=False,
                    int bit_depth=8):
        ''' add a variant to the bgen file on disk

        Args:
            varid: variant ID
            rsid: reference SNP ID
            chrom: chromosome the variant is on
            pos: nucleotide position of the variant
            alleles: list of allele strings
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
        cdef vector[string] _alleles = [x.encode('utf8') for x in alleles]

        if bit_depth < 1 or bit_depth > 32:
            raise ValueError(f'bit_depth must be between 1 and 32: {bit_depth}')
        
        # sanatize genotypes
        cdef double[:, :] geno_c
        cdef uint32_t n_samples, n_genos
        genotypes, n_samples, n_genos = self._validate_genotypes(genotypes)
        geno_c = self._make_contiguous(genotypes)
    
        # validate ploidy levels
        cdef int32_t min_ploidy, max_ploidy
        cdef uint8_t[:] ploidy_arr
        min_ploidy, max_ploidy, ploidy_arr = self._validate_ploidy(ploidy, n_samples)
        
        self._validate_layout1_data(_alleles, n_samples, n_genos, phased)
        
        # encode the genotypes before writing anything, so that a variant with
        # invalid genotypes cannot leave a partial variant in the bgen file
        cdef uint32_t geno_len = n_samples * n_genos
        if min_ploidy == max_ploidy:
            self.thisptr.encode_genotype_data(_alleles.size(), &geno_c[0, 0],
                                           geno_len, min_ploidy, phased, bit_depth)
        else:
            self.thisptr.encode_genotype_data(_alleles.size(), &geno_c[0, 0],
                                           geno_len, &ploidy_arr[0], min_ploidy, 
                                           max_ploidy, phased, bit_depth)
        
        var_offset = self.thisptr.write_variant_header(_varid, _rsid, _chrom, pos, _alleles, n_samples)
        end_offset = self.thisptr.write_genotype_data()
        
        self.indexer.add_variant(chrom, int(pos), rsid, alleles, var_offset, 
                                 end_offset - var_offset)

    def add_variant_direct(self, variant):
        ''' insert a BgenVar directly into the bgen file
        '''
        if not self.is_open:
            raise ValueError("bgen file is closed")

        chrom = variant.chrom
        pos = int(variant.pos)
        rsid = variant.rsid
        alleles = variant.alleles
        cdef vector[uint8_t] data = variant.copy_data()
        var_offset = self.thisptr.write_variant_direct(data)
        end_offset = var_offset + len(data)

        self.indexer.add_variant(chrom, pos, rsid, alleles, var_offset, end_offset - var_offset)

    def __enter__(self):
        return self
    
    def __exit__(self, exc_type, exc_value, traceback):
        self.close()
        return False
    
    def close(self):
        if self.is_open:
            del self.thisptr
        if self.indexer is not None:
            self.indexer.close()
        self.is_open = False
        self.indexer = None
        if _IS_WIN32 and time is not None:
            time.sleep(0.01)
