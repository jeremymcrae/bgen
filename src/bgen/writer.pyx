# cython: language_level=3, boundscheck=False, emit_linenums=True

import logging
from pathlib import Path
import sqlite3
import sys
import time

from libcpp cimport bool
from libcpp.string cimport string
from libcpp.vector cimport vector
from libc.stdint cimport uint8_t, uint16_t, uint32_t, uint64_t

from cython.operator cimport dereference as deref

import numpy as np

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

cdef extern from 'genotypes.h' namespace 'bgen':
    uint32_t get_max_probs(int &max_ploidy, int &n_alleles, bool &phased)

cdef extern from 'writer.h' namespace 'bgen':
    cdef cppclass CppBgenWriter:
        # declare class constructor and methods
        CppBgenWriter(string &path, uint32_t n_samples, string &free_data, 
                   uint32_t compression, uint32_t layout, vector[string] &samples) except +
        uint64_t write_variant_header(string &varid, string &rsid, string &chrom, 
                        uint32_t &pos, vector[string] &alleles, uint32_t _n_samples) except +
        uint64_t add_genotype_data(uint16_t n_alleles,
                         double *genotypes, uint32_t geno_len, uint8_t ploidy,
                         bool phased, uint8_t bit_depth) except +
        uint64_t add_genotype_data(uint16_t n_alleles,
                         double *genotypes, uint32_t geno_len, uint8_t *ploidy,
                         uint32_t min_ploidy, uint32_t max_ploidy,
                         bool phased, uint8_t bit_depth) except +

class Indexer:
    ''' class to automatically index bgen files as they are being constructed
    '''
    def __init__(self, bgen_path):
        index_path = Path(str(bgen_path + '.bgi'))
        if index_path.exists():
            index_path.unlink()
        self.conn = sqlite3.connect(index_path)
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
        query = '''INSERT INTO Variant VALUES (?, ?, ?, ?, ?, ?, ?, ?)'''
        params = (chrom, pos, rsid, len(alleles), alleles[0], alleles[1], offset, size)
        self.cur.execute(query, params)
    
    def close(self):
        self.conn.commit()
        if sys.platform == 'win32':
            time.sleep(0.01)

cdef class BgenWriter:
    ''' class to open bgen files from disk, and access variant data within
    '''
    cdef CppBgenWriter * thisptr
    cdef string path
    cdef bool is_open
    cdef object indexer
    def __cinit__(self, path, uint32_t n_samples, samples=[], compression='zstd',
                  layout=2, metadata=None):
        if isinstance(path, Path):
            path = str(path)
        
        if compression not in [None, 'zstd', 'zlib']:
            raise ValueError(f'compression type {compression} not one of zlib or zstd')
        
        cdef uint32_t compress_flag=0
        if compression == 'zlib':
            compress_flag = 1
        elif compression == 'zstd':
            compress_flag = 2
        
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
    
    def add_variant(self, varid, rsid, chrom, uint32_t pos, alleles, 
                    double[:,:] genotypes, ploidy=2, bool phased=False,
                    uint8_t bit_depth=8):
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
            bit_depth: interger from 1-32 (inclusive) for how many bits to store
                each genotype in.
        '''

        # re-define variables into cpp objects
        cdef string _varid = varid.encode('utf8')
        cdef string _rsid = rsid.encode('utf8')
        cdef string _chrom = chrom.encode('utf8')
        cdef vector[string] _alleles = [x.encode('utf8') for x in alleles]
        cdef uint32_t n_samples = len(genotypes)

        if not self.is_open:
            raise ValueError("bgen file is closed")
        var_offset = self.thisptr.write_variant_header(_varid, _rsid, _chrom, pos, _alleles, n_samples)

        # determine ploidy levels
        cdef uint32_t ploidy_n=0
        cdef uint8_t[:] ploidy_arr = np.array([], dtype=np.uint8)
        if isinstance(ploidy, int):
            ploidy_n = ploidy
        elif isinstance(ploidy, np.ndarray):
            ploidy_arr = ploidy
        else:
            raise ValueError('ploidy must be either integer, or numpy array of integers')
        
        cdef geno_len = genotypes.shape[0] * genotypes.shape[1]
        if len(ploidy_arr) == 0:
            end_offset = self.thisptr.add_genotype_data(_alleles.size(), &genotypes[0, 0], 
                                           geno_len, ploidy_n, phased, bit_depth)
        else:
            end_offset = self.thisptr.add_genotype_data(_alleles.size(), &genotypes[0, 0], 
                                           geno_len, &ploidy_arr[0], min(ploidy), 
                                           max(ploidy), phased, bit_depth)
        
        self.indexer.add_variant(chrom, int(pos), rsid, alleles, var_offset, 
                                 end_offset - var_offset)

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
        if sys.platform == 'win32':
            time.sleep(0.01)
