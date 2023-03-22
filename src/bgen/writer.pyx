# cython: language_level=3, boundscheck=False, emit_linenums=True

import logging
from pathlib import Path

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
        void write_variant_header(string &varid, string &rsid, string &chrom, 
                        uint32_t &pos, vector[string] &alleles, uint32_t _n_samples) except +
        void add_genotype_data(uint16_t n_alleles,
                         double *genotypes, uint32_t geno_len, uint8_t ploidy,
                         bool phased, uint8_t bit_depth) except +
        void add_genotype_data(uint16_t n_alleles,
                         double *genotypes, uint32_t geno_len, vector[uint8_t] &ploidy,
                         uint32_t min_ploidy, uint32_t max_ploidy,
                         bool phased, uint8_t bit_depth) except +

cdef class BgenWriter:
    ''' class to open bgen files from disk, and access variant data within
    '''
    cdef CppBgenWriter * thisptr
    cdef string path
    cdef bool is_open
    def __cinit__(self, path, uint32_t n_samples, string free_data, 
                  compression='zstd', layout=2, vector[string] samples=[]):
        if isinstance(path, Path):
            path = str(path)
        
        if compression not in [None, 'zstd', 'zlib']:
            raise ValueError(f'compression type {compression} not one of zlib or zstd')
        
        cdef uint32_t compress_flag=0
        if compression == 'zlib':
            compress_flag = 1
        elif compression == 'zstd':
            compress_flag = 2

        self.path = path.encode('utf8')
        
        logging.debug(f'opening CppBgenWriter from {self.path.decode("utf")}')
        self.thisptr = new CppBgenWriter(self.path, n_samples, free_data, compress_flag, layout, samples)
        self.is_open = True
    
    def __dealloc__(self):
        if self.is_open:
          del self.thisptr
        self.is_open = False
    
    def __repr__(self):
        return f'BgenFile("{self.path.decode("utf8")}", "{self.sample_path.decode("utf8")}")'
    
    def add_variant(self, string varid, string rsid, string chrom, uint32_t pos, 
                    vector[string] alleles, uint32_t n_samples, 
                    double[:,:] genotypes, vector[uint8_t] ploidy=[], 
                    bool phased=False, uint8_t bit_depth=8):
        ''' get list of samples in the bgen file
        '''
        if not self.is_open:
            raise ValueError("bgen file is closed")
        self.thisptr.write_variant_header(varid, rsid, chrom, pos, alleles, n_samples)

        cdef uint32_t ploidy_n
        cdef uint32_t ploidy_min, ploidy_max
        cdef geno_len = genotypes.shape[0] * genotypes.shape[1]
        if ploidy.size() == 0:
            ploidy_n = 2
            if phased:
                ploidy_n = len(genotypes[0, ]) // alleles.size()
            self.thisptr.add_genotype_data(alleles.size(), &genotypes[0, 0], 
                                           geno_len, ploidy_n, phased, bit_depth)
        else:
            self.thisptr.add_genotype_data(alleles.size(), &genotypes[0, 0], 
                                           geno_len, ploidy,  min(ploidy), 
                                           max(ploidy), phased, bit_depth)

    def __enter__(self):
        return self
    
    def __exit__(self, exc_type, exc_value, traceback):
        self.close()
        return False
    
    def close(self):
        self.__dealloc__()
