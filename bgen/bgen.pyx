# cython: language_level=3, boundscheck=False

from libcpp cimport bool
from libcpp.string cimport string
from libcpp.vector cimport vector
from libc.stdint cimport uint32_t, uint64_t

from cython.operator cimport dereference as deref

import numpy as np

cdef extern from "<iostream>" namespace "std":
    cdef cppclass istream:
        istream& read(char *, int) except +
        istream& seekg(long) except +

cdef extern from "<iostream>" namespace "std::ios_base":
    cdef cppclass open_mode:
        pass
    cdef open_mode binary

cdef extern from "<fstream>" namespace "std":
    cdef cppclass ifstream(istream):
        ifstream() except +
        ifstream(const string&) except +
        ifstream(const string&, open_mode) except +

cdef extern from 'variant.h' namespace 'bgen':
    cdef cppclass Variant:
        # declare class constructor and methods
        Variant(ifstream & handle, uint64_t & offset, int layout, int compression, int expected_n) except +
        Variant() except +
        vector[float] minor_allele_dosage()
        
        # declare public attributes
        string varid, rsid, chrom, minor_allele
        int pos
        long offset
        vector[string] alleles

cdef extern from 'samples.h' namespace 'bgen':
    cdef cppclass Samples:
        Samples(ifstream & handle, int n_samples) except +
        Samples(string path, int n_samples) except +
        Samples(int n_samples) except +
        Samples() except +
        
        # declare public attributes
        vector[string] samples

cdef extern from 'header.h' namespace 'bgen':
    cdef cppclass Header:
        Header(ifstream & handle) except +
        Header();
        uint32_t offset
        uint32_t nvariants
        uint32_t nsamples
        int compression
        int layout
        bool has_sample_ids

cdef extern from 'bgen.h' namespace 'bgen':
    cdef cppclass Bgen:
        # declare class constructor and methods
        Bgen(string path, string sample_path) except +
        Variant & operator[](int idx)
        Variant & get(int idx)
        void drop_variants(vector[int] indices)
        vector[string] varids()
        vector[string] rsids()
        vector[string] chroms()
        vector[uint32_t] positions()
        
        # declare public attributes
        vector[Variant] variants
        Samples samples
        Header header

cdef class IFStream:
    ''' basic cython implementation of std::ifstream, for easy pickling
    '''
    cdef ifstream * ptr
    cdef string path
    def __cinit__(self, string path):
        self.path = path
        self.ptr = new ifstream(path)
    
    def seekg(self, long offset):
        self.ptr.seekg(offset)
    
    def __str__(self):
        return self.path.decode('utf8')
    
    def __dealloc__(self):
        # make sure to clean up ifstream memory once finished
        del self.ptr
    
    def __reduce__(self):
        return (self.__class__, (self.path, ))

cdef class BgenHeader:
    cdef uint32_t offset
    cdef uint32_t nvariants
    cdef uint32_t nsamples
    cdef int compression
    cdef int layout
    cdef bool has_sample_ids
    def __cinit__(self, uint32_t offset, uint32_t nvariants, uint32_t nsamples,
            int compression, int layout, bool has_sample_ids):
        self.offset = offset
        self.nvariants = nvariants
        self.nsamples = nsamples
        self.compression = compression
        self.layout = layout
        self.has_sample_ids = has_sample_ids
    
    def __repr__(self):
        return f'BgenHeader(offset={self.offset}, nvariants={self.nvariants}, ' \
            f'nsamples={self.nsamples}, compression={self.compression}, ' \
            f'layout={self.layout}, has_sample_ids={self.has_sample_ids})'

cdef class BgenVar:
    cdef Variant thisptr
    cdef IFStream handle
    cdef uint64_t offset
    cdef int layout, compression, expected_n
    def __cinit__(self, IFStream handle, uint64_t offset, int layout, int compression, int expected_n):
        self.handle = handle
        self.offset = offset
        self.layout = layout
        self.compression = compression
        self.expected_n = expected_n
        
        self.thisptr = Variant(deref(self.handle.ptr), offset, layout, compression, expected_n)
    
    def __repr__(self):
       return f'BgenVar("{self.varid}", "{self.rsid}", "{self.chrom}", {self.pos}, {self.alleles})'
    
    def __str__(self):
       return f'{self.rsid} - {self.chrom}:{self.pos} {self.alleles}'
    
    def __reduce__(self):
        ''' enable pickling of a BgenVar object
        '''
        return (self.__class__, (self.handle, self.thisptr.offset, self.layout, self.compression, self.expected_n))
    
    @property
    def varid(self):
      return self.thisptr.varid.decode('utf8')
    @property
    def rsid(self):
        return self.thisptr.rsid.decode('utf8')
    @property
    def chrom(self):
        return self.thisptr.chrom.decode('utf8')
    @property
    def pos(self):
        return self.thisptr.pos
    @property
    def alleles(self):
        return [x.decode('utf8') for x in self.thisptr.alleles]
    @property
    def minor_allele(self):
        ''' get the minor allele of a biallelic variant
        '''
        return self.thisptr.minor_allele.decode('utf8')
    @property
    def minor_allele_dosage(self):
        ''' get the dosage for the minor allele for a biallelic variant
        
        In order for this to be fast, we need to get the dosage data as a vector,
        then get a memory view on that data, and finally return as a numpy array
        '''
        # get the vector data for the dosage
        cdef vector[float] dosage = self.thisptr.minor_allele_dosage()
        # convert to a memory view, see:
        # https://cython.readthedocs.io/en/latest/src/userguide/memoryviews.html#coercion-to-numpy
        cdef float[::1] arr = <float [:dosage.size()]>dosage.data()
        # return as a numpy array
        return np.asarray(arr)

cdef class BgenFile:
    cdef Bgen * thisptr
    cdef string path, sample_path
    cdef IFStream handle
    def __cinit__(self, str path, str sample_path=''):
        self.path = path.encode('utf8')
        self.sample_path = sample_path.encode('utf8')
        self.thisptr = new Bgen(self.path, self.sample_path)
        self.handle = IFStream(self.path)
    
    def __dealloc__(self):
        del self.thisptr
    
    def __repr__(self):
        return f'BgenFile("{self.path.decode("utf8")}", "{self.sample_path.decode("utf8")}")'
    
    def __iter__(self):
        for idx in range(len(self)):
            yield self[idx]
    
    def __len__(self):
        return self.thisptr.variants.size()
    
    def __getitem__(self, int idx):
        ''' pull out a Variant by index position
        '''
        if idx >= len(self) or idx < 0:
            raise IndexError(f'cannot get Variant at index: {idx}')
        
        cdef long offset = self.thisptr.variants[idx].offset
        return BgenVar(self.handle, offset, self.thisptr.header.layout,
          self.thisptr.header.compression, self.thisptr.header.nsamples)
    
    @property
    def header(self):
      ''' get header info from bgen file
      '''
      hdr = self.thisptr.header
      return BgenHeader(hdr.offset, hdr.nvariants, hdr.nsamples,
          hdr.compression, hdr.compression, hdr.has_sample_ids)
    
    @property
    def samples(self):
      ''' get list of samples in the bgen file
      '''
      samples = self.thisptr.samples.samples
      return [x.decode('utf8') for x in samples]
    
    def drop_variants(self, list indices):
        ''' drops variants from bgen by indices, for avoiding processing variants
        '''
        self.thisptr.drop_variants(indices)
    
    def varids(self):
      ''' get the varint IDs of all variants in the bgen file
      '''
      varids = self.thisptr.varids()
      return [x.decode('utf8') for x in varids]
    
    def rsids(self):
      ''' get the rsIDs of all variants in the bgen file
      '''
      rsids = self.thisptr.rsids()
      return [x.decode('utf8') for x in rsids]
    
    def chroms(self):
        ''' get the chromosomes of all variants in the bgen file
        '''
        chroms = self.thisptr.chroms()
        return [x.decode('utf8') for x in chroms]
    
    def positions(self):
        ''' get the positions of all variants in the bgen file
        '''
        return self.thisptr.positions()
  
