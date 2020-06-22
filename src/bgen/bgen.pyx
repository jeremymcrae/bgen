# cython: language_level=3, boundscheck=False, emit_linenums=True

import logging
from pathlib import Path

from libcpp cimport bool
from libcpp.string cimport string
from libcpp.vector cimport vector
from libc.stdint cimport uint8_t, uint32_t, uint64_t

from cython.operator cimport dereference as deref

import numpy as np

from bgen.index import Index

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
        float * minor_allele_dosage()
        float * probs_1d()
        int probs_per_sample()
        bool phased()
        uint8_t * ploidy()
        
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
        Bgen(string path, string sample_path, bool delay_parsing) except +
        void parse_all_variants()
        Variant & next_var() except +
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
    ''' holds information about the Bgen file, obtained from the intial header.
    '''
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
    ''' holds data for a Variant from a bgen file
    
    This constructs a new Variant, rather than using a object pointer, in order
    to make pickling the object easier.
    
    Initialization takes about 1e-5 seconds per variant, so we can only run
    through variants in a file at 100,000 variants per second at most (assuming
    no other work is being done).
    
    This shouldn't be a limitation in practise, since the rate limiting part is
    parsing genotype information, which runs at about 90 variants per second
    for files with 500,000 samples.
    '''
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
        
        # construct new Variant from the handle, offset and other file info
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
    def is_phased(self):
      return self.thisptr.phased()
    @property
    def ploidy(self):
        ''' get the ploidy for each sample
        '''
        cdef uint8_t * ploid = self.thisptr.ploidy()
        return np.copy(np.asarray(<uint8_t [:self.expected_n]>ploid))
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
        # get the vector data for the dosage. We use a cython memoryview to
        # quickly convert to a numpy array
        # https://cython.readthedocs.io/en/latest/src/userguide/memoryviews.html#coercion-to-numpy
        cdef float * dosage = self.thisptr.minor_allele_dosage()
        return np.asarray(<float [:self.expected_n]>dosage)
    @property
    def probabilities(self):
        ''' get the allelic probabilities for a variant
        '''
        cdef float * probs = self.thisptr.probs_1d()
        cdef int cols = self.thisptr.probs_per_sample()
        cdef int size = self.expected_n * cols
        if self.is_phased:
            ploidy = self.ploidy
            size = ploidy.sum() * cols
        
        arr = np.asarray(<float [:size]>probs)
        data = np.reshape(arr, (-1, cols))
        
        cdef int current = 0
        cdef int phase_width = data.shape[1]
        if self.is_phased:
            # phased data initially comes as one row per haploytpe. This is
            # reshaped to concatenate haplotype data into single row. Fill in a
            # new array from the old data row by row
            
            # create an empty array filled with nans
            ragged = np.empty((len(ploidy), ploidy.max() * cols))
            ragged.fill(np.nan)
            
            # fill in the empty array
            for i, x in enumerate(ploidy):
                for y in range(x):
                    start = y * phase_width
                    end = start + phase_width
                    ragged[i, start:end] = data[current]
                    current += 1
            
            data = ragged
        
        return data

cdef class BgenFile:
    ''' class to open bgen files from disk, and access variant data within
    '''
    cdef Bgen * thisptr
    cdef string path, sample_path
    cdef bool delay_parsing
    cdef IFStream handle
    cdef object index
    cdef bool is_open
    def __cinit__(self, path, sample_path='', bool delay_parsing=False):
        if isinstance(path, Path):
            path = str(path)
        if isinstance(sample_path, Path):
            sample_path = str(sample_path)
        
        delay_parsing |= self._check_for_index(path)
        
        self.path = path.encode('utf8')
        self.sample_path = sample_path.encode('utf8')
        self.delay_parsing = delay_parsing
        
        samp = '' if sample_path == '' else f', (samples={self.sample_path.decode("utf")})'
        logging.debug(f'opening BgenFile from {self.path.decode("utf")}{samp}')
        self.thisptr = new Bgen(self.path, self.sample_path, self.delay_parsing)
        self.handle = IFStream(self.path)
        self.is_open = True
    
    def __dealloc__(self):
        if self.is_open:
          del self.thisptr
          self.handle = None
          self.index = None
        
        self.is_open = False
    
    def __repr__(self):
        return f'BgenFile("{self.path.decode("utf8")}", "{self.sample_path.decode("utf8")}")'
    
    def __iter__(self):
        return self
    
    def __next__(self):
        ''' iterate through all variants in the bgen file
        '''
        # while True:
        try:
            offset = self.thisptr.next_var().offset
            return BgenVar(self.handle, offset, self.thisptr.header.layout,
                self.thisptr.header.compression, self.thisptr.header.nsamples)
        except IndexError:
            raise StopIteration
    
    def __len__(self):
      if not self.is_open:
          raise ValueError("bgen file is closed")
      
      length = self.thisptr.variants.size()
      if length > 0:
          return length
      else:
          return self.thisptr.header.nvariants
    
    def __getitem__(self, int idx):
        ''' pull out a Variant by index position
        '''
        if idx >= len(self) or idx < 0:
            raise IndexError(f'cannot get Variant at index: {idx}')
        
        # account for lazy loading variants from bgen
        if self.index is None and self.thisptr.variants.size() == 0:
            self.thisptr.parse_all_variants()
        
        cdef long offset
        offset = self.index.offset_by_index(idx) if self.index else self.thisptr.variants[idx].offset
        return BgenVar(self.handle, offset, self.thisptr.header.layout,
          self.thisptr.header.compression, self.thisptr.header.nsamples)
    
    def _check_for_index(self, bgen_path):
        ''' creates self.index if a bgenix index file is available
        '''
        index_path = Path(bgen_path + '.bgi')
        idx_exists = index_path.exists()
        self.index = Index(index_path) if idx_exists else None
        return idx_exists
    
    @property
    def header(self):
      ''' get header info from bgen file
      '''
      if not self.is_open:
          raise ValueError("bgen file is closed")
      
      hdr = self.thisptr.header
      return BgenHeader(hdr.offset, hdr.nvariants, hdr.nsamples,
          hdr.compression, hdr.compression, hdr.has_sample_ids)
    
    @property
    def samples(self):
      ''' get list of samples in the bgen file
      '''
      if not self.is_open:
          raise ValueError("bgen file is closed")
      
      samples = self.thisptr.samples.samples
      return [x.decode('utf8') for x in samples]
    
    def drop_variants(self, list indices):
        ''' drops variants from bgen by indices, for avoiding processing variants
        '''
        if not self.is_open:
            raise ValueError("bgen file is closed")
        
        if self.delay_parsing:
            self.thisptr.parse_all_variants()
        
        self.thisptr.drop_variants(indices)
    
    def fetch(self, chrom, start=None, stop=None):
        ''' fetches all variants within a genomic region
        
        Args:
            chrom: chromosome that variants must be on
            start: start nucleotide of region. If None, gets all variants on chromosome
            stop: end nucleotide of region. If None, gets variants with positions after start
        
        Yields:
            BgenVars for variants within the genome region
        '''
        if not self.is_open:
            raise ValueError('bgen file is closed')
        
        if not self.index:
            raise ValueError("can't fetch variants without index")
        
        for offset in self.index.fetch(chrom, start, stop):
            yield BgenVar(self.handle, offset, self.thisptr.header.layout,
                self.thisptr.header.compression, self.thisptr.header.nsamples)
    
    def with_rsid(self, rsid):
      ''' get BgenVar from file given an rsID
      '''
      if not self.is_open:
          raise ValueError('bgen file is closed')
      
      if self.index:
          offset = self.index.offset_by_rsid(rsid)
          return BgenVar(self.handle, offset, self.thisptr.header.layout,
              self.thisptr.header.compression, self.thisptr.header.nsamples)
      
      if not self.delay_parsing:
          idx = [i for i, x in enumerate(self.rsids) if x == rsid]
          if len(idx) == 0:
              raise ValueError(f'cannot find variant match for {rsid}')
          elif len(idx) > 1:
              raise ValueError(f'multiple variant matches for {rsid}')
          return self[idx]
      
      raise ValueError("can't get variant without fully loading the bgen, or indexing")
    
    def at_position(self, pos):
      ''' get BgenVar from file given a position
      '''
      if not self.is_open:
          raise ValueError('bgen file is closed')
      
      if self.index:
          offset = self.index.offset_by_pos(pos)
          return BgenVar(self.handle, offset, self.thisptr.header.layout,
              self.thisptr.header.compression, self.thisptr.header.nsamples)
      
      if not self.delay_parsing:
          idx = [i for i, x in enumerate(self.positions) if x == pos]
          if len(idx) == 0:
              raise ValueError(f'cannot find variant match at pos: {pos}')
          elif len(idx) > 1:
              raise ValueError(f'multiple variant matches at pos: {pos}')
          return self[idx]
      
      raise ValueError("can't get variant without fully loading the bgen, or indexing")
    
    def varids(self):
      ''' get the variant IDs of all variants in the bgen file
      '''
      if not self.is_open:
          raise ValueError("bgen file is closed")
      
      if self.index:
          raise ValueError("can't load varids when using an index file")
      
      varids = self.thisptr.varids()
      return [x.decode('utf8') for x in varids]
    
    def rsids(self):
      ''' get the rsIDs of all variants in the bgen file
      '''
      if not self.is_open:
          raise ValueError("bgen file is closed")
      
      if self.index:
          return self.index.rsids
      
      rsids = self.thisptr.rsids()
      return [x.decode('utf8') for x in rsids]
    
    def chroms(self):
        ''' get the chromosomes of all variants in the bgen file
        '''
        if not self.is_open:
            raise ValueError("bgen file is closed")
        
        if self.index:
            return self.index.chroms
        
        chroms = self.thisptr.chroms()
        return [x.decode('utf8') for x in chroms]
    
    def positions(self):
        ''' get the positions of all variants in the bgen file
        '''
        if not self.is_open:
            raise ValueError("bgen file is closed")
        
        if self.index:
            return self.index.positions
        
        return self.thisptr.positions()
    
    def __enter__(self):
        return self
    
    def __exit__(self, exc_type, exc_value, traceback):
        self._close()
        return False
    
    def _close(self):
        if self.is_open:
            del self.thisptr
            self.handle = None
            self.index = None
        
        self.is_open = False
