# cython: language_level=3, boundscheck=False, emit_linenums=True

import logging
from pathlib import Path
import sys
import warnings

from libcpp cimport bool
from libcpp.memory cimport shared_ptr, make_shared
from libcpp.string cimport string
from libcpp.vector cimport vector
from libc.stdint cimport uint8_t, uint32_t, uint64_t, uintptr_t
from libc.string cimport memcpy

from cython.operator cimport dereference as deref

import numpy as np

from bgen.index import Index

cdef extern from "<iostream>" namespace "std::ios_base":
    cdef cppclass open_mode:
        pass
    cdef open_mode binary
    cdef cppclass iostate:
        pass
    cdef iostate failbit

cdef extern from "<iostream>" namespace "std":
    cdef cppclass istream:
        istream() except +
        void setstate(iostate state) except +

cdef extern from 'variant.h' namespace 'bgen':
    cdef cppclass Variant:
        # declare class constructor and methods
        Variant(istream * handle, uint64_t & offset, int layout, int compression, int expected_n, bool is_stdin) except +
        Variant() except +
        void minor_allele_dosage(float * dosage) except +
        void alt_dosage(float * dosage) except +
        void probs_1d(float * dosage) except +
        int probs_per_sample() except +
        bool phased() except +
        uint8_t * ploidy() except +
        vector[uint8_t] copy_data() except +
        
        # declare public attributes
        string varid, rsid, chrom, minor_allele
        int pos
        long offset
        uint64_t next_variant_offset
        vector[string] alleles
        istream * handle

cdef extern from 'samples.h' namespace 'bgen':
    cdef cppclass Samples:
        Samples(istream * handle, int n_samples) except +
        Samples(string path, int n_samples) except +
        Samples(int n_samples) except +
        Samples() except +
        
        # declare public attributes
        vector[string] samples

cdef extern from 'header.h' namespace 'bgen':
    cdef cppclass Header:
        Header(istream * handle) except +
        Header();
        uint32_t offset
        uint32_t nvariants
        uint32_t nsamples
        int compression
        int layout
        string extra
        bool has_sample_ids

cdef extern from 'reader.h' namespace 'bgen':
    cdef cppclass CppBgenReader:
        # declare class constructor and methods
        CppBgenReader(string path, string sample_path, bool delay_parsing) except +
        void parse_all_variants()
        Variant & operator[](int idx)
        Variant & get(int idx)
        void drop_variants(vector[int] indices)
        vector[string] varids()
        vector[string] rsids()
        vector[string] chroms()
        vector[uint32_t] positions()
        
        # declare public attributes
        istream * handle
        vector[Variant] variants
        Samples samples
        Header header
        uint64_t offset

cdef extern from 'utils.h' namespace 'bgen':
    cdef struct Range:
        uint8_t _min
        uint8_t _max
    uint64_t fast_ploidy_sum(uint8_t * ploidy, uint32_t & size) except +
    Range fast_range(uint8_t * ploidy, uint32_t & size)

cdef class IStream:
    ''' basic cython implementation of std::istream, for easy pickling
    '''
    cdef istream * ptr
    def __cinit__(self, uint64_t ptr):
        self.ptr = <istream*>ptr
    
    def __str__(self):
        return f'std::istream at {<uint64_t>self.ptr}'
    
    def __dealloc__(self):
        pass
    
    def __reduce__(self):
        return (self.__class__, (<uint64_t>self.ptr, ))

cdef class OpenStatus:
    ''' class to share status of whether a bgen file is currently open
    
    This uses a shared pointer so we can store the status in a BgenReader
    object, but also pass that status to BgenVariant objects. When a bgen
    file is closed, the status changes for all loaded variants, in order
    to prevent reading more data from the closed file.
    
    This may not work for pickled BgenVariant objects (after unpickling), as
    the BgenReader object may have been closed after the BgenVariant was 
    pickled. In order to allow pickled BgenVariants, as part of unpickling 
    the BgenVariant, it creates a new OpenStatus object where status=True.
    '''
    cdef shared_ptr[bool] ptr
    def __cinit__(self):
        self.ptr = make_shared[bool](True)
    def __str__(self):
        return f'status={deref(self.ptr)}'
    def __eq__(self, other):
        return deref(self.ptr) == other
    def off(self):
        self.ptr = make_shared[bool](False)
    def __reduce__(self):
        return (self.__class__, ())

cdef class BgenHeader:
    ''' holds information about the Bgen file, obtained from the intial header.
    '''
    cdef uint32_t _offset
    cdef uint32_t _nvariants
    cdef uint32_t _nsamples
    cdef int _compression
    cdef int _layout
    cdef bool _has_sample_ids
    cdef string _metadata
    cdef object compress_formats
    def __cinit__(self, uint32_t offset, uint32_t nvariants, uint32_t nsamples,
            int compression, int layout, bool has_sample_ids, string metadata):
        self._offset = offset
        self._nvariants = nvariants
        self._nsamples = nsamples
        self._compression = compression
        self._layout = layout
        self._has_sample_ids = has_sample_ids
        self._metadata = metadata
        self.compress_formats = {0: None, 1: 'zlib', 2: 'zstd'}
    
    def __repr__(self):
        return f'BgenHeader(offset={self.offset}, nvariants={self.nvariants}, ' \
            f'nsamples={self.nsamples}, compression={self.compression}, ' \
            f'layout={self.layout}, has_sample_ids={self.has_sample_ids})'
    
    @property
    def offset(self): return self._offset
    @property
    def nsamples(self): return self._nsamples
    @property
    def nvariants(self): return self._nvariants
    @property
    def compression(self): return [None, 'zlib', 'zstd'][self._compression]
    @property
    def compression(self): return self.compress_formats[self._compression]
    @property
    def layout(self): return self._layout
    @property
    def has_sample_ids(self): return self._has_sample_ids
    @property
    def metadata(self): return self._metadata.decode('utf8')

cdef class BgenVar:
    ''' holds data for a Variant from a bgen file
    
    This constructs a new Variant, rather than using a object pointer, in order
    to make pickling the object easier.
    
    Initialization takes about 1e-5 seconds per variant, so we can only run
    through variants in a file at 100,000 variants per second at most (assuming
    no other work is being done).
    
    This shouldn't be a limitation in practise, since the rate limiting part is
    parsing genotype information, which runs at about 500 variants per second
    for files with 500,000 samples.
    '''
    cdef Variant * thisptr
    cdef IStream handle
    cdef uint64_t offset
    cdef int layout, compression, expected_n
    cdef bool is_stdin
    cdef OpenStatus is_open
    def __cinit__(self,
                  IStream handle,
                  uint64_t offset,
                  int layout,
                  int compression,
                  int expected_n,
                  bool is_stdin,
                  OpenStatus is_open,
                  ):
        self.handle = handle
        self.offset = offset
        self.layout = layout
        self.compression = compression
        self.expected_n = expected_n
        self.is_stdin = is_stdin
        self.is_open = is_open
        
        # construct new Variant from the handle, offset and other file info
        self.thisptr = new Variant(self.handle.ptr, offset, layout, compression, expected_n, is_stdin)
    
    def __repr__(self):
       return f'BgenVar("{self.varid}", "{self.rsid}", "{self.chrom}", {self.pos}, {self.alleles})'
    
    def __str__(self):
       return f'{self.rsid} - {self.chrom}:{self.pos} {self.alleles}'
    
    def __reduce__(self):
        ''' enable pickling of a BgenVar object
        '''
        warnings.warn("pickling BgenVar - make sure their BgenReader objects exist when unpickling", RuntimeWarning)
        return (self.__class__, (self.handle, self.thisptr.offset, self.layout, self.compression, self.expected_n, self.is_stdin, self.is_open))
    
    def __dealloc__(self):
        del self.thisptr
    
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
    def fileoffset(self):
        return self.offset
    @property
    def next_variant_offset(self):
        return self.thisptr.next_variant_offset
    @property
    def alleles(self):
        return [x.decode('utf8') for x in self.thisptr.alleles]
    cdef __check_closed(self):
        ''' make sure we cannot read from a closed bgen file
        
        If the bgen file is closed before we try to access the genotype data,
        then we should not be able to read from the file. However, this
        object holds a pointer to a istream object, and it does not know if
        the underlying memory is valid, so could try to access it anyway.
        
        We prevent this by checking a shared pointer held by all the BgenVars 
        which were opened by a given BgenReader. This shared pointer will have
        been set to false when the bgen file closed, so we can check that
        shared status, and set the istream iostate to fail. Later cpp calls which
        would read from the stale file handle will raise ValuErrors instead due
        to checking the failbit before attempting to read from the istream.
        
        However, this does not work for pickled BgenVars, since we cannot know
        whether the istream object exists or not.
        '''
        if not deref(self.is_open.ptr):
            self.thisptr.handle.setstate(failbit)
    @property
    def is_phased(self):
        return self.thisptr.phased()
    @property
    def ploidy(self):
        ''' get the ploidy for each sample
        '''
        self.__check_closed()
        cdef uint8_t * ploid = self.thisptr.ploidy()
        cdef uint64_t size = self.expected_n
        cdef uint8_t[::1] arr = np.empty(size, dtype=np.uint8, order='C')
        memcpy(&arr[0], ploid, size)
        return arr
    @property
    def minor_allele(self):
        ''' get the minor allele of a biallelic variant
        '''
        return self.thisptr.minor_allele.decode('utf8')
    @property
    def minor_allele_dosage(self):
        ''' dosage for the minor allele for a biallelic variant
        '''
        self.__check_closed()
        cdef float[:] dose = np.empty(self.expected_n, dtype=np.float32, order='C')
        self.thisptr.minor_allele_dosage(&dose[0])
        return np.asarray(dose)
    @property
    def alt_dosage(self):
        ''' dosage for the alt allele for a biallelic variant
        '''
        self.__check_closed()
        cdef float[:] dose = np.empty(self.expected_n, dtype=np.float32, order='C')
        self.thisptr.alt_dosage(&dose[0])
        return np.asarray(dose)
    @property
    def probabilities(self):
        ''' get the allelic probabilities for a variant
        '''
        self.__check_closed()
        cdef int cols = self.thisptr.probs_per_sample()
        cdef uint32_t n_samples = self.expected_n
        cdef uint64_t size = n_samples * cols
        cdef uint8_t[::1] ploidy
        if self.is_phased:
            ploidy = self.ploidy
            size = fast_ploidy_sum(&ploidy[0], n_samples) * cols
        
        cdef float[:] arr = np.empty(size, dtype=np.float32, order='C')
        self.thisptr.probs_1d(&arr[0])
        
        cdef int current = 0
        cdef int phase_width
        cdef Range ploidy_range
        if self.is_phased:
            ploidy_range = fast_range(&ploidy[0], n_samples)
            max_ploidy = ploidy_range._max
            min_ploidy = ploidy_range._min
            if min_ploidy == max_ploidy:
                # quickly reshape probs if ploidy is constant
                width = max_ploidy * cols
                data = np.reshape(arr, (-1, width))
            else:
                # phased data initially comes as one row per haploytpe. This is
                # reshaped to concatenate haplotype data into single row. Fill
                # in a new array from the old data row by row
                data = np.reshape(arr, (-1, cols))
                phase_width = data.shape[1]
                
                # create an empty array filled with nans
                ragged = np.empty((len(ploidy), max_ploidy * cols))
                ragged.fill(np.nan)
                
                # fill in the empty array
                for i, x in enumerate(ploidy):
                    for y in range(x):
                        start = y * phase_width
                        end = start + phase_width
                        ragged[i, start:end] = data[current]
                        current += 1
                
                data = ragged
        else:
            data = np.reshape(arr, (-1, cols))
        
        return data
    
    def copy_data(self):
        ''' get a copy of the data on disk for the variant.
        
        This can be used to quickly copy data from one bgen into another, if
        you have the same set of samples between the source and destination
        bgens. This primarily avoids decompressing, decoding, re-encoding, and
        compressing the genotype data.
        '''
        self.__check_closed()
        cdef vector[uint8_t] data = self.thisptr.copy_data()
        return data

cdef class BgenReader:
    ''' class to open bgen files from disk, and access variant data within
    '''
    cdef CppBgenReader * thisptr
    cdef string path, sample_path
    cdef bool delay_parsing, is_stdin
    cdef IStream handle
    cdef object index
    cdef OpenStatus is_open
    cdef uint64_t offset
    def __cinit__(self, path, sample_path='', bool delay_parsing=False):
        if isinstance(path, Path):
            path = str(path)
        if isinstance(sample_path, Path):
            sample_path = str(sample_path)
        self.is_stdin = self.__is_from_stdin(path)
        if self.is_stdin:
            delay_parsing = True
            path = '/dev/stdin'
        
        delay_parsing |= self._check_for_index(path)
        
        self.path = path.encode('utf8')
        self.sample_path = sample_path.encode('utf8')
        self.delay_parsing = delay_parsing
        
        samp = '' if sample_path == '' else f', (samples={self.sample_path.decode("utf")})'
        logging.debug(f'opening BgenFile from {self.path.decode("utf")}{samp}')
        self.thisptr = new CppBgenReader(self.path, self.sample_path, self.delay_parsing)
        self.handle = IStream(<uint64_t>self.thisptr.handle)
        self.is_open = OpenStatus()
        self.offset = self.thisptr.offset
    
    def __is_from_stdin(self, bgen_path):
        if bgen_path is sys.stdin:
            return True
        elif str(bgen_path) == '/dev/stdin':
            return True
        elif str(bgen_path) == '-':
            return True
        return False
    
    def __dealloc__(self):
        self.close()
    
    def __repr__(self):
        return f'BgenFile("{self.path.decode("utf8")}", "{self.sample_path.decode("utf8")}")'
    
    def __iter__(self):
        return self
    
    def __next__(self):
        ''' iterate through all variants in the bgen file
        '''
        if not self.is_open == True:
            raise ValueError('bgen file is closed')
        
        try:
            var = BgenVar(self.handle, self.offset, self.thisptr.header.layout,
                self.thisptr.header.compression, self.thisptr.header.nsamples, self.is_stdin,
                self.is_open)
            self.offset = var.next_variant_offset
            return var
        except IndexError:
            raise StopIteration
    
    def __len__(self):
      if not self.is_open == True:
          raise ValueError("bgen file is closed")
      
      length = self.thisptr.variants.size()
      if length > 0:
          return length
      else:
          return self.thisptr.header.nvariants
    
    def __getitem__(self, int idx):
        ''' pull out a Variant by index position
        '''
        if not self.is_open == True:
            raise ValueError('bgen file is closed')
        
        if idx >= len(self) or idx < 0:
            raise IndexError(f'cannot get Variant at index: {idx}')
        
        # account for lazy loading variants from bgen
        if self.index is None and self.thisptr.variants.size() == 0:
            self.thisptr.parse_all_variants()
        
        cdef long offset
        offset = self.index.offset_by_index(idx) if self.index else self.thisptr.variants[idx].offset
        return BgenVar(self.handle, offset, self.thisptr.header.layout,
          self.thisptr.header.compression, self.thisptr.header.nsamples,
          self.is_stdin, self.is_open)
    
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
      if not self.is_open == True:
          raise ValueError("bgen file is closed")
      
      hdr = self.thisptr.header
      return BgenHeader(hdr.offset, hdr.nvariants, hdr.nsamples,
          hdr.compression, hdr.layout, hdr.has_sample_ids, hdr.extra)
    
    @property
    def samples(self):
      ''' get list of samples in the bgen file
      '''
      if not self.is_open == True:
          raise ValueError("bgen file is closed")
      
      samples = self.thisptr.samples.samples
      return [x.decode('utf8') for x in samples]
    
    def drop_variants(self, list indices):
        ''' drops variants from bgen by indices, for a a =ing processing variants
        '''
        if not self.is_open == True:
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
        if not self.is_open == True:
            raise ValueError('bgen file is closed')
        
        if not self.index:
            raise ValueError("can't fetch variants without index")
        
        for offset in self.index.fetch(chrom, start, stop):
            yield BgenVar(self.handle, offset, self.thisptr.header.layout,
                self.thisptr.header.compression, self.thisptr.header.nsamples,
                self.is_stdin, self.is_open)
    
    def with_rsid(self, rsid):
      ''' get BgenVar from file given an rsID
      '''
      if not self.is_open == True:
          raise ValueError('bgen file is closed')
      
      if self.index:
          offsets = self.index.offset_by_rsid(rsid)
          return [BgenVar(self.handle, int(offset), self.thisptr.header.layout,
                          self.thisptr.header.compression, self.thisptr.header.nsamples,
                          self.is_stdin, self.is_open)
                  for offset in offsets]
      
      if not self.delay_parsing:
          idx = [i for i, x in enumerate(self.rsids()) if x == rsid]
          return [self[i] for i in idx]
      
      raise ValueError("can't get variant without fully loading the bgen, or indexing")
    
    def at_position(self, pos):
      ''' get BgenVar from file given a position
      '''
      if not self.is_open == True:
          raise ValueError('bgen file is closed')
      
      if self.index:
          offsets = self.index.offset_by_pos(pos)
          return [BgenVar(self.handle, int(offset), self.thisptr.header.layout,
                          self.thisptr.header.compression, self.thisptr.header.nsamples,
                          self.is_stdin, self.is_open) 
                  for offset in offsets]
      
      if not self.delay_parsing:
          idx = [i for i, x in enumerate(self.positions) if x == pos]
          return [self[i] for i in idx]
      
      raise ValueError("can't get variant without fully loading the bgen, or indexing")
    
    def varids(self):
      ''' get the variant IDs of all variants in the bgen file
      '''
      if not self.is_open == True:
          raise ValueError("bgen file is closed")
      
      if self.index:
          raise ValueError("can't load varids when using an index file")
      
      varids = self.thisptr.varids()
      return [x.decode('utf8') for x in varids]
    
    def rsids(self):
      ''' get the rsIDs of all variants in the bgen file
      '''
      if not self.is_open == True:
          raise ValueError("bgen file is closed")
      
      if self.index:
          return self.index.rsids
      
      rsids = self.thisptr.rsids()
      return [x.decode('utf8') for x in rsids]
    
    def chroms(self):
        ''' get the chromosomes of all variants in the bgen file
        '''
        if not self.is_open == True:
            raise ValueError("bgen file is closed")
        
        if self.index:
            return self.index.chroms
        
        chroms = self.thisptr.chroms()
        return [x.decode('utf8') for x in chroms]
    
    def positions(self):
        ''' get the positions of all variants in the bgen file
        '''
        if not self.is_open == True:
            raise ValueError("bgen file is closed")
        
        if self.index:
            return self.index.positions
        
        return self.thisptr.positions()
    
    def __enter__(self):
        return self
    
    def __exit__(self, exc_type, exc_value, traceback):
        self.close()
        return False
    
    def close(self):
        if self.is_open == True:
            del self.thisptr
            self.handle = None
            if self.index:
                self.index.close()
            self.index = None
            self.is_open.off()

BgenFile = BgenReader