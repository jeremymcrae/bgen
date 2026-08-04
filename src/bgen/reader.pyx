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

# Random access needs to seek to a variant's offset, which a stream cannot do, so the
# only way through a pipe is to iterate. Kept in one place because several methods hit
# the same wall and had been describing it differently, or blaming the file instead.
NO_RANDOM_ACCESS = ('cannot pick out single variants while reading from stdin, since '
                    'that needs to seek within the bgen. Iterate over the bgen to '
                    'reach the variants in order, or read from a file instead')

cdef extern from "<iostream>" namespace "std::ios_base":
    cdef cppclass open_mode:
        pass
    cdef open_mode binary
    cdef cppclass iostate:
        pass
    cdef iostate badbit

cdef extern from "<iostream>" namespace "std":
    cdef cppclass istream:
        istream() except +
        void setstate(iostate state) except +

cdef extern from "<memory>" namespace "std":
    cdef cppclass shared_ptr_istream "std::shared_ptr<std::istream>":
        shared_ptr_istream() except +
        shared_ptr_istream(const shared_ptr_istream &) except +
        istream * get()

cdef extern from 'utils.h' namespace 'bgen':
    shared_ptr_istream borrowed_stream(istream * handle) except +

cdef extern from 'variant.h' namespace 'bgen':
    cdef cppclass Variant:
        # declare class constructor and methods
        Variant(shared_ptr_istream handle, uint64_t & offset, int layout, int compression, int expected_n, bool is_stdin) except +
        Variant() except +
        void minor_allele_dosage(float * dosage) except +
        void alt_dosage(float * dosage) except +
        string get_minor_allele() except +
        void probs_1d(float * dosage) except +
        int probs_per_sample() except +
        bool phased() except +
        uint8_t * ploidy() except +
        vector[uint8_t] copy_data() except +
        
        # declare public attributes
        string varid, rsid, chrom
        int pos
        long offset
        uint64_t next_variant_offset
        vector[string] alleles
        shared_ptr_istream handle

cdef extern from 'samples.h' namespace 'bgen':
    cdef cppclass Samples:
        Samples(istream * handle, int n_samples, uint64_t file_size) except +
        Samples(string path, int n_samples) except +
        Samples(int n_samples) except +
        Samples() except +
        
        # declare public methods
        const vector[string]& get_samples() except +

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
        void parse_all_variants() except +
        Variant & operator[](int idx) except +
        Variant & get(int idx) except +
        void drop_variants(vector[int] indices) except +
        vector[string] varids() except +
        vector[string] rsids() except +
        vector[string] chroms() except +
        vector[uint32_t] positions() except +
        
        # declare public attributes
        shared_ptr_istream handle
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
    
    This holds a share of the bgen stream, so that the file stays open for as
    long as any BgenVar still needs to read from it, even if the BgenReader it
    came from has been closed.
    '''
    cdef shared_ptr_istream ptr
    def __cinit__(self, uint64_t ptr):
        # a stream passed in as a raw address (i.e. from unpickling) belongs to
        # some other reader, so hold it without taking ownership
        self.ptr = borrowed_stream(<istream*>ptr)
    
    def __str__(self):
        return f'std::istream at {<uint64_t>self.ptr.get()}'
    
    def __dealloc__(self):
        pass
    
    def __reduce__(self):
        return (self.__class__, (<uint64_t>self.ptr.get(), ))

cdef IStream wrap_stream(shared_ptr_istream ptr):
    ''' build an IStream which shares ownership of an already open stream
    '''
    cdef IStream stream = IStream.__new__(IStream, 0)
    stream.ptr = ptr
    return stream

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

# compression flags as stored in the bgen header, mapped to the names the
# BgenWriter accepts, so a variant copied between files can be checked. this is
# cdef so it stays private to the module, and shared by BgenHeader and BgenVar
cdef dict COMPRESSION_FORMATS = {0: None, 1: 'zlib', 2: 'zstd'}

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
    def __cinit__(self, uint32_t offset, uint32_t nvariants, uint32_t nsamples,
            int compression, int layout, bool has_sample_ids, string metadata):
        self._offset = offset
        self._nvariants = nvariants
        self._nsamples = nsamples
        self._compression = compression
        self._layout = layout
        self._has_sample_ids = has_sample_ids
        self._metadata = metadata
    
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
    def compression(self): return COMPRESSION_FORMATS[self._compression]
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
    cdef int _layout, _compression, expected_n
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
        self._layout = layout
        self._compression = compression
        self.expected_n = expected_n
        self.is_stdin = is_stdin
        self.is_open = is_open
        
        # construct new Variant from the handle, offset and other file info
        self.thisptr = new Variant(self.handle.ptr, offset, layout, compression, expected_n, is_stdin)
    
    def __repr__(self):
       return f'BgenVar("{self.varid}", "{self.rsid}", "{self.chrom}", {self.pos}, {self.alleles})'
    
    def __str__(self):
       return f'{self.rsid} - {self.chrom}:{self.pos} {self.alleles}'
    
    @property
    def layout(self):
        ''' bgen layout the variant data is encoded with (1 or 2)
        '''
        return self._layout
    
    @property
    def compression(self):
        ''' compression scheme of the variant data (None, 'zlib' or 'zstd')
        '''
        return COMPRESSION_FORMATS[self._compression]
    
    @property
    def n_samples(self):
        ''' number of samples the variant holds genotypes for
        '''
        return self.expected_n
    
    cdef tuple _init_args(self):
        ''' the arguments needed to rebuild an equivalent BgenVar
        '''
        return (self.handle, self.thisptr.offset, self._layout, self._compression,
                self.expected_n, self.is_stdin, self.is_open)
    
    def __reduce__(self):
        ''' enable pickling of a BgenVar object
        '''
        warnings.warn("pickling BgenVar - make sure their BgenReader objects exist when unpickling", RuntimeWarning)
        return (self.__class__, self._init_args())
    
    def __copy__(self):
        ''' copy a BgenVar without hitting the __reduce__ warning
        '''
        return self.__class__(*self._init_args())
    
    def __deepcopy__(self, memo):
        ''' deepcopy a BgenVar
        
        The bgen file is shared with the original rather than duplicated, since
        the point of a BgenVar is to read from an open bgen file.
        '''
        cdef BgenVar var = self.__class__(*self._init_args())
        memo[id(self)] = var
        return var
    
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
        shared status, and set the istream badbit. Later cpp calls which would
        read from the stale file handle will raise ValueErrors instead due to
        checking the badbit before attempting to read from the istream.
        
        The badbit is used rather than the failbit, since the failbit is also
        set by ordinary reads which run past the end of the file (e.g. when
        iteration reaches the last variant), and those are recoverable, whereas
        a closed bgen must stay unreadable.
        
        However, this does not work for pickled BgenVars, since we cannot know
        whether the istream object exists or not.
        '''
        if not deref(self.is_open.ptr):
            self.thisptr.handle.get().setstate(badbit)
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
        return np.asarray(arr)
    @property
    def minor_allele(self):
        ''' get the minor allele of a biallelic variant
        '''
        self.__check_closed()
        return self.thisptr.get_minor_allele().decode('utf8')
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
                ragged = np.empty((len(ploidy), max_ploidy * cols), dtype=np.float32)
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
        if self.is_stdin:
            raise ValueError('cannot copy variant data directly while reading from stdin')
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
    # variants returned by __next__ so far, to spot a truncated bgen which stops
    # short of the variant count in the header
    cdef uint64_t n_iterated
    def __cinit__(self, path, sample_path='', bool delay_parsing=False):
        if isinstance(path, Path):
            path = str(path)
        if isinstance(sample_path, Path):
            sample_path = str(sample_path)
        self.is_stdin = self.__is_from_stdin(path)
        if self.is_stdin:
            delay_parsing = True
            path = '/dev/stdin'
        
        if Path(path).exists() and Path(path).is_dir():
            raise ValueError(f'bgen path is for a folder: {path}')
        
        delay_parsing |= self._check_for_index(path)
        
        self.path = path.encode('utf8')
        self.sample_path = sample_path.encode('utf8')
        self.delay_parsing = delay_parsing
        
        samp = '' if sample_path == '' else f', (samples={self.sample_path.decode("utf")})'
        logging.debug(f'opening BgenFile from {self.path.decode("utf")}{samp}')
        self.thisptr = new CppBgenReader(self.path, self.sample_path, self.delay_parsing)
        self.handle = wrap_stream(self.thisptr.handle)
        self.is_open = OpenStatus()
        self.offset = self.thisptr.offset
        self.n_iterated = 0
    
    def __is_from_stdin(self, bgen_path):
        if bgen_path is sys.stdin:
            return True
        elif str(bgen_path) == '/dev/stdin':
            return True
        elif str(bgen_path) == '-':
            return True
        return False
    
    def __dealloc__(self):
        # a deallocated reader has nobody left to report errors to, so drop them
        # rather than emitting 'Exception ignored in __dealloc__' noise. Call
        # close() directly (or use a with block) to have errors raised.
        try:
            self.close()
        except Exception:
            pass
    
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
            self.n_iterated += 1
            return var
        except IndexError:
            # Running out of file before reaching the header's variant count means
            # the bgen is truncated. Without this, iteration would just stop early
            # and a truncated file would look like a complete, shorter one.
            if self.n_iterated < self.thisptr.header.nvariants:
                raise ValueError(f'bgen is truncated - the header lists '
                                 f'{self.thisptr.header.nvariants} variants, but only '
                                 f'{self.n_iterated} could be read')
            raise StopIteration
    
    def __len__(self):
      ''' number of variants in the bgen file
      
      With delay_parsing this is the count from the bgen header, which is all that
      is knowable without reading the file. A truncated bgen still names its
      original number of variants there, so this can be more than can actually be
      read - iterating, or getting a variant by index, raises in that case.
      '''
      if not self.is_open == True:
          raise ValueError("bgen file is closed")
      
      length = self.thisptr.variants.size()
      if length > 0:
          return length
      
      if self.index is not None:
          # variants are looked up through the index, so its count is the one that
          # matches what indexing this BgenReader can actually reach
          return len(self.index)
      
      return self.thisptr.header.nvariants
    
    def __getitem__(self, Py_ssize_t idx):
        ''' pull out a Variant by index position
        '''
        cdef Py_ssize_t orig_idx, size
        
        if not self.is_open == True:
            raise ValueError('bgen file is closed')
        
        orig_idx = idx
        size = len(self)
        if idx < 0:
            idx += size
        
        if idx >= size or idx < 0:
            raise IndexError(f'cannot get Variant at index: {orig_idx}')
        
        # Getting a variant by index has to seek to it, which a pipe cannot do. Say so
        # before parsing, since otherwise the seek fails after the walk and the error
        # blames the file for being truncated, which sends people looking for
        # corruption in a bgen that is perfectly fine.
        if self.is_stdin:
            raise ValueError(NO_RANDOM_ACCESS)
        
        # account for lazy loading variants from bgen
        if self.index is None and self.thisptr.variants.size() == 0:
            self.thisptr.parse_all_variants()
        
        cdef long offset
        offset = self.index.offset_by_index(idx) if self.index else self.thisptr.variants[idx].offset
        try:
            return BgenVar(self.handle, offset, self.thisptr.header.layout,
              self.thisptr.header.compression, self.thisptr.header.nsamples,
              self.is_stdin, self.is_open)
        except IndexError:
            # The index was in range, so running out of file means the bgen does
            # not hold the variant its index (or header) says it should. Reporting
            # the end of the file for a valid index just looks like a bug here.
            raise ValueError(f'bgen is truncated - could not read the variant at '
                             f'index {orig_idx}')
    
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
      
      samples = self.thisptr.samples.get_samples()
      return [x.decode('utf8') for x in samples]
    
    def drop_variants(self, list indices):
        ''' drops variants from bgen by indices, for avoiding processing variants
        
        Raises IndexError if any index is negative or beyond the last variant,
        and ValueError if any index is duplicated. Nothing is dropped unless
        every index is valid.
        
        .. deprecated::
            This does not drop variants consistently, and will be removed. Only
            len() reflects the dropped variants - iterating the BgenReader still
            yields every variant, calling rsids()/varids()/chroms()/positions()
            re-parses the bgen and so undoes the drop, and indexed access only
            respects the drop if no bgenix (.bgi) index file is present. Filter
            the variants in python instead, e.g. by skipping the unwanted ones
            while iterating, or by selecting indices with bfile[i].
        '''
        if not self.is_open == True:
            raise ValueError("bgen file is closed")
        
        # NOTE: no stacklevel here. A compiled cython method does not push a
        # python frame, so the default stacklevel=1 already attributes the
        # warning to the calling python code. Passing stacklevel=2 would skip
        # past the caller, report the location as 'sys:1', and stop the default
        # warning filters from displaying it at all.
        warnings.warn('BgenReader.drop_variants is deprecated and will be '
                      'removed - it does not drop variants consistently. '
                      'Filter variants in python instead.',
                      DeprecationWarning)
        
        if self.delay_parsing:
            self.thisptr.parse_all_variants()
        
        self.thisptr.drop_variants(indices)
    
    def fetch(self, chrom, start=None, stop=None):
        ''' fetches all variants within a genomic region
        
        Args:
            chrom: chromosome that variants must be on
            start: start nucleotide of region. If None, gets variants with
                positions up to stop, or all variants on the chromosome if stop
                is also None
            stop: end nucleotide of region. If None, gets variants with positions after start
        
        Yields:
            BgenVars for variants within the genome region
        '''
        if not self.is_open == True:
            raise ValueError('bgen file is closed')
        
        if not self.index:
            if self.is_stdin:
                raise ValueError(NO_RANDOM_ACCESS)
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
      
      if self.is_stdin:
          raise ValueError(NO_RANDOM_ACCESS)
      
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
          idx = [i for i, x in enumerate(self.positions()) if x == pos]
          return [self[i] for i in idx]
      
      if self.is_stdin:
          raise ValueError(NO_RANDOM_ACCESS)
      
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
        if not self.is_open == True:
            # this can also be a partially constructed reader, where __cinit__
            # raised before is_open was set, so still release the index
            self._close_index()
            return
        
        # mark the bgen as closed before releasing anything, so that an error
        # part way through cannot leave this reporting itself as open while
        # holding a freed pointer, and so a second close() cannot delete the
        # pointer twice
        self.is_open.off()
        try:
            del self.thisptr
            self.thisptr = NULL
            self.handle = None
        finally:
            self._close_index()
    
    def _close_index(self):
        ''' release the bgenix index, if one was opened
        '''
        index = self.index
        self.index = None
        if index is not None:
            index.close()

BgenFile = BgenReader