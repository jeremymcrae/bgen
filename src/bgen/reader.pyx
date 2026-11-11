# cython: language_level=3, boundscheck=False, emit_linenums=True

import logging
from pathlib import Path
import sys
import warnings

from libcpp cimport bool
from libcpp.string cimport string
from libcpp.vector cimport vector
from libc.stdint cimport uint8_t, uint32_t, uint64_t
from libc.string cimport memcpy

import numpy as np

from bgen.index import Index

# Random access needs to seek to a variant's offset, which a stream cannot do, so the
# only way through a pipe is to iterate. Kept in one place because several methods hit
# the same wall and had been describing it differently, or blaming the file instead.
NO_RANDOM_ACCESS = ('cannot pick out single variants while reading from stdin, since '
                    'that needs to seek within the bgen. Iterate over the bgen to '
                    'reach the variants in order, or read from a file instead')

cdef extern from "<iostream>" namespace "std::ios_base":
    cdef cppclass iostate:
        pass
    cdef iostate badbit

cdef extern from "<iostream>" namespace "std":
    cdef cppclass istream:
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
        Variant(shared_ptr_istream handle, uint64_t offset, int layout, int compression, int expected_n, bool is_stdin) except +
        void minor_allele_dosage(float * dosage) except +
        void alt_dosage(float * dosage) except +
        string get_minor_allele() except +
        void probs_1d(float * dosage) except +
        int probs_per_sample() except +
        bool phased() except +
        bool probs_above_max() except +
        uint8_t * ploidy() except +
        vector[uint8_t] copy_data() except +
        
        # declare public attributes
        string varid, rsid, chrom
        uint32_t pos
        uint64_t offset
        uint64_t next_variant_offset
        vector[string] alleles
        shared_ptr_istream handle

cdef extern from 'samples.h' namespace 'bgen':
    cdef cppclass Samples:
        const vector[string]& get_samples() except +

cdef extern from 'header.h' namespace 'bgen':
    cdef cppclass Header:
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
        void close_stream() except +
        uint64_t first_variant_offset()
        void parse_all_variants() except +
        void drop_variants(vector[int] indices) except +
        
        # declare public attributes
        shared_ptr_istream handle
        vector[Variant] variants
        Samples samples
        Header header

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
    
    A BgenReader and the BgenVariant objects opened from it hold the same
    OpenStatus. When a bgen file is closed, the status changes for all loaded
    variants, in order to prevent reading more data from the closed file.
    
    This may not work for pickled BgenVariant objects (after unpickling), as
    the BgenReader object may have been closed after the BgenVariant was 
    pickled. In order to allow pickled BgenVariants, as part of unpickling 
    the BgenVariant, it creates a new OpenStatus object where status=True.
    '''
    cdef bint status
    def __cinit__(self):
        self.status = True
    def __str__(self):
        return f'status={self.status}'
    def __eq__(self, other):
        return self.status == other
    def off(self):
        self.status = False
    def __reduce__(self):
        return (self.__class__, ())

# compression flags as stored in the bgen header, mapped to the names the
# BgenWriter accepts, so a variant copied between files can be checked. this is
# cdef so it stays private to the module, and shared by BgenHeader and BgenVar
cdef dict COMPRESSION_FORMATS = {0: None, 1: 'zlib', 2: 'zstd'}

cdef class BgenHeader:
    ''' holds information about the Bgen file, obtained from the intial header.
    '''
    cdef readonly uint32_t offset
    cdef readonly uint32_t nvariants
    cdef readonly uint32_t nsamples
    cdef int _compression
    cdef readonly int layout
    cdef readonly bool has_sample_ids
    cdef string _metadata
    def __cinit__(self, uint32_t offset, uint32_t nvariants, uint32_t nsamples,
            int compression, int layout, bool has_sample_ids, string metadata):
        self.offset = offset
        self.nvariants = nvariants
        self.nsamples = nsamples
        self._compression = compression
        self.layout = layout
        self.has_sample_ids = has_sample_ids
        self._metadata = metadata
    
    def __repr__(self):
        return f'BgenHeader(offset={self.offset}, nvariants={self.nvariants}, ' \
            f'nsamples={self.nsamples}, compression={self.compression}, ' \
            f'layout={self.layout}, has_sample_ids={self.has_sample_ids})'
    
    @property
    def compression(self): return COMPRESSION_FORMATS[self._compression]
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
        
        We prevent this by checking a status shared by all the BgenVars
        which were opened by a given BgenReader. This status will have
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
        if not self.is_open.status:
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
    def __warn_if_malformed(self):
        ''' warn if the decode just done found probabilities summing above the maximum

        The bgen spec has the probabilities of a sample sum to the largest value the bit
        depth can store, and the final value of each group is not stored but inferred as
        the remainder. Poorly encoded files can have negative remainders, so this warns if
        it happens.
        '''
        if self.thisptr.probs_above_max():
            varid = self.thisptr.varid.decode('utf8')
            rsid = self.thisptr.rsid.decode('utf8')
            logging.warning(f'variant {rsid}/{varid} stores genotype probabilities which '
                            f'sum to more than the bit depth allows, so this bgen is '
                            f'malformed. The probability inferred for the affected samples '
                            f'is not reliable, and may be negative or clamped to zero')
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
        # parse first, to validate the sample count before it sizes the array
        self.thisptr.probs_per_sample()
        cdef float[:] dose = np.empty(self.expected_n, dtype=np.float32, order='C')
        self.thisptr.minor_allele_dosage(&dose[0])
        self.__warn_if_malformed()
        return np.asarray(dose)
    @property
    def alt_dosage(self):
        ''' dosage for the alt allele for a biallelic variant
        '''
        self.__check_closed()
        # parse first, to validate the sample count before it sizes the array
        self.thisptr.probs_per_sample()
        cdef float[:] dose = np.empty(self.expected_n, dtype=np.float32, order='C')
        self.thisptr.alt_dosage(&dose[0])
        self.__warn_if_malformed()
        return np.asarray(dose)
    @property
    def probabilities(self):
        ''' get the allelic probabilities for a variant
        '''
        self.__check_closed()
        cdef int cols = self.thisptr.probs_per_sample()
        cdef uint32_t n_samples = self.expected_n
        cdef uint64_t size = <uint64_t> n_samples * cols
        cdef uint8_t * ploidy = NULL
        cdef bool phased = self.thisptr.phased()
        if phased:
            ploidy = self.thisptr.ploidy()
            size = fast_ploidy_sum(ploidy, n_samples) * cols
        
        cdef float[::1] arr = np.empty(size, dtype=np.float32, order='C')
        self.thisptr.probs_1d(&arr[0])
        self.__warn_if_malformed()

        cdef Range ploidy_range
        cdef float[:, ::1] ragged
        cdef uint64_t used = 0
        cdef size_t count
        cdef uint32_t i
        if phased:
            ploidy_range = fast_range(ploidy, n_samples)
            max_ploidy = ploidy_range._max
            min_ploidy = ploidy_range._min
            if min_ploidy == max_ploidy:
                # quickly reshape probs if ploidy is constant
                width = max_ploidy * cols
                data = np.reshape(arr, (-1, width))
            else:
                # phased data comes as one row per haplotype. Copy each sample's
                # haplotypes into a single row, leaving nan beyond its ploidy
                data = np.full((n_samples, max_ploidy * cols), np.nan, dtype=np.float32)
                ragged = data
                for i in range(n_samples):
                    count = <size_t> ploidy[i] * cols
                    if count > 0:
                        memcpy(&ragged[i, 0], &arr[used], count * sizeof(float))
                    used += count
        else:
            data = np.reshape(arr, (-1, cols))
        
        return data
    
    def copy_data(self):
        ''' get a copy of the data on disk for the variant, as bytes
        
        This can be used to quickly copy data from one bgen into another, if
        you have the same set of samples between the source and destination
        bgens. This primarily avoids decompressing, decoding, re-encoding, and
        compressing the genotype data.
        '''
        if self.is_stdin:
            raise ValueError('cannot copy variant data directly while reading from stdin')
        self.__check_closed()
        cdef vector[uint8_t] data = self.thisptr.copy_data()
        if data.size() == 0:
            return b''
        return (<const char *> data.data())[:data.size()]

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
    def __cinit__(self, path, sample_path='', bool delay_parsing=True):
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
        self.offset = self.thisptr.first_variant_offset()
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
    
    cdef int _check_open(self) except -1:
        ''' raise an error if the bgen file has been closed
        '''
        if self.is_open is None or not self.is_open.status:
            raise ValueError('bgen file is closed')
        return 0
    
    cdef BgenVar _variant_at(self, uint64_t offset):
        ''' open the variant which starts at a byte offset in the bgen
        '''
        return BgenVar(self.handle, offset, self.thisptr.header.layout,
                       self.thisptr.header.compression, self.thisptr.header.nsamples,
                       self.is_stdin, self.is_open)
    
    def __iter__(self):
        return self
    
    def __next__(self):
        ''' iterate through all variants in the bgen file
        '''
        self._check_open()
        
        cdef BgenVar var
        try:
            var = self._variant_at(self.offset)
            # read the next offset off the C++ variant, rather than through the python
            # property, which would box it into a python int just to unbox it
            self.offset = var.thisptr.next_variant_offset
            self.n_iterated += 1
            return var
        except IndexError:
            # running out of file before the header's variant count means the bgen
            # is truncated
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
      self._check_open()
      
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
        
        self._check_open()
        
        orig_idx = idx
        size = len(self)
        if idx < 0:
            idx += size
        
        if idx >= size or idx < 0:
            raise IndexError(f'cannot get Variant at index: {orig_idx}')
        
        # getting a variant by index has to seek to it, which a pipe cannot do
        if self.is_stdin:
            raise ValueError(NO_RANDOM_ACCESS)
        
        # account for lazy loading variants from bgen
        if self.index is None and self.thisptr.variants.size() == 0:
            self.thisptr.parse_all_variants()
        
        cdef uint64_t offset
        offset = self.index.offset_by_index(idx) if self.index else self.thisptr.variants[idx].offset
        try:
            return self._variant_at(offset)
        except IndexError:
            # the index was in range, so the bgen lacks a variant that its index (or
            # header) lists
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
      self._check_open()
      
      hdr = self.thisptr.header
      return BgenHeader(hdr.offset, hdr.nvariants, hdr.nsamples,
          hdr.compression, hdr.layout, hdr.has_sample_ids, hdr.extra)
    
    @property
    def samples(self):
      ''' get list of samples in the bgen file
      '''
      self._check_open()
      
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
        self._check_open()
        
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
        self._check_open()
        
        if not self.index:
            if self.is_stdin:
                raise ValueError(NO_RANDOM_ACCESS)
            raise ValueError("can't fetch variants without index")
        
        for offset in self.index.fetch(chrom, start, stop):
            yield self._variant_at(offset)
    
    cdef list _lookup(self, key, bint by_rsid):
        ''' get the BgenVars whose rsID, or else position, matches a key
        '''
        self._check_open()
        
        if self.index:
            if by_rsid:
                offsets = self.index.offset_by_rsid(key)
            else:
                offsets = self.index.offset_by_pos(key)
            return [self._variant_at(int(offset)) for offset in offsets]
        
        if not self.delay_parsing:
            values = self.rsids() if by_rsid else self.positions()
            return [self[i] for i, x in enumerate(values) if x == key]
        
        if self.is_stdin:
            raise ValueError(NO_RANDOM_ACCESS)
        
        raise ValueError("can't get variant without fully loading the bgen, or indexing")
    
    def with_rsid(self, rsid):
      ''' get BgenVar from file given an rsID
      '''
      return self._lookup(rsid, True)
    
    def at_position(self, pos):
      ''' get BgenVar from file given a position
      '''
      return self._lookup(pos, False)
    
    def varids(self):
      ''' get the variant IDs of all variants in the bgen file
      '''
      self._check_open()
      
      if self.index:
          raise ValueError("can't load varids when using an index file")
      
      self.thisptr.parse_all_variants()
      return [self.thisptr.variants[i].varid.decode('utf8')
              for i in range(self.thisptr.variants.size())]
    
    def rsids(self):
      ''' get the rsIDs of all variants in the bgen file
      '''
      self._check_open()
      
      if self.index:
          return self.index.rsids
      
      self.thisptr.parse_all_variants()
      return [self.thisptr.variants[i].rsid.decode('utf8')
              for i in range(self.thisptr.variants.size())]
    
    def chroms(self):
        ''' get the chromosomes of all variants in the bgen file
        '''
        self._check_open()
        
        if self.index:
            return self.index.chroms
        
        self.thisptr.parse_all_variants()
        return [self.thisptr.variants[i].chrom.decode('utf8')
                for i in range(self.thisptr.variants.size())]
    
    def positions(self):
        ''' get the positions of all variants in the bgen file
        '''
        self._check_open()
        
        if self.index:
            return self.index.positions
        
        self.thisptr.parse_all_variants()
        return [self.thisptr.variants[i].pos for i in range(self.thisptr.variants.size())]
    
    def __enter__(self):
        return self
    
    def __exit__(self, exc_type, exc_value, traceback):
        self.close()
        return False
    
    def close(self):
        if self.is_open is None or not self.is_open.status:
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
            # release the file itself, not just this reader's reference to it. Any
            # Variant taken from this reader shares the stream and may outlive the
            # reader, which would otherwise leave the file open, and windows will
            # not let an open file be deleted
            if self.thisptr != NULL:
                self.thisptr.close_stream()
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