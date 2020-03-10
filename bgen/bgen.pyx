# cython: language_level=3, boundscheck=False

from libcpp.string cimport string
from libcpp.vector cimport vector
from libc.stdint cimport uint32_t

from cython.operator cimport dereference as deref

cdef extern from "<iostream>" namespace "std":
    cdef cppclass istream:
        istream& read(char *, int) except +

cdef extern from "<iostream>" namespace "std::ios_base":
    cdef cppclass open_mode:
        pass
    cdef open_mode binary

cdef extern from "<fstream>" namespace "std":
    cdef cppclass ifstream(istream):
        ifstream(const char*) except +
        ifstream(const char*, open_mode) except +

cdef extern from 'variant.h' namespace 'bgen':
    cdef cppclass Variant:
        # declare class constructor and methods
        Variant(ifstream & handle, int layout, int compression, int expected_n) except +
        Variant() except +
        vector[float] alt_dosage()
        
        # declare public attributes
        string varid, rsid, chrom
        int pos
        vector[string] alleles

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

cdef class BgenVar:
    cdef string _varid
    cdef string _rsid
    cdef string _chrom
    cdef int _pos
    cdef vector[string] _alleles
    cdef vector[float] _alt_dosage
    def __cinit__(self, string varid, string rsid, string chrom, int pos,
            vector[string] alleles, vector[float] alt_dosage):
        self._varid = varid
        self._rsid = rsid
        self._chrom = chrom
        self._pos = pos
        self._alleles = alleles
        self._alt_dosage = alt_dosage
    
    def __repr__(self):
        return f'BgenVar("{self.varid}", "{self.rsid}", "{self.chrom}", x{self.pos}, {self.alleles})'
    
    def __str__(self):
        return f'{self.rsid} - {self.chrom}:{self.pos} {self.alleles}'
    
    @property
    def varid(self):
      return self._varid.decode('utf8')
    @property
    def rsid(self):
        return self._rsid.decode('utf8')
    @property
    def chrom(self):
        return self._chrom.decode('utf8')
    @property
    def pos(self):
        return self._pos
    @property
    def alleles(self):
        return [x.decode('utf8') for x in self._alleles]
    @property
    def alt_dosage(self):
        return self._alt_dosage

cdef class BgenFile:
    cdef Bgen * thisptr
    cdef string path, sample_path
    def __cinit__(self, str path, str sample_path=''):
        self.path = path.encode('utf8')
        self.sample_path = sample_path.encode('utf8')
        self.thisptr = new Bgen(self.path, self.sample_path)
    
    def __dealloc__(self):
        del self.thisptr
    
    def __repr__(self):
        return f'BgenFile("{self.path.decode("utf8")}", "{self.sample_path.decode("utf8")}")'
    
    def __iter__(self):
        length = self.thisptr.variants.size()
        for idx in range(length):
            yield self[idx]
    
    def __getitem__(self, int idx):
        ''' pull out a Variant by index position
        '''
        variant = self.thisptr.get(idx)
        # print(variant.varid)
        return BgenVar(variant.varid, variant.rsid, variant.chrom, variant.pos,
            variant.alleles, variant.alt_dosage())
    
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
  
