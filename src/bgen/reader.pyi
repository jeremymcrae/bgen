import os
from typing import Any, IO, Iterator

import numpy as np
from numpy.typing import NDArray

class IStream:
    ''' basic cython implementation of std::istream, for easy pickling
    '''
    def __new__(cls, ptr: int) -> IStream: ...
    def __str__(self) -> str: ...
    def __reduce__(self) -> tuple[Any, ...]: ...

class OpenStatus:
    ''' class to share status of whether a bgen file is currently open
    '''
    def __new__(cls) -> OpenStatus: ...
    def __str__(self) -> str: ...
    def __eq__(self, other: object, /) -> bool: ...
    def off(self) -> None: ...
    def __reduce__(self) -> tuple[Any, ...]: ...

class BgenHeader:
    ''' holds information about the Bgen file, obtained from the intial header.
    '''
    def __new__(cls,
                offset: int,
                nvariants: int,
                nsamples: int,
                compression: int,
                layout: int,
                has_sample_ids: bool,
                metadata: bytes,
                ) -> BgenHeader: ...
    def __repr__(self) -> str: ...
    @property
    def offset(self) -> int: ...
    @property
    def nsamples(self) -> int: ...
    @property
    def nvariants(self) -> int: ...
    @property
    def compression(self) -> str | None:
        ''' compression scheme - one of None, 'zlib' or 'zstd'
        '''
        ...
    @property
    def layout(self) -> int: ...
    @property
    def has_sample_ids(self) -> bool: ...
    @property
    def metadata(self) -> str: ...

class BgenVar:
    ''' holds data for a Variant from a bgen file
    '''
    def __new__(cls,
                handle: IStream,
                offset: int,
                layout: int,
                compression: int,
                expected_n: int,
                is_stdin: bool,
                is_open: OpenStatus,
                ) -> BgenVar: ...
    def __repr__(self) -> str: ...
    def __str__(self) -> str: ...
    def __reduce__(self) -> tuple[Any, ...]: ...
    @property
    def varid(self) -> str: ...
    @property
    def rsid(self) -> str: ...
    @property
    def chrom(self) -> str: ...
    @property
    def pos(self) -> int: ...
    @property
    def fileoffset(self) -> int: ...
    @property
    def next_variant_offset(self) -> int: ...
    @property
    def alleles(self) -> list[str]: ...
    @property
    def is_phased(self) -> bool: ...
    @property
    def ploidy(self) -> NDArray[np.uint8]:
        ''' get the ploidy for each sample
        '''
        ...
    @property
    def minor_allele(self) -> str:
        ''' get the minor allele of a biallelic variant
        '''
        ...
    @property
    def minor_allele_dosage(self) -> NDArray[np.float32]:
        ''' dosage for the minor allele for a biallelic variant
        '''
        ...
    @property
    def alt_dosage(self) -> NDArray[np.float32]:
        ''' dosage for the alt allele for a biallelic variant
        '''
        ...
    @property
    def probabilities(self) -> NDArray[np.floating[Any]]:
        ''' get the allelic probabilities for a variant
        '''
        ...
    def copy_data(self) -> list[int]:
        ''' get a copy of the data on disk for the variant
        '''
        ...

class BgenReader:
    ''' class to open bgen files from disk, and access variant data within
    '''
    def __new__(cls,
                path: str | os.PathLike[str] | IO[Any],
                sample_path: str | os.PathLike[str] = '',
                delay_parsing: bool = False,
                ) -> BgenReader: ...
    def __repr__(self) -> str: ...
    def __iter__(self) -> Iterator[BgenVar]: ...
    def __next__(self) -> BgenVar: ...
    def __len__(self) -> int: ...
    def __getitem__(self, idx: int) -> BgenVar:
        ''' pull out a Variant by index position
        '''
        ...
    def __enter__(self) -> BgenReader: ...
    def __exit__(self, exc_type: Any, exc_value: Any, traceback: Any) -> bool: ...
    @property
    def header(self) -> BgenHeader:
        ''' get header info from bgen file
        '''
        ...
    @property
    def samples(self) -> list[str]:
        ''' get list of samples in the bgen file
        '''
        ...
    def drop_variants(self, indices: list[int]) -> None:
        ''' drops variants from bgen by indices, for avoiding processing variants
        '''
        ...
    def fetch(self,
              chrom: str,
              start: int | None = None,
              stop: int | None = None,
              ) -> Iterator[BgenVar]:
        ''' fetches all variants within a genomic region
        '''
        ...
    def with_rsid(self, rsid: str) -> list[BgenVar]:
        ''' get BgenVars from file given an rsID
        '''
        ...
    def at_position(self, pos: int) -> list[BgenVar]:
        ''' get BgenVars from file given a position
        '''
        ...
    def varids(self) -> list[str]:
        ''' get the variant IDs of all variants in the bgen file
        '''
        ...
    def rsids(self) -> list[str]:
        ''' get the rsIDs of all variants in the bgen file
        '''
        ...
    def chroms(self) -> list[str]:
        ''' get the chromosomes of all variants in the bgen file
        '''
        ...
    def positions(self) -> list[int]:
        ''' get the positions of all variants in the bgen file
        '''
        ...
    def close(self) -> None: ...

BgenFile = BgenReader
