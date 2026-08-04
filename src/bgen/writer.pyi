import os
from typing import Any, Iterable, Sequence

import numpy as np
from numpy.typing import ArrayLike, NDArray

from bgen.reader import BgenVar

class Indexer:
    ''' class to automatically index bgen files as they are being constructed
    '''
    conn: Any
    cur: Any
    def __init__(self, bgen_path: str | os.PathLike[str]) -> None: ...
    def create_tables(self) -> None: ...
    def add_variant(self,
                    chrom: str,
                    pos: int,
                    rsid: str,
                    alleles: Sequence[str],
                    offset: int,
                    size: int,
                    ) -> None: ...
    def close(self) -> None: ...

class BgenWriter:
    ''' class to write bgen files to disk
    '''
    def __new__(cls,
                path: str | os.PathLike[str],
                n_samples: int,
                samples: Iterable[str] | None = ...,
                compression: str | None = 'zstd',
                layout: int = 2,
                metadata: str | None = None,
                ) -> BgenWriter: ...
    def __repr__(self) -> str: ...
    def __enter__(self) -> BgenWriter: ...
    def __exit__(self, exc_type: Any, exc_value: Any, traceback: Any) -> bool: ...
    def add_variant(self,
                    varid: str,
                    rsid: str,
                    chrom: str,
                    pos: int,
                    alleles: Sequence[str],
                    genotypes: ArrayLike,
                    ploidy: int | NDArray[np.integer[Any]] = 2,
                    phased: bool = False,
                    bit_depth: int = 8,
                    ) -> None:
        ''' add a variant to the bgen file on disk

        Args:
            varid: variant ID
            rsid: reference SNP ID
            chrom: chromosome the variant is on
            pos: nucleotide position of the variant
            alleles: list of allele strings. Duplicates are allowed, but warned about
            genotypes: numpy array of genotype proabilities, ordered as per the
                bgen samples.
            ploidy: integer for constant ploidy, or numpy array of ploidy values per
                sample, in same order as genotypes
            phased: whether the genotypes are for phased data or not
            bit_depth: integer from 1-32 (inclusive) for how many bits to store
                each genotype in. Stored probabilities step by 1/(2**bit_depth - 1),
                so the default of 8 keeps two decimal places. Depths that lose more
                than that warn. Depths above 24 only cost space, since probabilities
                are read back as float32.
        '''
        ...
    def add_variant_direct(self, variant: BgenVar) -> None:
        ''' insert a BgenVar directly into the bgen file
        '''
        ...
    def _check_compatible(self, variant: BgenVar) -> None:
        ''' check a variant's data can be copied into this file unchanged
        '''
        ...
    def close(self) -> None: ...
