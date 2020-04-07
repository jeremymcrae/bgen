
from pathlib import Path

import numpy as np

def load_gen_data():
    ''' load data from "example.gen" for comparison with the bgen files
    '''
    class GenVar:
        ''' store data for easy comparison with BgenVars
        '''
        def __init__(self, chrom, varid, rsid, pos, ref, alt, probs):
            self.chrom = chrom
            self.varid = varid
            self.rsid = rsid
            self.pos = int(pos)
            self.alleles = [ref, alt]
            self.probabilities = probs
        def __eq__(self, other):
            return self.chrom == other.chrom and self.pos == other.pos and \
                self.varid == other.varid and self.rsid == other.rsid and \
                set(self.alleles) == set(other.alleles)
    
    variants = []
    path = Path(__file__).parent /  "data" / "example.gen"
    with open(path, 'rt') as gen:
        for line in gen:
            chrom, varid, rsid, pos, ref, alt, *probs = line.strip('\n').split()
            probs = np.array(list(map(float, probs)))
            probs = np.reshape(probs, (-1, 3))
            nonzero = (probs == 0.0).all(axis=1)
            probs[nonzero] = float('nan')
            var = GenVar(chrom, varid, rsid, pos, ref, alt, probs)
            variants.append(var)
    return variants

def epsilon(bit_depth):
    ''' get the max difference expected given the bit depth used for the probs
    '''
    return 2 / (2 ** (bit_depth - 1))

def array_delta(truth, parsed):
    ''' get the max absolute difference between two numpy arrays of equal shape
    '''
    deltas = truth - parsed
    deltas = deltas.flatten()
    deltas = deltas[~np.isnan(deltas)]
    return abs(deltas).max()

def arrays_equal(truth, parsed, bit_depth):
    ''' check that two arrays are sufficiently equal
    '''
    eps_abs = 3.2e-5
    delta = array_delta(truth, parsed)
    return delta < epsilon(bit_depth) or delta < eps_abs
