
from pathlib import Path
import re

import numpy as np

def can_cap_memory():
    ''' whether this platform can actually cap a child's address space

    macOS defines RLIMIT_AS (aliased to RLIMIT_RSS in its headers), so the constant
    existing does not mean setrlimit() on it works. The kernel there rejects the call
    with EINVAL, which surfaces as "current limit exceeds maximum limit". Requiring
    /proc/self/status too, rather than probing setrlimit directly, also rules out
    sizing the cap as an absolute headroom, which is tight enough to hang OpenBLAS's
    thread init on a machine whose imports alone take more than that.
    '''
    try:
        import resource
    except ImportError:
        return False
    return hasattr(resource, 'RLIMIT_AS') and Path('/proc/self/status').exists()

def cap_after_import(headroom):
    ''' source for a child process to cap its own address space

    Returns code to run *after* the child's imports, so that the cap applies only to the
    allocation under test. A child that caps first hangs rather than failing: OpenBLAS
    sizes its thread arenas by core count, and when RLIMIT_AS is too tight to mmap them,
    its init retries forever. Importing bgen already needs several hundred MB, so how
    much room a fixed cap leaves depends on the machine, and a machine with more cores
    than the one a cap was chosen on crosses the line and hangs.

    Sizing the cap from what the imports actually took removes that dependency, leaving
    headroom bytes for the allocation under test. Only linux publishes VmSize, and
    elsewhere the cap falls back to being an absolute one.
    '''
    return f'''
import resource
used = 0
try:
    with open('/proc/self/status') as handle:
        for line in handle:
            if line.startswith('VmSize'):
                used = int(line.split()[1]) * 1024
except OSError:
    pass
resource.setrlimit(resource.RLIMIT_AS, (used + {headroom}, used + {headroom}))
'''

class GenVar:
    ''' store data for easy comparison with BgenVars
    '''
    def __init__(self, chrom, varid, rsid, pos, alleles, probs):
        self.chrom = chrom
        self.varid = varid
        self.rsid = rsid
        self.pos = int(pos)
        self.alleles = alleles
        self.probabilities = probs
    def __repr__(self):
        return f'{self.rsid} - {self.chrom}:{self.pos} {self.alleles}'
    def __eq__(self, other):
        return self.chrom == other.chrom and self.pos == other.pos and \
            self.varid == other.varid and self.rsid == other.rsid and \
            set(self.alleles) == set(other.alleles)

def load_gen_data():
    ''' load data from "example.gen" for comparison with the bgen files
    '''
    variants = []
    path = Path(__file__).parent /  "data" / "example.gen"
    with open(path, 'rt') as gen:
        for line in gen:
            chrom, varid, rsid, pos, ref, alt, *probs = line.strip('\n').split()
            probs = np.array(list(map(float, probs)))
            probs = np.reshape(probs, (-1, 3))
            nonzero = (probs == 0.0).all(axis=1)
            probs[nonzero] = float('nan')
            var = GenVar(chrom, varid, rsid, pos, [ref, alt], probs)
            variants.append(var)
    return variants

def parse_vcf_samples(format, samples):
    ''' parses sample data from VCF into probabilities and ploidy
    '''
    samples = [dict(zip(format.split(':'), x.split(':'))) for x in samples]
    keys = {'GT': re.compile('[/|]'), 'GP': re.compile(','), 'HP': re.compile(',')}
    for x in samples:
        for k in x:
            x[k] = re.split(keys[k], x[k])
    
    probs = [x['GP'] if 'GP' in format else x['HP'] for x in samples]
    max_len = max(len(x) for x in probs)
    
    for i in range(len(probs)):
        if len(probs[i]) < max_len:
            probs[i] += ['nan'] * (max_len - len(probs[i]))
    
    ploidy = np.array([len(x['GT']) for x in samples])
    return np.array(probs, dtype=float), ploidy

def load_vcf_data():
    '''load data from 'complex.vcf' for comparison with the complex bgen files
    '''
    variants = []
    path = Path(__file__).parent / 'data' / 'complex.vcf'
    with open(path, 'rt') as vcf:
        for line in vcf:
            if line.startswith('#'):
                continue
            chrom, pos, varid, ref, alts, _, _, _, format, *samples = line.strip('\n').split('\t')
            varid = varid.split(',')
            if len(varid) > 1:
                rsid, varid = varid
            else:
                varid, rsid = '', varid[0]
            probs, ploidy = parse_vcf_samples(format, samples)
            var = GenVar(chrom, varid, rsid, pos, [ref] + alts.split(','), probs)
            var.ploidy = ploidy
            variants.append(var)
    return variants

def load_haps_data():
    ''' load data from 'haplotypes.haps' for comparison with haplotypes.bgen
    '''
    variants = []
    path = Path(__file__).parent / 'data' / 'haplotypes.haps'
    with open(path, 'rt') as haps:
        for line in haps:
            chrom, varid, rsid, pos, ref, alt, *probs = line.strip('\n').split(' ')
            probs = [[1.0, 0.0] if x == '0' else [0.0, 1.0] for x in probs]
            probs = [probs[pos:pos + 2] for pos in range(0, len(probs), 2)]
            probs = [x[0] + x[1] for x in probs]
            var = GenVar(chrom, varid, rsid, pos, [ref, alt], np.array(probs, dtype=np.int8))
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
