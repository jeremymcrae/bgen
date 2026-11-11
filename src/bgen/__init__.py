from importlib.metadata import version

__version__ = version('bgen')

from bgen.reader import BgenReader, BgenVar
from bgen.writer import BgenWriter

# listed explicitly so type checkers treat these as re-exported from bgen
__all__ = ['BgenReader', 'BgenVar', 'BgenWriter']
