from importlib.metadata import version

__name__ = 'bgen'
__version__ = version(__name__)

from bgen.reader import BgenReader, BgenVar
from bgen.writer import BgenWriter

# listed explicitly so type checkers treat these as re-exported from bgen
__all__ = ['BgenReader', 'BgenVar', 'BgenWriter']
