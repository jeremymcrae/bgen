from importlib.metadata import version

__name__ = 'bgen'
__version__ = version(__name__)

from bgen.reader import BgenReader, BgenVar
from bgen.writer import BgenWriter
