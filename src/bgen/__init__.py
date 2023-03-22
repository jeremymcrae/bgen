from pkg_resources import get_distribution

__version__ = get_distribution('bgen').version

from bgen.reader import BgenReader
from bgen.writer import BgenWriter
