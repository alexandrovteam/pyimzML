__version__ = u'1.1.2'

from pyimzml.ImzMLParser import ImzMLParser
from pyimzml.ImzMLWriter import ImzMLWriter
from pyimzml.compression import NoCompression, ZlibCompression

__all__ = ["ImzMLParser", "ImzMLWriter", "NoCompression", "ZlibCompression"]
