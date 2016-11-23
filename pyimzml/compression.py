import zlib

class NoCompression(object):
    """
    No compression.
    """
    def __init__(self):
        pass

    def rounding(self, data):
        return data

    def compress(self, bytes):
        return bytes

    def decompress(self, bytes):
        return bytes

    name = "no compression"

class ZlibCompression(object):
    """
    Zlib compression with optional rounding of values.
    Rounding helps the compression, but is lossy.

    :param round_amt:
        Number of digits after comma. None means no rounding.
    """
    def __init__(self, round_amt=None):
        self.round_amt = round_amt

    def rounding(self, data):
        if self.round_amt is not None:
            return [round(x, self.round_amt) for x in data]
        return data

    def compress(self, bytes):
        return zlib.compress(bytes)

    def decompress(self, bytes):
        return zlib.decompress(bytes)

    name = "zlib compression"
