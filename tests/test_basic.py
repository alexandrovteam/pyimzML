import unittest
from context import getspectrum
import pyimzml.ImzMLParser as imzmlp

class ImzMLParser(unittest.TestCase):
    def test_bisect(self):
        mzs = [100., 201.89, 201.99, 202.0, 202.01, 202.10000001, 400.]
        test_mz = 202.0
        test_tol = 0.1
        ix_l, ix_u = imzmlp._bisect_spectrum(mzs, test_mz, test_tol)
        assert ix_l == 2
        assert ix_u == 4
        assert ix_l <= ix_u
        assert mzs[ix_l] >= test_mz - test_tol
        assert mzs[ix_u] <= test_mz + test_tol

if __name__ == '__main__':
    unittest.main()