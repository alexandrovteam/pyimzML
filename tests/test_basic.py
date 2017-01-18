import unittest
from context import getspectrum
import pyimzml.ImzMLParser as imzmlp

class ImzMLParser(unittest.TestCase):
    def test_bisect(self):
        mzs, intensities = getspectrum(200, 1000, 1000)
        test_mz = 321.0
        test_tol = 0.1
        ix_l, ix_u = imzmlp._bisect_spectrum(mzs, test_mz, test_tol)
        assert mzs[ix_l] >= test_mz - test_tol
        assert mzs[ix_u] <= test_mz + test_tol

if __name__ == '__main__':
    unittest.main()