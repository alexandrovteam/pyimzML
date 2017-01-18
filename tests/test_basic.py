import unittest
from context import getspectrum
import pyimzml.ImzMLParser as imzmlp

class ImzMLParser(unittest.TestCase):
    def bisect_test(self):
        mzs, intensities = getspectrum(200, 1000, 1000)
        test_mzs = [321, 150]
        test_tols = [0.1, 0.1]
        ix_l, ix_u = imzmlp._bisect_spectrum(mzs, test_mz, test_tol)
        assert mzs[ix_l] >= test_mz - test_tol
        assert mzs[ix_u] <= test_mz + test_tol

if __name__ == '__main__':
    unittest.main()