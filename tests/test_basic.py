import unittest
import numpy as np
from context import getspectrum
import pyimzml.ImzMLParser as imzmlp
import pyimzml.ImzMLWriter as imzmlw

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
   		test_mz = 401.0
        test_tol = 0.1
        ix_l, ix_u = imzmlp._bisect_spectrum(mzs, test_mz, test_tol)
        print(ix_l, ix_u, len(mzs))
        assert ix_l == 7
        assert ix_u == 7
        assert ix_l <= ix_u


class ImzMLWriter(unittest.TestCase):
    def test_simple_write(self):
        mzs = np.linspace(100,1000,20)
        ints = np.random.rand(mzs.shape[0])
        coords = [1,1,1]
        with imzmlw.ImzMLWriter("test.mzML", mode="processed") as imzml:
            imzml.addSpectrum(mzs, ints, coords=coords)


    def test_two_spectra_write(self):
        with imzmlw.ImzMLWriter("test.mzML", mode="processed") as imzml:
            mzs = np.linspace(100,1000,20)
            ints = np.random.rand(mzs.shape[0])
            coords = [1,1,1]
            imzml.addSpectrum(mzs, ints, coords=coords)
            ints = np.random.rand(mzs.shape[0])
            coords = [1, 2, 1]
            imzml.addSpectrum(mzs, ints, coords=coords)

    def test_userParam_write(self):
        with imzmlw.ImzMLWriter("test.mzML", mode="processed") as imzml:
            mzs = np.linspace(100,1000,20)
            ints = np.random.rand(mzs.shape[0])
            coords = [1,1,1]
            pos = [50, 100, 0]
            userParams = [{'name': 'xCoord', "value": str(pos[0])}, {'name': 'yCoord', 'value': str(pos[1])}]
            imzml.addSpectrum(mzs, ints, coords=coords, userParams=userParams)
            ints = np.random.rand(mzs.shape[0])
            coords = [1, 2, 1]
            pos = [50, 200, 0]
            userParams = [{'name': 'xCoord', "value": str(pos[0])}, {'name': 'yCoord', 'value': str(pos[1])}]
            imzml.addSpectrum(mzs, ints, coords=coords, userParams=userParams)

if __name__ == '__main__':
    unittest.main()