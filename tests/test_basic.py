import pickle
import unittest

import numpy as np
from pathlib import Path
from .context import getspectrum
import pyimzml.ImzMLParser as imzmlp
import pyimzml.ImzMLWriter as imzmlw

# Example files from https://ms-imaging.org/wp/imzml/example-files-test/
CONTINUOUS_IMZML_PATH = str(Path(__file__).parent / 'data/Example_Continuous.imzML')
CONTINUOUS_IBD_PATH = str(Path(__file__).parent / 'data/Example_Continuous.ibd')
PROCESSED_IMZML_PATH = str(Path(__file__).parent / 'data/Example_Processed.imzML')
PROCESSED_IBD_PATH = str(Path(__file__).parent / 'data/Example_Processed.ibd')
PARSE_LIB_TEST_CASES = ['lxml', 'ElementTree']
DATA_TEST_CASES = [
    ('Continuous', CONTINUOUS_IMZML_PATH, CONTINUOUS_IBD_PATH),
    ('Processed', PROCESSED_IMZML_PATH, PROCESSED_IBD_PATH),
]
ALL_TEST_CASES = [(parse_lib, data_name, imzml_path, ibd_path)
                  for parse_lib in PARSE_LIB_TEST_CASES
                  for data_name, imzml_path, ibd_path in DATA_TEST_CASES]


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

    def test_getspectrum(self):
        for parse_lib, data_name, imzml_path, ibd_path in ALL_TEST_CASES:
            with self.subTest(parse_lib=parse_lib, data=data_name),\
                 imzmlp.ImzMLParser(imzml_path, parse_lib=parse_lib) as parser:

                mzs, ints = parser.getspectrum(4)

                assert len(parser.coordinates) == 9
                assert mzs.dtype == np.float32
                assert ints.dtype == np.float32
                assert len(mzs) == 8399
                assert len(ints) == 8399
                assert np.all(mzs > 100.0)
                assert np.all(mzs < 800.0)
                assert np.all(ints >= 0.0)
                assert np.all(ints < 3.0)

    def test_files_instead_of_paths(self):
        for parse_lib, data_name, imzml_path, ibd_path in ALL_TEST_CASES:
            with self.subTest(parse_lib=parse_lib, data=data_name),\
                 open(imzml_path, 'rb') as imzml_file,\
                 open(ibd_path, 'rb') as ibd_file,\
                 imzmlp.ImzMLParser(imzml_file, parse_lib=parse_lib, ibd_file=ibd_file) as parser:

                mzs, ints = parser.getspectrum(4)

                assert len(parser.coordinates) == 9
                assert len(mzs) > 0
                assert len(ints) > 0


class PortableSpectrumReader(unittest.TestCase):
    def test_read_file(self):
        spectrum_idx = 4
        for parse_lib, data_name, imzml_path, ibd_path in ALL_TEST_CASES:
            with self.subTest(parse_lib=parse_lib, data=data_name),\
                 imzmlp.ImzMLParser(imzml_path, parse_lib=parse_lib) as normal_parser,\
                 open(ibd_path, 'rb') as ibd_file:

                normal_mzs, normal_ints = normal_parser.getspectrum(spectrum_idx)

                detached_parser = imzmlp.ImzMLParser(imzml_path, parse_lib=parse_lib, ibd_file=None)
                portable_reader = detached_parser.portable_spectrum_reader()
                # Pickle and unpickle to ensure it survives for its intended use case
                portable_reader = pickle.loads(pickle.dumps(portable_reader))
                portable_mzs, portable_ints = portable_reader.read_spectrum_from_file(ibd_file, spectrum_idx)

                assert np.all(normal_mzs == portable_mzs)
                assert np.all(normal_ints == portable_ints)


class ImzMLWriter(unittest.TestCase):
    def test_simple_write(self):
        mzs = np.linspace(100,1000,20)
        ints = np.random.rand(mzs.shape[0])
        coords = [1,1,1]
        with imzmlw.ImzMLWriter("test.mzML", mode="processed") as imzml:
            imzml.addSpectrum(mzs, ints, coords=coords)

if __name__ == '__main__':
    unittest.main()