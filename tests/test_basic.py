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

    def test_parse_metadata(self):
        for parse_lib, data_name, imzml_path, ibd_path in ALL_TEST_CASES:
            with self.subTest(parse_lib=parse_lib, data=data_name),\
                 imzmlp.ImzMLParser(imzml_path, parse_lib=parse_lib) as parser:
                md = parser.metadata
                # fileDescription section
                assert md.file_description['MS:1000579'] == True
                assert 'ibd SHA-1' in md.file_description
                assert len(md.file_description.source_files) == 1
                assert md.file_description.source_files['sf1']['Thermo RAW format'] == True
                assert md.file_description.source_files['sf1'].attrs['name'] == 'Example.raw'
                assert len(md.file_description.contacts) == 1

                # referenceableParamGroupList section
                assert len(md.referenceable_param_groups) == 4
                assert md.referenceable_param_groups['scan1']['increasing m/z scan']

                # sampleList section
                assert len(md.samples) == 1
                assert md.samples['sample1']['sample number'] == '1'

                # softwareList section
                assert len(md.softwares) == 2
                assert md.softwares['Xcalibur']['Xcalibur']

                # scanSettingsList section
                assert len(md.scan_settings) == 1
                assert md.scan_settings['scansettings1']['pixel size (x)'] == 100.0

                # instrumentConfigurationList section
                assert len(md.instrument_configurations) == 1
                ic = md.instrument_configurations['LTQFTUltra0']
                assert ic.param_by_name['instrument serial number'] == 'none'
                assert len(ic.components) == 3
                assert ic.components[0].type == 'source'
                assert ic.components[1].type == 'analyzer'
                assert ic.components[2].type == 'detector'
                assert ic.software_ref == 'Xcalibur'

                # dataProcessingList section
                assert len(md.data_processings) == 2
                assert md.data_processings['XcaliburProcessing'].methods[0].attrs['softwareRef'] == 'Xcalibur'
                assert md.data_processings['XcaliburProcessing'].methods[0]['low intensity data point removal']

    def test_parse_full_spectrum_metadata(self):
        for parse_lib, data_name, imzml_path, ibd_path in ALL_TEST_CASES:
            with self.subTest(parse_lib=parse_lib, data=data_name),\
                 imzmlp.ImzMLParser(imzml_path, parse_lib=parse_lib, include_spectra_metadata='full') as parser:
                assert len(parser.spectrum_full_metadata) == len(parser.coordinates)
                spectrum = parser.spectrum_full_metadata[0]
                assert spectrum['ms level'] == 0  # comes from referenceable param group
                assert spectrum['total ion current'] > 100
                assert spectrum.scan_list_params['no combination']
                assert spectrum.scans[0].attrs['instrumentConfigurationRef'] == 'LTQFTUltra0'
                assert spectrum.scans[0]['position x'] == 1
                assert 'm/z array' in spectrum.binary_data_arrays[0]
                assert 'intensity array' in spectrum.binary_data_arrays[1]

    def test_parse_partial_spectrum_metadata(self):
        TIC, POS_X, EXT_LEN, INVALID = 'MS:1000285', 'IMS:1000050', 'IMS:1000104', 'INVALID'
        ACCESSIONS = [TIC, POS_X, EXT_LEN, INVALID]
        for parse_lib, data_name, imzml_path, ibd_path in ALL_TEST_CASES:
            with self.subTest(parse_lib=parse_lib, data=data_name),\
                 imzmlp.ImzMLParser(imzml_path, parse_lib=parse_lib, include_spectra_metadata=ACCESSIONS) as parser:

                assert len(parser.spectrum_metadata_fields[TIC]) == len(parser.coordinates)
                assert len(parser.spectrum_metadata_fields[POS_X]) == len(parser.coordinates)
                assert len(parser.spectrum_metadata_fields[EXT_LEN]) == len(parser.coordinates)
                assert len(parser.spectrum_metadata_fields[INVALID]) == len(parser.coordinates)

                assert all(tic > 100 for tic in parser.spectrum_metadata_fields[TIC])
                assert all(isinstance(pos_x, int) for pos_x in parser.spectrum_metadata_fields[POS_X])
                assert all(isinstance(ext_len, int) for ext_len in parser.spectrum_metadata_fields[EXT_LEN])
                assert all(invalid is None for invalid in parser.spectrum_metadata_fields[INVALID])



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