# -*- coding: utf-8 -*-

# Copyright 2015 Dominik Fay
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from bisect import bisect_left, bisect_right
import sys
import re
from pathlib import Path

from warnings import warn
import numpy as np

from pyimzml.metadata import Metadata, SpectrumData
from pyimzml.ontology.ontology import convert_cv_param

PRECISION_DICT = {"32-bit float": 'f', "64-bit float": 'd', "32-bit integer": 'i', "64-bit integer": 'l'}
SIZE_DICT = {'f': 4, 'd': 8, 'i': 4, 'l': 8}
INFER_IBD_FROM_IMZML = object()
XMLNS_PREFIX = "{http://psi.hupo.org/ms/mzml}"

param_group_elname = "referenceableParamGroup"
data_processing_elname = "dataProcessing"
instrument_confid_elname = "instrumentConfiguration"


def choose_iterparse(parse_lib=None):
    if parse_lib == 'ElementTree':
        from xml.etree.ElementTree import iterparse
    elif parse_lib == 'lxml':
        from lxml.etree import iterparse
    else:
        from xml.etree.ElementTree import iterparse
    return iterparse


def _get_cv_param(elem, accession, deep=False, convert=False):
    base = './/' if deep else ''
    node = elem.find('%s%scvParam[@accession="%s"]' % (base, XMLNS_PREFIX, accession))
    if node is not None:
        if convert:
            return convert_cv_param(accession, node.get('value'))
        return node.get('value')


class ImzMLParser:
    """
    Parser for imzML 1.1.0 files (see specification here:
    http://imzml.org/download/imzml/specifications_imzML1.1.0_RC1.pdf).

    Iteratively reads the .imzML file into memory while pruning the per-spectrum metadata (everything in
    <spectrumList> elements) during initialization. Returns a spectrum upon calling getspectrum(i). The binary file
    is read in every call of getspectrum(i). Use enumerate(parser.coordinates) to get all coordinates with their
    respective index. Coordinates are always 3-dimensional. If the third spatial dimension is not present in
    the data, it will be set to zero.

    The global metadata fields in the imzML file are stored in parser.metadata.
    Spectrum-specific metadata fields are not stored by default due to avoid memory issues,
    use the `include_spectra_metadata` parameter if spectrum-specific metadata is needed.
    """

    def __init__(
            self,
            filename,
            parse_lib=None,
            ibd_file=INFER_IBD_FROM_IMZML,
            include_spectra_metadata=None,
    ):
        """
        Opens the two files corresponding to the file name, reads the entire .imzML
        file and extracts required attributes. Does not read any binary data, yet.

        :param filename:
            name of the XML file. Must end with .imzML. Binary data file must be named equally but ending with .ibd
            Alternatively an open file or Buffer Protocol object can be supplied, if ibd_file is also supplied
        :param parse_lib:
            XML-parsing library to use: 'ElementTree' or 'lxml', the later will be used if argument not provided
        :param ibd_file:
            File or Buffer Protocol object for the .ibd file. Leave blank to infer it from the imzml filename.
            Set to None if no data from the .ibd file is needed (getspectrum calls will not work)
        :param include_spectra_metadata:
            None, 'full', or a list/set of accession IDs.
            If 'full' is given, parser.spectrum_full_metadata will be populated with a list of
                complex objects containing the full metadata for each spectrum.
            If a list or set is given, parser.spectrum_metadata_fields will be populated with a dict mapping
                accession IDs to lists. Each list will contain the values for that accession ID for
                each spectrum. Note that for performance reasons, this mode only searches the
                spectrum itself for the value. It won't check any referenced referenceable param
                groups if the accession ID isn't present in the spectrum metadata.
        """
        # ElementTree requires the schema location for finding tags (why?) but
        # fails to read it from the root element. As this should be identical
        # for all imzML files, it is hard-coded here and prepended before every tag
        self.sl = "{http://psi.hupo.org/ms/mzml}"
        # maps each imzML number format to its struct equivalent
        self.precisionDict = dict(PRECISION_DICT)
        # maps each number format character to its amount of bytes used
        self.sizeDict = dict(SIZE_DICT)
        self.filename = filename
        self.mzOffsets = []
        self.intensityOffsets = []
        self.mzLengths = []
        self.intensityLengths = []
        # list of all (x,y,z) coordinates as tuples.
        self.coordinates = []
        self.root = None
        self.metadata = None
        self.polarity = None
        if include_spectra_metadata == 'full':
            self.spectrum_full_metadata = []
        elif include_spectra_metadata is not None:
            include_spectra_metadata = set(include_spectra_metadata)
            self.spectrum_metadata_fields = {
                k: [] for k in include_spectra_metadata
            }

        self.mzGroupId = self.intGroupId = self.mzPrecision = self.intensityPrecision = None
        self.iterparse = choose_iterparse(parse_lib)
        self.__iter_read_spectrum_meta(include_spectra_metadata)
        if ibd_file is INFER_IBD_FROM_IMZML:
            # name of the binary file
            ibd_filename = self._infer_bin_filename(self.filename)
            self.m = open(ibd_filename, "rb")
        else:
            self.m = ibd_file

        # Dict for basic imzML metadata other than those required for reading
        # spectra. See method __readimzmlmeta()
        self.imzmldict = self.__readimzmlmeta()
        self.imzmldict['max count of pixels z'] = np.asarray(self.coordinates)[:,2].max()

    @staticmethod
    def _infer_bin_filename(imzml_path):
        imzml_path = Path(imzml_path)
        ibd_path = [f for f in imzml_path.parent.glob('*')
                    if re.match(r'.+\.ibd', str(f), re.IGNORECASE) and f.stem == imzml_path.stem][0]
        return str(ibd_path)

    # system method for use of 'with ... as'
    def __enter__(self):
        return self

    # system method for use of 'with ... as'
    def __exit__(self, exc_t, exc_v, trace):
        if self.m is not None:
            self.m.close()

    def __iter_read_spectrum_meta(self, include_spectra_metadata):
        """
        This method should only be called by __init__. Reads the data formats, coordinates and offsets from
        the .imzML file and initializes the respective attributes. While traversing the XML tree, the per-spectrum
        metadata is pruned, i.e. the <spectrumList> element(s) are left behind empty.

        Supported accession values for the number formats: "MS:1000521", "MS:1000523", "IMS:1000141" or
        "IMS:1000142". The string values are "32-bit float", "64-bit float", "32-bit integer", "64-bit integer".
        """
        mz_group = int_group = None
        slist = None
        elem_iterator = self.iterparse(self.filename, events=("start", "end"))

        if sys.version_info > (3,):
            _, self.root = next(elem_iterator)
        else:
            _, self.root = elem_iterator.next()

        is_first_spectrum = True

        for event, elem in elem_iterator:
            if elem.tag == self.sl + "spectrumList" and event == "start":
                self.__process_metadata()
                slist = elem
            elif elem.tag == self.sl + "spectrum" and event == "end":
                self.__process_spectrum(elem, include_spectra_metadata)
                if is_first_spectrum:
                    self.__read_polarity(elem)
                    is_first_spectrum = False
                slist.remove(elem)
        self.__fix_offsets()

    def __fix_offsets(self):
        # clean up the mess after morons who use signed 32-bit where unsigned 64-bit is appropriate
        def fix(array):
            fixed = []
            delta = 0
            prev_value = float('nan')
            for value in array:
                if value < 0 and prev_value >= 0:
                    delta += 2**32
                fixed.append(value + delta)
                prev_value = value
            return fixed

        self.mzOffsets = fix(self.mzOffsets)
        self.intensityOffsets = fix(self.intensityOffsets)

    def __process_metadata(self):
        if self.metadata is None:
            self.metadata = Metadata(self.root)
            for param_id, param_group in self.metadata.referenceable_param_groups.items():
                if 'm/z array' in param_group.param_by_name:
                    self.mzGroupId = param_id
                    for name, dtype in self.precisionDict.items():
                        if name in param_group.param_by_name:
                            self.mzPrecision = dtype
                if 'intensity array' in param_group.param_by_name:
                    self.intGroupId = param_id
                    for name, dtype in self.precisionDict.items():
                        if name in param_group.param_by_name:
                            self.intensityPrecision = dtype
            if not hasattr(self, 'mzPrecision'):
                raise RuntimeError("Could not determine m/z precision")
            if not hasattr(self, 'intensityPrecision'):
                raise RuntimeError("Could not determine intensity precision")

    def __process_spectrum(self, elem, include_spectra_metadata):
        arrlistelem = elem.find('%sbinaryDataArrayList' % self.sl)
        mz_group = None
        int_group = None
        for e in arrlistelem:
            ref = e.find('%sreferenceableParamGroupRef' % self.sl).attrib["ref"]
            if ref == self.mzGroupId:
                mz_group = e
            elif ref == self.intGroupId:
                int_group = e
        self.mzOffsets.append(int(_get_cv_param(mz_group, 'IMS:1000102')))
        self.mzLengths.append(int(_get_cv_param(mz_group, 'IMS:1000103')))
        self.intensityOffsets.append(int(_get_cv_param(int_group, 'IMS:1000102')))
        self.intensityLengths.append(int(_get_cv_param(int_group, 'IMS:1000103')))
        scan_elem = elem.find('%sscanList/%sscan' % (self.sl, self.sl))
        x = _get_cv_param(scan_elem, 'IMS:1000050')
        y = _get_cv_param(scan_elem, 'IMS:1000051')
        z = _get_cv_param(scan_elem, 'IMS:1000052')
        if z is not None:
            self.coordinates.append((int(x), int(y), int(z)))
        else:
            self.coordinates.append((int(x), int(y), 1))

        if include_spectra_metadata == 'full':
            self.spectrum_full_metadata.append(
                SpectrumData(elem, self.metadata.referenceable_param_groups)
            )
        elif include_spectra_metadata:
            for param in include_spectra_metadata:
                value = _get_cv_param(elem, param, deep=True, convert=True)
                self.spectrum_metadata_fields[param].append(value)

    def __read_polarity(self, elem):
        # It's too slow to always check all spectra, so first check the referenceable_param_groups
        # in the header to see if they indicate the polarity. If not, try to detect it from
        # the first spectrum's full metadata.
        # LIMITATION: This won't detect "mixed" polarity if polarity is only specified outside the
        # referenceable_param_groups.
        param_groups = self.metadata.referenceable_param_groups.values()
        spectrum_metadata = SpectrumData(elem, self.metadata.referenceable_param_groups)
        has_positive = (
            any('positive scan' in group for group in param_groups)
            or 'positive scan' in spectrum_metadata
        )
        has_negative = (
            any('negative scan' in group for group in param_groups)
            or 'negative scan' in spectrum_metadata
        )
        if has_positive and has_negative:
            self.polarity = 'mixed'
        elif has_positive:
            self.polarity = 'positive'
        elif has_negative:
            self.polarity = 'negative'

    def __readimzmlmeta(self):
        """
        DEPRECATED - use self.metadata instead, as it has much greater detail and allows for
        multiple scan settings / instruments.

        This method should only be called by __init__. Initializes the imzmldict with frequently used metadata from
        the .imzML file.

        :return d:
            dict containing above mentioned meta data
        :rtype:
            dict
        :raises Warning:
            if an xml attribute has a number format different from the imzML specification
        """
        d = {}
        scan_settings_list_elem = self.root.find('%sscanSettingsList' % self.sl)
        instrument_config_list_elem = self.root.find('%sinstrumentConfigurationList' % self.sl)
        scan_settings_params = [
            ("max count of pixels x", "IMS:1000042"),
            ("max count of pixels y", "IMS:1000043"),
            ("max dimension x", "IMS:1000044"),
            ("max dimension y", "IMS:1000045"),
            ("pixel size x", "IMS:1000046"),
            ("pixel size y", "IMS:1000047"),
            ("matrix solution concentration", "MS:1000835"),
        ]
        instrument_config_params = [
            ("wavelength", "MS:1000843"),
            ("focus diameter x", "MS:1000844"),
            ("focus diameter y", "MS:1000845"),
            ("pulse energy", "MS:1000846"),
            ("pulse duration", "MS:1000847"),
            ("attenuation", "MS:1000848"),
        ]

        for name, accession in scan_settings_params:
            try:
                val = _get_cv_param(scan_settings_list_elem, accession, deep=True, convert=True)
                if val is not None:
                    d[name] = val
            except ValueError:
                warn(Warning('Wrong data type in XML file. Skipped attribute "%s"' % name))

        for name, accession in instrument_config_params:
            try:
                val = _get_cv_param(instrument_config_list_elem, accession, deep=True, convert=True)
                if val is not None:
                    d[name] = val
            except ValueError:
                warn(Warning('Wrong data type in XML file. Skipped attribute "%s"' % name))
        return d

    def get_physical_coordinates(self, i):
        """
        For a pixel index i, return the real-world coordinates in nanometers.

        This is equivalent to multiplying the image coordinates of the given pixel with the pixel size.

        :param i: the pixel index
        :return: a tuple of x and y coordinates.
        :rtype: Tuple[float]
        :raises KeyError: if the .imzML file does not specify the attributes "pixel size x" and "pixel size y"
        """
        try:
            pixel_size_x = self.imzmldict["pixel size x"]
            pixel_size_y = self.imzmldict["pixel size y"]
        except KeyError:
            raise KeyError("Could not find all pixel size attributes in imzML file")
        image_x, image_y = self.coordinates[i][:2]
        return image_x * pixel_size_x, image_y * pixel_size_y

    def getspectrum(self, index):
        """
        Reads the spectrum at specified index from the .ibd file.

        :param index:
            Index of the desired spectrum in the .imzML file

        Output:

        mz_array: numpy.ndarray
            Sequence of m/z values representing the horizontal axis of the desired mass
            spectrum
        intensity_array: numpy.ndarray
            Sequence of intensity values corresponding to mz_array
        """
        mz_bytes, intensity_bytes = self.get_spectrum_as_string(index)
        mz_array = np.frombuffer(mz_bytes, dtype=self.mzPrecision)
        intensity_array = np.frombuffer(intensity_bytes, dtype=self.intensityPrecision)
        return mz_array, intensity_array

    def get_spectrum_as_string(self, index):
        """
        Reads m/z array and intensity array of the spectrum at specified location
        from the binary file as a byte string. The string can be unpacked by the struct
        module. To get the arrays as numbers, use getspectrum

        :param index:
            Index of the desired spectrum in the .imzML file
        :rtype: Tuple[str, str]

        Output:

        mz_string:
            string where each character represents a byte of the mz array of the
            spectrum
        intensity_string:
            string where each character represents a byte of the intensity array of
            the spectrum
        """
        offsets = [self.mzOffsets[index], self.intensityOffsets[index]]
        lengths = [self.mzLengths[index], self.intensityLengths[index]]
        lengths[0] *= self.sizeDict[self.mzPrecision]
        lengths[1] *= self.sizeDict[self.intensityPrecision]
        self.m.seek(offsets[0])
        mz_string = self.m.read(lengths[0])
        self.m.seek(offsets[1])
        intensity_string = self.m.read(lengths[1])
        return mz_string, intensity_string

    def portable_spectrum_reader(self):
        """
        Builds a PortableSpectrumReader that holds the coordinates list and spectrum offsets in the .ibd file
        so that the .ibd file can be read without opening the .imzML file again.

        The PortableSpectrumReader can be safely pickled and unpickled, making it useful for reading the spectra
        in a distributed environment such as PySpark or PyWren.
        """
        return PortableSpectrumReader(self.coordinates,
                                      self.mzPrecision, self.mzOffsets, self.mzLengths,
                                      self.intensityPrecision, self.intensityOffsets, self.intensityLengths)


def getionimage(p, mz_value, tol=0.1, z=1, reduce_func=sum):
    """
    Get an image representation of the intensity distribution
    of the ion with specified m/z value.

    By default, the intensity values within the tolerance region are summed.

    :param p:
        the ImzMLParser (or anything else with similar attributes) for the desired dataset
    :param mz_value:
        m/z value for which the ion image shall be returned
    :param tol:
        Absolute tolerance for the m/z value, such that all ions with values
        mz_value-|tol| <= x <= mz_value+|tol| are included. Defaults to 0.1
    :param z:
        z Value if spectrogram is 3-dimensional.
    :param reduce_func:
        the bahaviour for reducing the intensities between mz_value-|tol| and mz_value+|tol| to a single value. Must
        be a function that takes a sequence as input and outputs a number. By default, the values are summed.

    :return:
        numpy matrix with each element representing the ion intensity in this
        pixel. Can be easily plotted with matplotlib
    """
    tol = abs(tol)
    im = np.zeros((p.imzmldict["max count of pixels y"], p.imzmldict["max count of pixels x"]))
    for i, (x, y, z_) in enumerate(p.coordinates):
        if z_ == 0:
            UserWarning("z coordinate = 0 present, if you're getting blank images set getionimage(.., .., z=0)")
        if z_ == z:
            mzs, ints = map(lambda x: np.asarray(x), p.getspectrum(i))
            min_i, max_i = _bisect_spectrum(mzs, mz_value, tol)
            im[y - 1, x - 1] = reduce_func(ints[min_i:max_i+1])
    return im


def browse(p):
    """
    Create a per-spectrum metadata browser for the parser.
    Usage::

        # get a list of the instrument configurations used in the first pixel
        instrument_configurations = browse(p).for_spectrum(0).get_ids("instrumentConfiguration")

    Currently, ``instrumentConfiguration``, ``dataProcessing`` and ``referenceableParamGroup`` are supported.

    For browsing all spectra iteratively, you should by all means use **ascending** indices. Doing otherwise can result
    in quadratic runtime. The following example shows how to retrieve all unique instrumentConfigurations used::

        browser = browse(p)
        all_config_ids = set()
        for i, _ in enumerate(p.coordinates):
            all_config_ids.update(browser.for_spectrum(i).get_ids("instrumentConfiguration"))

    This is a list of ids with which you can find the corresponding ``<instrumentConfiguration>`` tag in the xml tree.

    :param p: the parser
    :return: the browser
    """
    return _ImzMLMetaDataBrowser(p.root, p.filename, p.sl)


def _bisect_spectrum(mzs, mz_value, tol):
    ix_l, ix_u = bisect_left(mzs, mz_value - tol), bisect_right(mzs, mz_value + tol) - 1
    if ix_l == len(mzs):
        return len(mzs), len(mzs)
    if ix_u < 1:
        return 0, 0
    if ix_u == len(mzs):
        ix_u -= 1
    if mzs[ix_l] < (mz_value - tol):
        ix_l += 1
    if mzs[ix_u] > (mz_value + tol):
        ix_u -= 1
    return ix_l, ix_u


class _ImzMLMetaDataBrowser(object):
    def __init__(self, root, fn, sl):
        self._root = root
        self._sl = sl
        self._fn = fn
        self._iter, self._previous, self._list_elem = None, None, None
        self.iterparse = choose_iterparse()

    def for_spectrum(self, i):
        if self._previous is None or i <= self._previous:
            self._iter = self.iterparse(self._fn, events=("start", "end"))
        for event, s in self._iter:
            if s.tag == self._sl + "spectrumList" and event == "start":
                self._list_elem = s
            elif s.tag == self._sl + "spectrum" and event == "end":
                self._list_elem.remove(s)
                if s.attrib["index"] == str(i):
                    self._previous = i
                    return _SpectrumMetaDataBrowser(self._root, self._sl, s)


class _SpectrumMetaDataBrowser(object):
    def __init__(self, root, sl, spectrum):
        self._root = root
        self._sl = sl
        self._spectrum = spectrum

    def get_ids(self, element):
        param_methods = {
            param_group_elname: self._find_referenceable_param_groups,
            data_processing_elname: self._find_data_processing,
            instrument_confid_elname: self._find_instrument_configurations,
        }
        try:
            return param_methods[element]()
        except KeyError as e:
            raise ValueError("Unsupported element: " + str(element))

    def _find_referenceable_param_groups(self):
        param_group_refs = self._spectrum.findall("%sreferenceableParamGroupRef" % self._sl)
        ids = map(lambda g: g.attrib["ref"], param_group_refs)
        return ids

    def _find_instrument_configurations(self):
        ids = None
        scan_list = self._spectrum.find("%sscanList" % self._sl)
        if scan_list:
            scans = scan_list.findall("%sscan[@instrumentConfigurationRef]" % self._sl)
            ids = map(lambda s: s.attrib["instrumentConfigurationRef"], scans)
        if not ids:
            run = self._root.find("%srun")
            try:
                return [run.attrib["defaultInstrumentConfigurationRef"]]
            except KeyError as _:
                return list()
        else:
            return ids

    def _find_data_processing(self):
        try:
            return self._spectrum.attrib["dataProcessingRef"]
        except KeyError as _:
            spectrum_list = self._root.find("%srun/%sspectrumList" % tuple(2 * [self._sl]))
            try:
                return [spectrum_list.attrib["defaultDataProcessingRef"]]
            except KeyError as _:
                return []


class PortableSpectrumReader(object):
    """
    A pickle-able class for holding the minimal set of data required for reading,
    without holding any references to open files that wouldn't survive pickling.
    """

    def __init__(self, coordinates, mzPrecision, mzOffsets, mzLengths,
                 intensityPrecision, intensityOffsets, intensityLengths):
        self.coordinates = coordinates
        self.mzPrecision = mzPrecision
        self.mzOffsets = mzOffsets
        self.mzLengths = mzLengths
        self.intensityPrecision = intensityPrecision
        self.intensityOffsets = intensityOffsets
        self.intensityLengths = intensityLengths

    def read_spectrum_from_file(self, file, index):
        """
        Reads the spectrum at specified index from the .ibd file.

        :param file:
            File or file-like object for the .ibd file
        :param index:
            Index of the desired spectrum in the .imzML file

        Output:

        mz_array: numpy.ndarray
            Sequence of m/z values representing the horizontal axis of the desired mass
            spectrum
        intensity_array: numpy.ndarray
            Sequence of intensity values corresponding to mz_array
        """
        file.seek(self.mzOffsets[index])
        mz_bytes = file.read(self.mzLengths[index] * SIZE_DICT[self.mzPrecision])
        file.seek(self.intensityOffsets[index])
        intensity_bytes = file.read(self.intensityLengths[index] * SIZE_DICT[self.intensityPrecision])

        mz_array = np.frombuffer(mz_bytes, dtype=self.mzPrecision)
        intensity_array = np.frombuffer(intensity_bytes, dtype=self.intensityPrecision)

        return mz_array, intensity_array
