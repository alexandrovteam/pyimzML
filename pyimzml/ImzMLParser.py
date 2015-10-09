# -*- coding: UTF-8 -*-

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

from bisect import bisect_left

try:
    from lxml.etree import ElementTree
except ImportError:
    try:
        from xml.etree import cElementTree as ElementTree
    except ImportError:
        from xml.etree import ElementTree
import mmap
import struct
from warnings import warn

import numpy as np


class ImzMLParser:
    """
    Parser for imzML 1.1.0 files (see specification here:
    http://imzml.org/download/imzml/specifications_imzML1.1.0_RC1.pdf).

    Reads the .imzML file into memory once and entirely during initialization and
    returns a spectrum upon calling getspectrum(i). The binary file is read in
    every call of getspectrum(i). Use enumerate(parser.coordinates) to get all
    coordinates with their respective index. Relevant meta data is stored in
    parser.imzmldict
    """

    def __init__(self, filename):
        """
        Opens the two files corresponding to the file name, reads the entire .imzML
        file and extracts required attributes. Does not read any binary data, yet.

        :param filename:
            name of the XML file. Must end with .imzML. Binary data file must be named equally but ending with .ibd
        """
        # custom map sizes are currently not supported, therefore mapsize is hardcoded.
        mapsize = 0
        # ElementTree requires the schema location for finding tags (why?) but
        # fails to read it from the root element. As this should be identical
        # for all imzML files, it is hard-coded here and prepended before every tag
        self.sl = "{http://psi.hupo.org/ms/mzml}"
        # maps each imzML number format to its struct equivalent
        self.precisionDict = {"32-bit float": 'f', "64-bit float": 'd', "32-bit integer": 'i', "64-bit integer": 'l'}
        # maps each number format character to its amount of bytes used
        self.sizeDict = {'f': 4, 'd': 8, 'i': 4, 'l': 8}
        self.filename = filename
        self.mzOffsets = []
        self.intensityOffsets = []
        self.mzLengths = []
        self.intensityLengths = []
        # list of all (x,y) or (x,y,z) coordinates as tuples.
        self.coordinates = []
        self.root = ElementTree.parse(filename).getroot()
        # name of the binary file
        bin_filename = self.filename[:-5] + "ibd"
        f = open(bin_filename, "rb")
        # try memory mapping if possible
        try:
            self.m = mmap.mmap(f.fileno(), mapsize, access=mmap.ACCESS_READ)
        # if mmap failed to map file to memory
        except EnvironmentError:
            print ("Memory mapping failed. Using regular file pointer")
            self.m = f
        self.mzPrecision, self.intensityPrecision = self.readformats()
        self.readspectrummeta()
        # Dict for basic imzML metadata other than those required for reading
        # spectra. See method readimzmlmeta()
        self.imzmldict = self.readimzmlmeta()

    # system method for use of 'with ... as'
    def __enter__(self):
        return self

    # system method for use of 'with ... as'
    def __exit__(self, exc_t, exc_v, trace):
        self.m.close()

    def readformats(self):
        """
        This method should only be called by __init__.
        Reads the number formatting strings from the xml tree stored in the name
        attribute of cvParams with accession values "MS:1000521", "MS:1000523",
        "IMS:1000141" or "IMS:1000142". The string values are "32-bit float",
        "64-bit float", "32-bit integer", "64-bit integer".

        Output:

        mz_precision:
            number format for mz array as string specified above
        intensity_precision:
            number format for intensity array as string specified above

        :rtype: Tuple[str, str]
        :raises RuntimeError:
            if the format string could not be found.
        """
        valid_accession_strings = ("MS:1000521", "MS:1000523", "IMS:1000141", "IMS:1000142")

        param_group_list = self.root.find('%sreferenceableParamGroupList' % self.sl)
        param_groups = param_group_list.findall('%sreferenceableParamGroup' % self.sl)
        mz_group = None
        int_group = None
        mz_precision = None
        int_precision = None
        for group in param_groups:
            for param in group:
                if param.attrib["name"] == "m/z array":
                    mz_group = group
                    self.mzGroupId = group.attrib["id"]
                elif param.attrib["name"] == "intensity array":
                    int_group = group
                    self.intGroupId = group.attrib["id"]
        for s in valid_accession_strings:
            param = mz_group.find('%scvParam[@accession="%s"]' % (self.sl, s))
            if param is not None:
                mz_precision = self.precisionDict[param.attrib["name"]]
                break
        for s in valid_accession_strings:
            param = int_group.find('%scvParam[@accession="%s"]' % (self.sl, s))
            if param is not None:
                int_precision = self.precisionDict[param.attrib["name"]]
                break
        if (mz_precision is None) or (int_precision is None):
            raise RuntimeError("Unsupported number format: mz = %s, int = %s" % (mz_precision, int_precision))
        return mz_precision, int_precision

    def readspectrummeta(self):
        """
        This method should only be called by __init__.
        Initializes the attributes mzOffsets, mzLengths, intensityOffsets,
        coordinates and intensityLenghts with their respective values read from the xml tree.
        """
        runelem = self.root.find('%srun' % self.sl)
        speclistelem = runelem.find('%sspectrumList' % self.sl)
        for spectrumElem in list(speclistelem):
            arrlistelem = spectrumElem.find('%sbinaryDataArrayList' % self.sl)
            elist = list(arrlistelem)
            elist_sorted = [None, None]
            for e in elist:
                ref = e.find('%sreferenceableParamGroupRef' % self.sl).attrib["ref"]
                if ref == self.mzGroupId:
                    elist_sorted[0] = e
                elif ref == self.intGroupId:
                    elist_sorted[1] = e
            mz_offset_elem = elist_sorted[0].find('%scvParam[@accession="IMS:1000102"]' % self.sl)
            self.mzOffsets.append(int(mz_offset_elem.attrib["value"]))
            mz_length_elem = elist_sorted[0].find('%scvParam[@accession="IMS:1000103"]' % self.sl)
            self.mzLengths.append(int(mz_length_elem.attrib["value"]))
            intensity_offset_elem = elist_sorted[1].find('%scvParam[@accession="IMS:1000102"]' % self.sl)
            self.intensityOffsets.append(int(intensity_offset_elem.attrib["value"]))
            intensity_length_elem = elist_sorted[1].find('%scvParam[@accession="IMS:1000103"]' % self.sl)
            self.intensityLengths.append(int(intensity_length_elem.attrib["value"]))

            scan_elem = spectrumElem.find('%sscanList/%sscan' % (self.sl, self.sl))
            x = scan_elem.find('%scvParam[@accession="IMS:1000050"]' % self.sl).attrib["value"]
            y = scan_elem.find('%scvParam[@accession="IMS:1000051"]' % self.sl).attrib["value"]
            try:
                z = scan_elem.find('%scvParam[@accession="IMS:1000052"]' % self.sl).attrib["value"]
                self.coordinates.append((int(x), int(y), int(z)))
            except AttributeError:
                self.coordinates.append((int(x), int(y)))

    def readimzmlmeta(self):
        """
        This method should only be called by __init__. Initializes the imzmldict with frequently used metadata from
        the .imzML file.

        This method reads only a subset of the available meta information and may be extended in the future. The keys
        are named similarly to the imzML names. Currently supported keys: "max dimension x", "max dimension y",
        "pixel size x", "pixel size y", "matrix solution concentration", "wavelength", "focus diameter x",
        "focus diameter y", "pulse energy", "pulse duration", "attenuation".

        If a key is not found in the XML tree, it will not be in the dict either.

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
        supportedparams1 = [("max count of pixels x", int), ("max count of pixels y", int),
                            ("max dimension x", int), ("max dimension y", int), ("pixel size x", float),
                            ("pixel size y", float), ("matrix solution concentration", float)]
        supportedparams2 = [("wavelength", float),
                            ("focus diameter x", float), ("focus diameter y", float), ("pulse energy", float),
                            ("pulse duration", float), ("attenuation", float)]
        supportedaccessions1 = [("IMS:1000042", "value"), ("IMS:1000043", "value"),
                                ("IMS:1000044", "value"), ("IMS:1000045", "value"),
                                ("IMS:1000046", "value"), ("IMS:1000047", "value"), ("MS:1000835", "value")]
        supportedaccessions2 = [("MS:1000843", "value"), ("MS:1000844", "value"),
                                ("MS:1000845", "value"), ("MS:1000846", "value"), ("MS:1000847", "value"),
                                ("MS:1000848", "value")]
        for i in range(len(supportedparams1)):
            acc, attr = supportedaccessions1[i]
            elem = scan_settings_list_elem.find('.//%scvParam[@accession="%s"]' % (self.sl, acc))
            if elem is None:
                break
            name, T = supportedparams1[i]
            try:
                d[name] = T(elem.attrib[attr])
            except ValueError:
                warn(Warning('Wrong data type in XML file. Skipped attribute "%s"' % name))
        for i in range(len(supportedparams2)):
            acc, attr = supportedaccessions2[i]
            elem = instrument_config_list_elem.find('.//%scvParam[@accession="%s"]' % (self.sl, acc))
            if elem is None:
                break
            name, T = supportedparams2[i]
            try:
                d[name] = T(elem.attrib[attr])
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
        image_x, image_y = self.coordinates[i]
        return image_x * pixel_size_x, image_y * pixel_size_y

    def getspectrum(self, index):
        """
        Reads the spectrum at specified index from the .ibd file.

        :param index:
            Index of the desired spectrum in the .imzML file

        Output:

        mz_array:
            Sequence of m/z values representing the horizontal axis of the desired mass
            spectrum
        intensity_array:
            Sequence of intensity values corresponding to mz_array
        """
        mz_string, intensity_string = self.get_spectrum_as_string(index)
        mz_fmt = '<' + str(int(len(mz_string) / self.sizeDict[self.mzPrecision])) + self.mzPrecision
        intensity_fmt = '<' + str(
            int(len(intensity_string) / self.sizeDict[self.intensityPrecision])) + self.intensityPrecision
        mz_array = struct.unpack(mz_fmt, mz_string)
        intensity_array = struct.unpack(intensity_fmt, intensity_string)
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


def getionimage(p, mz_value, tol=0.1, z=None, reduce_func=sum):
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
    if z:
        im = np.zeros((p.imzmldict["max count of pixels y"], p.imzmldict["max count of pixels x"], z))
        for i, (x, y, z) in enumerate(p.coordinates):
            mzs, ints = p.getspectrum(i)
            min_i, max_i = _bisect_spectrum(mzs, mz_value, tol)
            im[y - 1, x - 1, z - 1] = reduce_func(ints[min_i:max_i])
    else:
        im = np.zeros((p.imzmldict["max count of pixels y"], p.imzmldict["max count of pixels x"]))
        for i, (x, y,) in enumerate(p.coordinates):
            mzs, ints = p.getspectrum(i)
            min_i, max_i = _bisect_spectrum(mzs, mz_value, tol)
            im[y - 1, x - 1] = reduce_func(ints[min_i:max_i])
    return im


def _bisect_spectrum(mzs, mz_value, tol):
    return bisect_left(mzs, mz_value - tol), bisect_left(mzs, mz_value + tol) + 1
