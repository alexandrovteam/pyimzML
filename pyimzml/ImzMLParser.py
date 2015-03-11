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

from bisect import bisect_left, bisect_right
# import ElementTree
from xml.etree import ElementTree
import mmap
import struct
from warnings import warn

"""
Parser for imzML 1.1.0 files (see specification here:
http://imzml.org/download/imzml/specifications_imzML1.1.0_RC1.pdf).

Reads the .imzML file into memory once and entirely during initialization and
returns a spectrum upon calling getspectrum(i). The binary file is read in
every call of getspectrum(i). Use enumerate(parser.coordinates) to get all
coordinates with their respective index. Relevant meta data is stored in
parser.imzmldict
"""
class ImzMLParser:
    """
    Opens the two files corresponding to the file name, reads the entire .imzML
    file and extracts required attributes. Does not read any binary data, yet.

    input
    filename:
        name of the XML file. Must end with .imzML. Binary data file must be
        named equally ending with .ibd
    """
    # def __init__(self, filename, mapsize=0):
    def __init__(self, filename):
        # custom map sizes are currently not supported, therefore mapsize is hardcoded.
        mapsize = 0
        # ElementTree requires the schema location for finding tags (why?) but
        # fails to read it from the root element. As this should be identical
        # for all imzML files, it is hard-coded here and prepended before every tag
        self.sl = "{http://psi.hupo.org/ms/mzml}"
        # maps each imzML number format to its struct equivalent
        self.precisionDict = {"32-bit float" : 'f', "64-bit float" : 'd', "32-bit integer" : 'i', "64-bit integer" : 'l'}
        # maps each number format character to its amount of bytes used
        self.sizeDict = {'f' : 4, 'd' : 8, 'i' : 4, 'l' : 8}
        self.filename = filename
        self.mzOffsets = []
        self.intensityOffsets = []
        self.mzLengths = []
        self.intensityLengths = []
        # list of all (x,y) or (x,y,z) coordinates as tuples.
        self.coordinates = []
        self.root = ElementTree.parse(filename).getroot()
        # name of the binary file
        binFilename = self.filename[:-5] + "ibd"
        f = open(binFilename, "rb")
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
    # system method for use of 'with ... as' (not working yet, do not use)
    def __enter__(self):
        pass
    # system method for use of 'with ... as' (not working yet, do not use)
    def __exit__(self, exc_t, exc_v, trace):
        self.m.close()

    """
    This method should only be called by __init__
    Reads the number formatting strings from the xml tree stored in the name
    attribute of cvParams with accession values "MS:1000521", "MS:1000523",
    "IMS:1000141" or "IMS:1000142". The string values are "32-bit float",
    "64-bit float", "32-bit integer", "64-bit integer".

    output
    mzPrecision:
        number format for mz array as string specified above
    intensityPrecision:
        number format for intensity array as string specified above

    Exceptions
    RuntimeError:
        if the format string could not be found.
    """
    def readformats(self):
        validAccessionStrings = ("MS:1000521", "MS:1000523", "IMS:1000141", "IMS:1000142")

        paramGroupList = self.root.find('%sreferenceableParamGroupList' % (self.sl))
        paramGroups = paramGroupList.findall('%sreferenceableParamGroup' % (self.sl))
        mzGroup = None
        intGroup = None
        mzPrecision = None
        intPrecision = None
        for group in paramGroups:
            for param in group:
                if param.attrib["name"] == "m/z array":
                    mzGroup = group
                    self.mzGroupId = group.attrib["id"]
                elif param.attrib["name"] == "intensity array":
                    intGroup = group
                    self.intGroupId = group.attrib["id"]
        for s in validAccessionStrings:
            param = mzGroup.find('%scvParam[@accession="%s"]' % (self.sl, s))
            if param is not None:
                mzPrecision = self.precisionDict[param.attrib["name"]]
                break
        for s in validAccessionStrings:
            param = intGroup.find('%scvParam[@accession="%s"]' % (self.sl, s))
            if param is not None:
                intPrecision = self.precisionDict[param.attrib["name"]]
                break
        if (mzPrecision is None) or (intPrecision is None):
            raise RuntimeError("Unsupported number format: mz = %s, int = %s" % (mzPrecision, intPrecision))
        return mzPrecision, intPrecision

    """
    This method should only be called by __init__.
    Initializes the attributes mzOffsets, mzLengths, intensityOffsets,
    coordinates and intensityLenghts with their respective values read from the xml tree.
    """
    def readspectrummeta(self):
        runElem = self.root.find('%srun' % (self.sl))
        specListElem = runElem.find('%sspectrumList' % (self.sl))
        for spectrumElem in list(specListElem):
            arrListElem = spectrumElem.find('%sbinaryDataArrayList' % (self.sl))
            eList = list(arrListElem)
            eListSorted = [None, None]
            for e in eList:
                ref = e.find('%sreferenceableParamGroupRef' % (self.sl)).attrib["ref"]
                if ref == self.mzGroupId:
                    eListSorted[0] = e
                elif ref == self.intGroupId:
                    eListSorted[1] = e
            mzOffsetElem = eListSorted[0].find('%scvParam[@accession="IMS:1000102"]' % (self.sl))
            self.mzOffsets.append(int(mzOffsetElem.attrib["value"]))
            mzLengthElem = eListSorted[0].find('%scvParam[@accession="IMS:1000103"]' % (self.sl))
            self.mzLengths.append(int(mzLengthElem.attrib["value"]))
            intensityOffsetElem = eListSorted[1].find('%scvParam[@accession="IMS:1000102"]' % (self.sl))
            self.intensityOffsets.append(int(intensityOffsetElem.attrib["value"]))
            intensityLengthElem = eListSorted[1].find('%scvParam[@accession="IMS:1000103"]' % (self.sl))
            self.intensityLengths.append(int(intensityLengthElem.attrib["value"]))

            scanElem = spectrumElem.find('%sscanList/%sscan' % (self.sl, self.sl))
            x = scanElem.find('%scvParam[@accession="IMS:1000050"]' % self.sl).attrib["value"]
            y = scanElem.find('%scvParam[@accession="IMS:1000051"]' % self.sl).attrib["value"]
            try:
                z = scanElem.find('%scvParam[@accession="IMS:1000052"]' % self.sl).attrib["value"]
                self.coordinates.append((int(x),int(y),int(z)))
            except AttributeError:
                self.coordinates.append((int(x),int(y)))

    """
    This method should only be called by __init__.
    Initializes the imzmldict with frequently used metadata ftom the .imzML file.
    This dict keys are incomplete and may be extended in the future. The keys are named
    similarly to the imzML names. Currently supported keys:
    "max dimension x", "max dimension y", "pixel size x", "pixel size y",
    "matrix solution concentration", "wavelength", "focus diameter x",
    "focus diameter y", "pulse energy", "pulse duration", "attenuation".

    If a key is not found in the XML tree, it will not be in the dict either.

    Output
    d:
        Dict containing above mentioned meta data

    Warnings
    Warning:
        if an xml attribute has a number format different from the imzML specificateion
    """
    def readimzmlmeta(self):
        d = {}
        scanSettingsListElem = self.root.find('%sscanSettingsList' % (self.sl))
        instrumentConfigListElem = self.root.find('%sinstrumentConfigurationList' % (self.sl))
        supportedParams1 = [("max count of pixels x", int), ("max count of pixels y",int),
        ("max dimension x",int), ("max dimension y",int), ("pixel size x",float),
        ("pixel size y",float), ("matrix solution concentration",float)]
        supportedParams2 = [("wavelength",float),
        ("focus diameter x",float), ("focus diameter y",float), ("pulse energy",float),
        ("pulse duration",float), ("attenuation",float)]
        supportedAccessions1 = [("IMS:1000042", "value"), ("IMS:1000043", "value"),
        ("IMS:1000044", "value"), ("IMS:1000045", "value"),
        ("IMS:1000046", "value"), ("IMS:1000047", "value"), ("MS:1000835", "value")]
        supportedAccessions2 = [("MS:1000843", "value"), ("MS:1000844", "value"),
        ("MS:1000845", "value"), ("MS:1000846", "value"), ("MS:1000847", "value"),
        ("MS:1000848", "value")]
        for i in range(len(supportedParams1)):
            acc, attr = supportedAccessions1[i]
            elem = scanSettingsListElem.find('.//%scvParam[@accession="%s"]' % (self.sl, acc))
            if elem is None:
                break
            name, T = supportedParams1[i]
            try:
                d[name] = T(elem.attrib[attr])
            except ValueError:
                warn(Warning('Wrong data type in XML file. Skipped attribute "%s"' % name))
        for i in range(len(supportedParams2)):
            acc, attr = supportedAccessions2[i]
            elem = instrumentConfigListElem.find('.//%scvParam[@accession="%s"]' % (self.sl, acc))
            if elem is None:
                break
            name, T = supportedParams2[i]
            try:
                d[name] = T(elem.attrib[attr])
            except ValueError:
                warn(Warning('Wrong data type in XML file. Skipped attribute "%s"' % name))
        return d

    """
    Reads the spectrum at specified index from the .ibd file.

    input
    index:
        Index of the desired spectrum in the .imzML file

    output
    mzArray:
        List of m/z values representing the horizontal axis of the desired mass
        spectrum
    intensityArray:
        List of intensity values corresponding to mzArray
    """
    def getspectrum(self, index):
        mzString, intensityString = self.get_spectrum_as_string(index)
        mzFmt = '<'+str(int(len(mzString)/self.sizeDict[self.mzPrecision]))+self.mzPrecision
        intensityFmt = '<'+str(int(len(intensityString)/self.sizeDict[self.intensityPrecision]))+self.intensityPrecision
        mzArray = struct.unpack(mzFmt, mzString)
        intensityArray = struct.unpack(intensityFmt, intensityString)
        return mzArray, intensityArray

    """
    Reads m/z array and intensity array of the spectrum at specified location
    from the binary file as a string. The string can be unpacked by the struct
    module. To get the arrays as numbers, use getspectrum

    input
    index:
        Index of the desired spectrum in the .imzML file

    output
    mzString:
        string where each character represents a byte of the mz array of the
        spectrum
    intensityString:
        string where each character represents a byte of the intensity array of
        the spectrum
    """
    def get_spectrum_as_string(self, index):
        offsets = [self.mzOffsets[index], self.intensityOffsets[index]]
        lengths = [self.mzLengths[index], self.intensityLengths[index]]
        lengths[0] = lengths[0]*self.sizeDict[self.mzPrecision]
        lengths[1] = lengths[1]*self.sizeDict[self.intensityPrecision]
        mzString = ""
        intensityString = ""
        self.m.seek(offsets[0])
        mzString = self.m.read(lengths[0])
        self.m.seek(offsets[1])
        intensityString = self.m.read(lengths[1])
        return mzString, intensityString

    """
    Get an image representation of the intensity distribution
    of the ion with specified m/z value.

    input
    mzValue:
        m/z value for which the ion image shall be returned
    tol:
        Absolute tolerance for the m/z value, such that all ions with values
        mzValue-|tol| <= x <= mzValue+|tol| are included. Defaults to 0.1
    z:
        z Value if spectrogram is 3-dimensional.

    output
        numpy matrix with each element representing the ion intensity in this
        pixel. Can be easily plotted with matplotlib
    """
    def getionimage(self, mzValue, tol=0.1, z=None):
        # import numpy here locally, in order to have the rest of the code
        # working if numpy is not installed
        import numpy as np
        tol = abs(tol)
        im = None
        if z:
            im = np.zeros((self.imzmldict["max count of pixels y"], self.imzmldict["max count of pixels x"], z))
            for i, (x,y,z) in enumerate(self.coordinates):
                mzA, intA = self.getspectrum(i)
                minI = bisect_left(mzA, mzValue-tol)
                maxI = bisect_left(mzA, mzValue+tol)+1
                im[y-1,x-1,z-1] = max(intA[minI:maxI])
        else:
            im = np.zeros((self.imzmldict["max count of pixels y"], self.imzmldict["max count of pixels x"]))
            for i, (x,y,) in enumerate(self.coordinates):
                mzA, intA = self.getspectrum(i)
                minI = bisect_left(mzA, mzValue-tol)
                maxI = bisect_left(mzA, mzValue+tol)+1
                im[y-1,x-1] = max(intA[minI:maxI])
        return im
