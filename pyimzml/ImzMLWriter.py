from __future__ import print_function

import os
import numpy as np
import uuid
import hashlib
import sys
import getopt
from collections import namedtuple, OrderedDict, defaultdict

from wheezy.template import Engine, CoreExtension, DictLoader

from pyimzml.compression import NoCompression, ZlibCompression

IMZML_TEMPLATE = """\
@require(uuid, sha1sum, mz_data_type, int_data_type, run_id, spectra, mode, obo_codes, obo_names, mz_compression, int_compression, polarity, spec_type, scan_direction, scan_pattern, scan_type, line_scan_direction)
<?xml version="1.0" encoding="ISO-8859-1"?>
<mzML xmlns="http://psi.hupo.org/ms/mzml" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://psi.hupo.org/ms/mzml http://psidev.info/files/ms/mzML/xsd/mzML1.1.0_idx.xsd" version="1.1">
  <cvList count="2">
    <cv uri="http://psidev.cvs.sourceforge.net/*checkout*/psidev/psi/psi-ms/mzML/controlledVocabulary/psi-ms.obo" fullName="Proteomics Standards Initiative Mass Spectrometry Ontology" id="MS" version="3.65.0"/>
    <cv uri="http://obo.cvs.sourceforge.net/*checkout*/obo/obo/ontology/phenotype/unit.obo" fullName="Unit Ontology" id="UO" version="12:10:2011"/>
  </cvList>

  <fileDescription>
    <fileContent>
      <cvParam cvRef="MS" accession="MS:1000579" name="MS1 spectrum" value=""/>
      @if spec_type=='centroid':
      <cvParam cvRef="MS" accession="MS:1000127" name="centroid spectrum" value=""/>
      @elif spec_type=='profile':
      <cvParam cvRef="MS" accession="MS:1000128" name="profile spectrum" value=""/>
      @end
      <cvParam cvRef="IMS" accession="IMS:@obo_codes[mode]" name="@mode" value=""/>
      <cvParam cvRef="IMS" accession="IMS:1000080" name="universally unique identifier" value="@uuid"/>
      <cvParam cvRef="IMS" accession="IMS:1000091" name="ibd SHA-1" value="@sha1sum"/>
    </fileContent>
  </fileDescription>

  <referenceableParamGroupList count="4">
    <referenceableParamGroup id="mzArray">
      <cvParam cvRef="MS" accession="MS:@obo_codes[mz_compression]" name="@mz_compression" value=""/>
      <cvParam cvRef="MS" accession="MS:1000514" name="m/z array" unitCvRef="MS" unitAccession="MS:1000040" unitName="m/z"/>
      <cvParam cvRef="MS" accession="MS:@obo_codes[mz_data_type]" name="@mz_data_type" value=""/>
      <cvParam cvRef="IMS" accession="IMS:1000101" name="external data" value="true"/>
    </referenceableParamGroup>
    <referenceableParamGroup id="intensityArray">
      <cvParam cvRef="MS" accession="MS:@obo_codes[int_data_type]" name="@int_data_type" value=""/>
      <cvParam cvRef="MS" accession="MS:1000515" name="intensity array" unitCvRef="MS" unitAccession="MS:1000131" unitName="number of detector counts"/>
      <cvParam cvRef="MS" accession="MS:@obo_codes[int_compression]" name="@int_compression" value=""/>
      <cvParam cvRef="IMS" accession="IMS:1000101" name="external data" value="true"/>
    </referenceableParamGroup>
    <referenceableParamGroup id="scan1">
      <cvParam cvRef="MS" accession="MS:1000093" name="increasing m/z scan"/>
      <cvParam cvRef="MS" accession="MS:1000512" name="filter string" value=""/>
    </referenceableParamGroup>
    <referenceableParamGroup id="spectrum1">
      <cvParam cvRef="MS" accession="MS:1000579" name="MS1 spectrum" value=""/>
      <cvParam cvRef="MS" accession="MS:1000511" name="ms level" value="0"/>
      @if spec_type=='centroid':
      <cvParam cvRef="MS" accession="MS:1000127" name="centroid spectrum" value=""/>
      @elif spec_type=='profile':
      <cvParam cvRef="MS" accession="MS:1000128" name="profile spectrum" value=""/>
      @end
      @if polarity=='positive':
      <cvParam cvRef="MS" accession="MS:1000130" name="positive scan" value=""/>
      @elif polarity=='negative':
      <cvParam cvRef="MS" accession="MS:1000129" name="negative scan" value=""/>
      @end
    </referenceableParamGroup>
  </referenceableParamGroupList>

  <softwareList count="1">
    <software id="pyimzml" version="0.0001">
      <cvParam cvRef="MS" accession="MS:1000799" name="custom unreleased software tool" value="pyimzml exporter"/>
    </software>
  </softwareList>

  <scanSettingsList count="1">
    <scanSettings id="scanSettings1">
      <cvParam cvRef="IMS" accession="IMS:@obo_codes[scan_direction]" name="@obo_names[scan_direction]"/>
      <cvParam cvRef="IMS" accession="IMS:@obo_codes[scan_pattern]" name="@obo_names[scan_pattern]"/>
      <cvParam cvRef="IMS" accession="IMS:@obo_codes[scan_type]" name="@obo_names[scan_type]"/>
      <cvParam cvRef="IMS" accession="IMS:@obo_codes[line_scan_direction]" name="@obo_names[line_scan_direction]"/>
      <cvParam cvRef="IMS" accession="IMS:1000042" name="max count of pixels x" value="@{(max(s.coords[0] for s in spectra))!!s}"/>
      <cvParam cvRef="IMS" accession="IMS:1000043" name="max count of pixels y" value="@{(max(s.coords[1] for s in spectra))!!s}"/>
    </scanSettings>
  </scanSettingsList>

  <instrumentConfigurationList count="1">
    <instrumentConfiguration id="IC1">
    </instrumentConfiguration>
  </instrumentConfigurationList>

  <dataProcessingList count="1">
    <dataProcessing id="export_from_pyimzml">
      <processingMethod order="0" softwareRef="pyimzml">
        <cvParam cvRef="MS" accession="MS:1000530" name="file format conversion" value="Output to imzML"/>
      </processingMethod>
    </dataProcessing>
  </dataProcessingList>

  <run defaultInstrumentConfigurationRef="IC1" id="@run_id">
    <spectrumList count="@{len(spectra)!!s}" defaultDataProcessingRef="export_from_pyimzml">
      @for index, s in enumerate(spectra):
      <spectrum defaultArrayLength="0" id="spectrum=@{(index+1)!!s}" index="@{(index+1)!!s}">
        <referenceableParamGroupRef ref="spectrum1"/>
        <cvParam cvRef="MS" accession="MS:1000528" name="lowest observed m/z" value="@{s.mz_min!!s}" unitCvRef="MS" unitAccession="MS:1000040" unitName="m/z"/>
        <cvParam cvRef="MS" accession="MS:1000527" name="highest observed m/z" value="@{s.mz_max!!s}" unitCvRef="MS" unitAccession="MS:1000040" unitName="m/z"/>
        <cvParam cvRef="MS" accession="MS:1000504" name="base peak m/z" value="@{s.mz_base!!s}" unitCvRef="MS" unitAccession="MS:1000040" unitName="m/z"/>
        <cvParam cvRef="MS" accession="MS:1000505" name="base peak intensity" value="@{s.int_base!!s}" unitCvRef="MS" unitAccession="MS:1000131" unitName="number of counts"/>
        <cvParam cvRef="MS" accession="MS:1000285" name="total ion current" value="@{s.int_tic!!s}"/>
        <scanList count="1">
          <cvParam accession="MS:1000795" cvRef="MS" name="no combination"/>
          <scan instrumentConfigurationRef="instrumentConfiguration0">
            <referenceableParamGroupRef ref="scan1"/>
            <cvParam accession="IMS:1000050" cvRef="IMS" name="position x" value="@{s.coords[0]!!s}"/>
            <cvParam accession="IMS:1000051" cvRef="IMS" name="position y" value="@{s.coords[1]!!s}"/>
            @if len(s.coords) == 3:
            <cvParam accession="IMS:1000052" cvRef="IMS" name="position z" value="@{s.coords[2]!!s}"/>
            @end
            @if s.userParams:
                @for up in s.userParams:
                <userParam name="@up['name']" value="@up['value']"/> 
                @end
            @end
          </scan>
        </scanList>
        <binaryDataArrayList count="2">
          <binaryDataArray encodedLength="0">
            <referenceableParamGroupRef ref="mzArray"/>
            <cvParam accession="IMS:1000103" cvRef="IMS" name="external array length" value="@{s.mz_len!!s}"/>
            <cvParam accession="IMS:1000104" cvRef="IMS" name="external encoded length" value="@{s.mz_enc_len!!s}"/>
            <cvParam accession="IMS:1000102" cvRef="IMS" name="external offset" value="@{s.mz_offset!!s}"/>
            <binary/>
          </binaryDataArray>
          <binaryDataArray encodedLength="0">
            <referenceableParamGroupRef ref="intensityArray"/>
            <cvParam accession="IMS:1000103" cvRef="IMS" name="external array length" value="@{s.int_len!!s}"/>
            <cvParam accession="IMS:1000104" cvRef="IMS" name="external encoded length" value="@{s.int_enc_len!!s}"/>
            <cvParam accession="IMS:1000102" cvRef="IMS" name="external offset" value="@{s.int_offset!!s}"/>
            <binary/>
          </binaryDataArray>
        </binaryDataArrayList>
      </spectrum>
      @end
    </spectrumList>
  </run>
</mzML>
"""

class _MaxlenDict(OrderedDict):
    def __init__(self, *args, **kwargs):
        self.maxlen = kwargs.pop('maxlen', None)
        OrderedDict.__init__(self, *args, **kwargs)

    def __setitem__(self, key, value):
        if self.maxlen is not None and len(self) >= self.maxlen:
            self.popitem(0) #pop oldest
        OrderedDict.__setitem__(self, key, value)

_Spectrum = namedtuple('_Spectrum', 'coords mz_len mz_offset mz_enc_len int_len int_offset int_enc_len mz_min mz_max mz_base int_base int_tic userParams') #todo: change named tuple to dict and parse xml template properly (i.e. remove hardcoding so parameters can be optional)

class ImzMLWriter(object):
    """
        Create an imzML+ibd file.

        :param output_filename:
            is used to make the base name by removing the extension (if any).
            two files will be made by adding ".ibd" and ".imzML" to the base name
        :param intensity_dtype:
            The numpy data type to use for saving intensity values
        :param mz_dtype:
            The numpy data type to use for saving mz array values
        :param mode:

            * "continuous" mode will save the first mz array only
            * "processed" mode save every mz array separately
            * "auto" mode writes only mz arrays that have not already been written
        :param intensity_compression:
            How to compress the intensity data before saving
            must be an instance of :class:`~pyimzml.compression.NoCompression` or :class:`~pyimzml.compression.ZlibCompression`
        :param mz_compression:
            How to compress the mz array data before saving
    """
    def __init__(self, output_filename,
                 mz_dtype=np.float64, intensity_dtype=np.float32, mode="auto", spec_type="centroid",
                 scan_direction="top_down", line_scan_direction="line_left_right", scan_pattern="one_way", scan_type="horizontal_line", 
                 mz_compression=NoCompression(), intensity_compression=NoCompression(),
                 polarity=None):

        self.mz_dtype = mz_dtype
        self.intensity_dtype = intensity_dtype
        self.mode = mode
        self.spec_type = spec_type
        self.mz_compression = mz_compression
        self.intensity_compression = intensity_compression
        self.run_id = os.path.splitext(output_filename)[0]
        self.filename = self.run_id + ".imzML"
        self.ibd_filename = self.run_id + ".ibd"
        self.xml = open(self.filename, 'w')
        self.ibd = open(self.ibd_filename, 'wb+')
        self.sha1 = hashlib.sha1()
        self.uuid = uuid.uuid4()
        
        self.scan_direction = scan_direction
        self.scan_pattern = scan_pattern
        self.scan_type = scan_type
        self.line_scan_direction = line_scan_direction

        self._write_ibd(self.uuid.bytes)

        self.wheezy_engine = Engine(loader=DictLoader({'imzml': IMZML_TEMPLATE}), extensions=[CoreExtension()])
        self.imzml_template = self.wheezy_engine.get_template('imzml')
        self.spectra = []
        self.first_mz = None
        self.hashes = defaultdict(list)  # mz_hash -> list of mz_data (disk location)
        self.lru_cache = _MaxlenDict(maxlen=10)  # mz_array (as tuple) -> mz_data (disk location)
        self._setPolarity(polarity)

    @staticmethod
    def _np_type_to_name(dtype):
        if dtype.__name__.startswith('float'):
            return "%s-bit float" % dtype.__name__[5:]
        elif dtype.__name__.startswith('int'):
            return "%s-bit integer" % dtype.__name__[3:]

    def _setPolarity(self, polarity):
        if polarity:
            if polarity.lower() in ['positive', 'negative']:
                self.polarity = polarity.lower()
            else:
                raise ValueError('value for polarity must be one of "positive", "negative". Received: {}'.format(polarity))
        else:
            self.polarity = ""

    def _write_xml(self):
        spectra = self.spectra
        mz_data_type = self._np_type_to_name(self.mz_dtype)
        int_data_type = self._np_type_to_name(self.intensity_dtype)
        obo_codes = {"32-bit integer": "1000519", 
                     "16-bit float": "1000520",
                     "32-bit float": "1000521",
                     "64-bit integer": "1000522",
                     "64-bit float": "1000523",
                     "continuous": "1000030",
                     "processed": "1000031",
                     "zlib compression": "1000574",
                     "no compression": "1000576",
                     "line_bottom_up": "1000492",
                     "line_left_right": "1000491",
                     "line_right_left": "1000490",
                     "line_top_down": "1000493",
                     "bottom_up": "1000400",
                     "left_right": "1000402",
                     "right_left": "1000403",
                     "top_down": "1000401",
                     "meandering": "1000410",
                     "one_way": "1000411",
                     "random_access": "1000412",
                     "horizontal_line": "1000480",
                     "vertical_line": "1000481"}
        obo_names = {"line_bottom_up": "linescan bottom up",
                     "line_left_right": "linescan left right",
                     "line_right_left": "linescan right left",
                     "line_top_down": "linescan top down",
                     "bottom_up": "bottom up",
                     "left_right": "left right",
                     "right_left": "right left",
                     "top_down": "top down",
                     "meandering": "meandering",
                     "one_way": "one way",
                     "random_access": "random access",
                     "horizontal_line": "horizontal line scan",
                     "vertical_line": "vertical line scan"}
        
        uuid = ("{%s}" % self.uuid).upper()
        sha1sum = self.sha1.hexdigest().upper()
        run_id = self.run_id
        if self.mode == 'auto':
            mode = "processed" if len(self.lru_cache) > 1 else "continuous"
        else:
            mode = self.mode
        spec_type = self.spec_type
        mz_compression = self.mz_compression.name
        int_compression = self.intensity_compression.name
        polarity = self.polarity
        scan_direction = self.scan_direction
        scan_pattern = self.scan_pattern
        scan_type = self.scan_type
        line_scan_direction = self.line_scan_direction
        
        self.xml.write(self.imzml_template.render(locals()))

    def _write_ibd(self, bytes):
        self.ibd.write(bytes)
        self.sha1.update(bytes)
        return len(bytes)

    def _encode_and_write(self, data, dtype=np.float32, compression=NoCompression()):
        data = np.asarray(data, dtype=dtype)
        offset = self.ibd.tell()
        bytes = data.tobytes()
        bytes = compression.compress(bytes)
        return offset, data.shape[0], self._write_ibd(bytes)

    def _read_mz(self, mz_offset, mz_len, mz_enc_len):
        '''reads a mz array from the currently open ibd file'''
        self.ibd.seek(mz_offset)
        data = self.ibd.read(mz_enc_len)
        self.ibd.seek(0, 2)
        data = self.mz_compression.decompress(data)
        return tuple(np.fromstring(data, dtype=self.mz_dtype))

    def _get_previous_mz(self, mzs):
        '''given an mz array, return the mz_data (disk location)
        if the mz array was not previously written, write to disk first'''
        mzs = tuple(mzs)  # must be hashable
        if mzs in self.lru_cache:
            return self.lru_cache[mzs]

        # mz not recognized ... check hash
        mz_hash = "%s-%s-%s" % (hash(mzs), sum(mzs), len(mzs))
        if mz_hash in self.hashes:
            for mz_data in self.hashes[mz_hash]:
                test_mz = self._read_mz(*mz_data)
                if mzs == test_mz:
                    self.lru_cache[test_mz] = mz_data
                    return mz_data
        # hash not recognized
        # must be a new mz array ... write it, add it to lru_cache and hashes
        mz_data = self._encode_and_write(mzs, self.mz_dtype, self.mz_compression)
        self.hashes[mz_hash].append(mz_data)
        self.lru_cache[mzs] = mz_data
        return mz_data

    def addSpectrum(self, mzs, intensities, coords, userParams=[]):
        """
        Add a mass spectrum to the file.

        :param mz:
            mz array
        :param intensities:
            intensity array
        :param coords:

            * 2-tuple of x and y position OR
            * 3-tuple of x, y, and z position

            note some applications want coords to be 1-indexed
        """
        # must be rounded now to allow comparisons to later data
        # but don't waste CPU time in continuous mode since the data will not be used anyway
        if self.mode != "continuous" or self.first_mz is None:
            mzs = self.mz_compression.rounding(mzs)
        intensities = self.intensity_compression.rounding(intensities)

        if self.mode == "continuous":
            if self.first_mz is None:
                self.first_mz = self._encode_and_write(mzs, self.mz_dtype, self.mz_compression)
            mz_data = self.first_mz
        elif self.mode == "processed":
            mz_data = self._encode_and_write(mzs, self.mz_dtype, self.mz_compression)
        elif self.mode == "auto":
            mz_data = self._get_previous_mz(mzs)
        else:
            raise TypeError("Unknown mode: %s" % self.mode)
        mz_offset, mz_len, mz_enc_len = mz_data

        int_offset, int_len, int_enc_len = self._encode_and_write(intensities, self.intensity_dtype, self.intensity_compression)
        mz_min = np.min(mzs)
        mz_max = np.max(mzs)
        ix_max = np.argmax(intensities)
        mz_base = mzs[ix_max]
        int_base = intensities[ix_max]
        int_tic = np.sum(intensities)
        s = _Spectrum(coords, mz_len, mz_offset, mz_enc_len, int_len, int_offset, int_enc_len, mz_min, mz_max, mz_base, int_base, int_tic, userParams)
        self.spectra.append(s)

    def close(self):  # 'close' is a more common use for this
        """
        Writes the XML file and closes all files.
        Will be called automatically if ``with``-pattern is used.
        """
        self.finish()

    def finish(self):
        '''alias of close()'''
        self.ibd.close()
        self._write_xml()
        self.xml.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_t, exc_v, trace):
        if exc_t is None:
            self.finish()
        else:
            self.ibd.close()
            self.xml.close()

def _main(argv):
    from pyimzml.ImzMLParser import ImzMLParser
    inputfile = ''
    outputfile = ''
    try:
        opts, args = getopt.getopt(argv,"hi:o:",["ifile=","ofile="])
    except getopt.GetoptError:
        print('test.py -i <inputfile> -o <outputfile>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('test.py -i <inputfile> -o <outputfile>')
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputfile = arg
        elif opt in ("-o", "--ofile"):
            outputfile = arg
    if inputfile == '':
        print('test.py -i <inputfile> -o <outputfile>')
        raise IOError('input file not specified')
    if outputfile=='':
        outputfile=inputfile+'.imzML'
    imzml = ImzMLParser(inputfile)
    spectra = []
    with ImzMLWriter(outputfile, mz_dtype=np.float32, intensity_dtype=np.float32) as writer:
        for i, coords in enumerate(imzml.coordinates):
            mzs, intensities = imzml.getspectrum(i)
            writer.addSpectrum(mzs, intensities, coords)
            spectra.append((mzs, intensities, coords))

    imzml = ImzMLParser(outputfile)
    spectra2 = []
    for i, coords in enumerate(imzml.coordinates):
        mzs, intensities = imzml.getspectrum(i)
        spectra2.append((mzs, intensities, coords))

    print(spectra[0] == spectra2[0])

if __name__ == '__main__':
    _main(sys.argv[1:])
