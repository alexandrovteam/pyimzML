from lxml.objectify import ElementMaker
import lxml.etree
import os
import numpy as np
import uuid
import hashlib

from wheezy.template.engine import Engine
from wheezy.template.ext.core import CoreExtension
from wheezy.template.loader import DictLoader

IMZML_TEMPLATE = """\
@require(uuid, sha1sum, mz_data_type, int_data_type, run_id, spectra)
<?xml version="1.0" encoding="ISO-8859-1"?>
<mzML xmlns="http://psi.hupo.org/ms/mzml" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://psi.hupo.org/ms/mzml http://psidev.info/files/ms/mzML/xsd/mzML1.1.0_idx.xsd" version="1.1">
  <cvList count="2">
    <cv uri="http://psidev.cvs.sourceforge.net/*checkout*/psidev/psi/psi-ms/mzML/controlledVocabulary/psi-ms.obo" fullName="Proteomics Standards Initiative Mass Spectrometry Ontology" id="MS" version="3.65.0"/>
    <cv uri="http://obo.cvs.sourceforge.net/*checkout*/obo/obo/ontology/phenotype/unit.obo" fullName="Unit Ontology" id="UO" version="12:10:2011"/>
  </cvList>

  <fileDescription>
    <fileContent>
      <cvParam cvRef="MS" accession="MS:1000579" name="MS1 spectrum" value=""/>
      <cvParam cvRef="IMS" accession="IMS:1000031" name="processed" value=""/>
      <cvParam cvRef="IMS" accession="IMS:1000080" name="universally unique identifier" value="@uuid"/>
      <cvParam cvRef="IMS" accession="IMS:1000091" name="ibd SHA-1" value="@sha1sum"/>
    </fileContent>
  </fileDescription>

  <referenceableParamGroupList count="2">
    <referenceableParamGroup id="mzArray">
      <cvParam cvRef="MS" accession="MS:1000514" name="m/z array" value=""/>
      <cvParam cvRef="MS" accession="MS:1000523" name="@mz_data_type" value=""/>
      <cvParam cvRef="MS" accession="MS:1000576" name="no compression" value=""/>
    </referenceableParamGroup>
    <referenceableParamGroup id="intensityArray">
      <cvParam cvRef="MS" accession="MS:1000515" name="intensity array" value=""/>
      <cvParam cvRef="MS" accession="MS:1000523" name="@int_data_type" value=""/>
      <cvParam cvRef="MS" accession="MS:1000576" name="no compression" value=""/>
    </referenceableParamGroup>
  </referenceableParamGroupList>

  <softwareList count="1">
    <software id="pyimzml" version="0.0001">
      <cvParam cvRef="MS" accession="MS:xxx" name="pyimzml" value="pyimzml exporter"/>
    </software>
  </softwareList>
  <scanSettingsList count="1">
    <scanSettings id="scanSettings1">
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
        <referenceableParamGroupRef ref="spectrum"/>
        <scanList count="1">
          <cvParam accession="MS:1000795" cvRef="MS" name="no combination"/>
          <scan instrumentConfigurationRef="instrumentConfiguration0">
            <cvParam accession="IMS:1000050" cvRef="IMS" name="position x" value="@{s.coords[0]!!s}"/>
            <cvParam accession="IMS:1000051" cvRef="IMS" name="position y" value="@{s.coords[1]!!s}"/>
            @if len(s.coords) == 3:
            <cvParam accession="IMS:1000052" cvRef="IMS" name="position z" value="@{s.coords[2]!!s}"/>
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

class ImzMLWriter(object):
    def __init__(self, output_filename, mz_dtype=np.float64, intensity_dtype=np.float32):
        self.mz_dtype = mz_dtype
        self.intensity_dtype = intensity_dtype
        self.run_id = os.path.splitext(output_filename)[0]
        self.filename =  self.run_id + ".imzML"
        self.ibd_filename = self.run_id + ".ibd"
        self.xml = open(self.filename, 'w')
        self.ibd = open(self.ibd_filename, 'w')
        self.sha1 = hashlib.sha1()
        self.ibd_offset = 0
        self.uuid = uuid.uuid4()

        self._write_ibd(self.uuid.bytes_le)

        self.wheezy_engine = Engine(loader=DictLoader({'imzml': IMZML_TEMPLATE}), extensions=[CoreExtension()])
        self.imzml_template = self.wheezy_engine.get_template('imzml')

        self.spectra = []

        from collections import namedtuple
        self.Spectrum = namedtuple('Spectrum', ['coords', 'mz_len', 'mz_offset', 'mz_enc_len', 'int_len', 'int_offset', 'int_enc_len'])

    def _write_xml(self):
        spectra = self.spectra
        mz_data_type = "%d-bit float" % (np.dtype(self.mz_dtype).itemsize * 8)
        int_data_type = "%d-bit float" % (np.dtype(self.intensity_dtype).itemsize * 8)
        uuid = '{' + self.uuid.urn.split(':')[-1].upper() + '}'
        sha1sum = self.sha1.hexdigest().upper()
        run_id = self.run_id
    
        self.xml.write(self.imzml_template.render(locals()))

    def addSpectrum(self, mzs, intensities, coords):
        mzs, intensities = np.asarray(mzs, dtype=self.mz_dtype), np.asarray(intensities, dtype=self.intensity_dtype)

        mz_offset = self.ibd_offset
        mz_len = mzs.shape[0]
        mzs_bytes = mzs.tobytes()
        self._write_ibd(mzs_bytes)
        int_offset = self.ibd_offset
        mz_enc_len = int_offset - mz_offset
        int_len = intensities.shape[0]
        intensities_bytes = intensities.tobytes()
        self._write_ibd(intensities_bytes)
        int_enc_len = self.ibd_offset - int_offset

        s = self.Spectrum(mz_offset=mz_offset, mz_len=mz_len,
                          int_offset=int_offset, int_len=int_len, coords=coords,
                          int_enc_len=int_enc_len, mz_enc_len=mz_enc_len)

        self.spectra.append(s)

    def _write_ibd(self, bytes):
        self.ibd.write(bytes)
        self.sha1.update(bytes)
        self.ibd_offset += len(bytes)

    def finish(self):
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

def main():
    import sys
    sys.path.append("/home/lomereiter/github/PyImzML")
    from pyimzml.ImzMLParser import ImzMLParser
    imzml = ImzMLParser("/home/lomereiter/Downloads/S042_Processed.imzML")
    spectra = []
    with ImzMLWriter("/tmp/test.imzML", mz_dtype=np.float32, intensity_dtype=np.float32) as writer:
        for i, coords in enumerate(imzml.coordinates):
            mzs, intensities = imzml.getspectrum(i)
            writer.addSpectrum(mzs, intensities, coords)
            spectra.append((mzs, intensities, coords))

    imzml = ImzMLParser("/tmp/test.imzML")
    spectra2 = []
    for i, coords in enumerate(imzml.coordinates):
        mzs, intensities = imzml.getspectrum(i)
        spectra2.append((mzs, intensities, coords))

    print spectra[0] == spectra2[0]

if __name__ == '__main__':
    main()
