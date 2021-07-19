## 1.5.1 (2021-07-??)
* Fix code that causes `SyntaxWarning` in Python 3.8+
* Change `ImzmlWriter` to output "linescan left right" instead of "line scan left right", to match the ontology

## 1.5.0 (2021-07-19)
* Handle mismatched accession for "positive scan"
* Default `ImzMLParser` to ElementTree if no `parse_lib` is specified
* Add `ImzMLParser.polarity` field

## 1.4.1 (2020-10-26)
* Fixed new modules missing from package

## 1.4.0 (2020-10-23)
* Add support for parsing all ImzML metadata
    * Global metadata is always included through `ImzMLParser.metadata`
    * Per-spectrum metadata requires `include_spectra_metadata='full'` 
      or `include_spectra_metadata=[... list of accessions]` to be passed to ImzMLParser.
* Handle mismatched accessions for specifying data types of binary arrays

## 1.3.0 (2019-05-24)
* Add `PortableSpectrumReader`, which holds the minimal subset of `ImzMLParser` needed to read m/z and intensity
  data from the .ibd file, and is able to be pickled. 
  
## 1.2.6 (2019-04-23)
* Changed `ImzMLParser.getspectrum` to return NumPy arrays instead of Python lists

## 1.2.5 (2019-04-10)
* Added `parse_lib` parameter to `ImzMLParser`, allowing ElementTree to be used instead of lxml

## 1.2.4 (2019-01-23)
* Support `MS:1000519` and `MS:1000522` accessions for specifying integer binary data types 

## 1.2.3 (2018-07-02)
* Support `ImzMLParser` detecting .ibd files with a case-insensitive search
