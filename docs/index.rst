Welcome to pyimzML documentation!
===================================

.. toctree::
   :maxdepth: 2

This package provides a parser of imzML format as well as a simple imzML writer.

Typical usage pattern is as follows:

.. code-block:: python

    from pyimzml.ImzMLParser import ImzMLParser

    p = ImzMLParser('Example.imzML')
    for idx, (x,y,z) in enumerate(p.coordinates):
        mzs, intensities = p.getspectrum(idx)
        # ...

    from pyimzml.ImzMLWriter import ImzMLWriter

    with ImzMLWriter('output.imzML') as w:
        for mzs, intensities, coords in my_spectra:
            # writes data to the .ibd file
            writer.addSpectrum(mzs, intensities, coords)
    # at this point imzML file is written and files are closed


.. _api:

API Reference
=============

.. automodule:: pyimzml.ImzMLParser
    :members:
    :undoc-members:
    :show-inheritance:

.. automodule:: pyimzml.ImzMLWriter
    :members:
    :undoc-members:
    :show-inheritance:

.. automodule:: pyimzml.compression
    :members:
    :undoc-members:
    :show-inheritance:
