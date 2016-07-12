pyimzML
=======

.. image:: https://readthedocs.org/projects/pyimzml/badge/?version=latest
    :target: http://pyimzml.readthedocs.org/en/latest/?badge=latest
    :alt: Documentation Status

Description
-----------
A parser for the imzML format used in imaging mass spectrometry. See specification
`here  <http://imzml.org/download/imzml/specifications_imzML1.1.0_RC1.pdf>`_.
Designed for use with imzML version 1.1.0. Outputs data as python lists and dicts.

The parser is developed by `Alexandrov Team <http://www.embl.de/research/units/scb/alexandrov/index.html>`_ at EMBL Heidelberg.

* Developer: Dominik Fay
* Contact and feedback: `Team contact <http://www.embl.de/research/units/scb/alexandrov/contact/index.html>`_

Installation
------------
pyimzML is available on `PyPI <https://pypi.python.org/pypi/pyimzML>`_. pyimzML
should be installed with pip using one of these three options:

* ``$ pip install pyimzml`` will install pyimzML from PyPI (easiest).
* ``$ pip install git+git://github.com/alexandrovteam/pyimzML.git`` will install pyimzML from github.
* Download the source tarball from `PyPI <https://pypi.python.org/pypi/pyimzML>`_ and ``$ pip install pyimzml-x-x-x.tar.gz``

**Dependency Notes**

* If numpy is not preinstalled, pip will build it from source before installing pyimzML. This can take several minutes.
* pyimzML has an optional dependency to `lxml <http://lxml.de/index.html>`_. If lxml is not installed, pyimzML will instead use the built-in cElementTree or ElementTree package.

**Testing**

To test your installation of pyimzML, you can download sample data from `imzml.org <http://imzml.org/index.php?option=com_content&view=article&id=186&Itemid=68>`_ and go through the following use cases.

Usage
-----
**Init**

First, we instanciate the parser. This will parse the entire .imzML file and
store the xml tree in memory

.. code-block:: python

    from pyimzml.ImzMLParser import ImzMLParser

    p = ImzMLParser('Example.imzML')

Note that the binary file must have the same name as the .imzML file and must
end with \'.ibd\'

**Spectra**

Now we want to get the the entire list of spectra (for later conversion into a
different format, for example). Spectra are accessed by their index in the
.imzML file passed as an argument. Each call of ``getspectrum`` will result in
two reading accesses in the binary file: One for the m/z array, one for the
intensity array. The coordinate of an index can be looked up in the attribute
``coordinates``. As spectral data images are not neccessarily rectangular, the
following generalized approach is recommended

.. code-block:: python

    for i, (x,y,z) in enumerate(p.coordinates):
        p.getspectrum(i)

Each spectrum is a tuple of two lists, the m/z array and the intensity array,
and can be plotted easily

.. code-block:: python

    import matplotlib.pyplot as plt

    mzA, intA = p.getspectrum(1)
    plt.plot(mzA, intA, 'o-')

**Ion images**

Now, we want an image displaying the intensity of an ion throughout the
spectrogram. Therefore, we use the getionimage function, which returns a numpy
matrix that can be displayed as an image by ``matplotlib``

.. code-block:: python

    from pyimzml.ImzMLParser import getionimage

    # pick a base peak
    mzA, intA = p.getspectrum(1)
    peakMz = mzA[intA.index(max(intA))]

    # show the image
    im = getionimage(p, peakMz)
    plt.imshow(im).set_interpolation('nearest')
    plt.colorbar()
    plt.show()

**Metadata**

*pyimzML* has limited support for the metadata embedded in the imzML file. For some general metadata, you can use
the parser's ``ìmzmldict`` attribute. You can find the exact list of supported metadata in the documentation of the
``__readimzmlmeta`` method.

Furthermore, *pyimzML* offers some per-spectrum metadata through its ``browse`` function. Follow the example below:

.. code-block:: python

    # get a list of the instrument configurations used in the first pixel
    instrument_configurations = browse(p).for_spectrum(0).get_ids("instrumentConfiguration")

For browsing all spectra iteratively, you should by all means use **ascending** indices. Doing otherwise can result
in quadratic runtime. The following example shows how to retrieve all unique instrumentConfigurations used:

.. code-block:: python

    browser = browse(p)
    all_config_ids = set()
    for i, _ in enumerate(p.coordinates):
        all_config_ids.update(browser.for_spectrum(i).get_ids("instrumentConfiguration"))

This is a list of ids with which you can find the corresponding ``<instrumentConfiguration>`` tag in the xml tree.
Currently, ``instrumentConfiguration``, ``dataProcessing`` and ``referenceableParamGroup`` are supported.

.. code-block:: python

    # use the ids to find the instrumentConfiguration elements in the ElementTree and get their cvParams
    cv_params = []
    for ref in instrument_configurations:
        # find the corresponding instrumentConfiguration element
        config = p.root.find('%sinstrumentConfigurationList/%sinstrumentConfiguration[@id="%s"]' % (p.sl, p.sl, ref))
        # get its cvParams
        cv_params += [param.attrib for param in config.findall("%scvParam" % p.sl)]
    print cv_params

