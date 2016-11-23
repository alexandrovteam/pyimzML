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

* pyimzML has an optional dependency to `lxml <http://lxml.de/index.html>`_. If lxml is not installed, pyimzML will instead use the built-in cElementTree or ElementTree package.

**Testing**

To test your installation of pyimzML, you can download sample data from `imzml.org <http://imzml.org/index.php?option=com_content&view=article&id=186&Itemid=68>`_ and go through the following use cases.

Documentation
-------------

Documentation is available on `ReadTheDocs <http://pyimzml.readthedocs.org/en/latest/?badge=latest>`
