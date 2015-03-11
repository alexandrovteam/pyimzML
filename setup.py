from setuptools import setup

setup(
    name='pyimzML',
    version='1.0.0',
    description="Parser for conversion of imzML 1.1.0 files",
    long_description="""
Parser for conversion of imzML 1.1.0 files.
See specification here: http://imzml.org/download/imzml/specifications_imzML1.1.0_RC1.pdf.
Outputs data as python lists and dicts.

All Python versions from 2.7 on are supported, but Python 3.3 or newer is
recommended for performance reasons.""",
    # The project's main homepage.
    url='https://github.com/pyIMS/pyimzML',
    # Author details
    author='Dominik Fay',
    author_email='dominik.fay@embl.de',

    license='Apache 2.0',
    # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: Apache Software License',
        # 'Programming Language :: Python :: 2',
        # 'Programming Language :: Python :: 2.6', # trouble with ElementTree
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.2',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
    ],
    keywords='bioinformatics imaging mass spectrometry parser imzML',

    # You can just specify the packages manually here if your project is
    # simple. Or you can use find_packages().
    # packages=find_packages(),
    packages=['pyimzml'],

    setup_requires=['numpy'],
    install_requires=['numpy'],
)
