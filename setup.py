###############################################################################
#####   Setup.py - installs all the required packages and dependancies    #####
###############################################################################

import pathlib
from setuptools import setup, find_packages

#get path to README file
HERE = pathlib.Path(__file__).parent
README = (HERE / 'README.md').read_text()

#protPy package metadata
__name__ = 'protpy'
__version__ = "1.2.0"
__description__ = "Python package for generating various biochemical, physiochemical and structural descriptors/features of protein sequences."
__author__ = 'AJ McKenna, https://github.com/amckenna41'
__authorEmail__ = 'amckenna41@qub.ac.uk'
__maintainer__ = "AJ McKenna"
__license__ = 'MIT'
__url__ = 'https://github.com/amckenna41/protPy'
__download_url__ = "https://github.com/amckenna41/protPy/archive/refs/heads/main.zip"
__status__ = "Production"
__keywords__ = ["bioinformatics", "protein engineering", "python", "pypi", "machine learning", \
                "aaindex", "protein descriptors", "physiochemical descriptors", "biochemical descriptors"
                "structural descriptors", "pySAR"]
__test_suite__ = "tests"

setup(name=__name__,
      version=__version__,
      description=__description__,
      long_description = README,
      long_description_content_type = "text/markdown",
      author=__license__,
      author_email=__authorEmail__,
      maintainer=__maintainer__,
      license=__license__,
      url=__url__,
      download_url=__download_url__,
      keywords=__keywords__,
      classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Environment :: Console',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Healthcare Industry',
        'Intended Audience :: Information Technology',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: Implementation :: PyPy',
        'Programming Language :: Python :: 3 :: Only',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Medical Science Apps.',
        'Topic :: Scientific/Engineering :: Mathematics',
        'Topic :: Software Development :: Libraries :: Python Modules'
      ],
      install_requires=[
        "aaindex",
        "numpy",
        "pandas",
        "varname",
      ],
      test_suite=__test_suite__,
      # packages=find_packages(), #create Manifest file to ignore results folder in dist
      packages=find_packages(),
      include_package_data=True,
      zip_safe=False)