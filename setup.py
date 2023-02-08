###############################################################################
#####   Setup.py - installs all the required packages and dependancies    #####
###############################################################################

import pathlib
from setuptools import setup, find_packages
import sys
import protpy

#ensure python version is greater than 3
if (sys.version_info[0] < 3):
    sys.exit('Python 3 is the minimum version requirement.')

#get path to README file
HERE = pathlib.Path(__file__).parent
README = (HERE / 'README.md').read_text()

setup(name=protpy.__name__,
      version=protpy.__version__,
      description=protpy.__description__,
      long_description = README,
      long_description_content_type = "text/markdown",
      author=protpy.__license__,
      author_email=protpy.__authorEmail__,
      maintainer=protpy.__maintainer__,
      license=protpy.__license__,
      url=protpy.__url__,
      download_url=protpy.__download_url__,
      keywords=protpy.__keywords__,
      classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Environment :: Console',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Healthcare Industry',
        'Intended Audience :: Information Technology',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3.7',
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
        "varname"
      ],
      test_suite=protpy.__test_suite__,
      # packages=find_packages(), #create Manifest file to ignore results folder in dist
      packages=find_packages(),
      include_package_data=True,
      zip_safe=False)