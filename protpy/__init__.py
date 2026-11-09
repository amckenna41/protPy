from .autocorrelation import *
from .composition import *
from .conjoint_triad import *
from .ctd import *
from .sequence_order import *
from .scales import *
from .physicochemical import *
from .complexity import *
from .spectral import *
from .structural import *
from .utils import get_descriptor_info, list_descriptors, print_descriptor_info
from .encode import all_descriptors, encode_dataset, read_sequences_from_file

#canonical amino acids available in the package
from ._constants import amino_acids

#protPy package metadata
__name__ = 'protpy'
__version__ = "1.6.0"
__description__ = "Python package for generating various biochemical, physicochemical and structural descriptors/features of protein sequences."
__author__ = 'AJ McKenna: https://github.com/amckenna41'
__authorEmail__ = 'amckenna41@qub.ac.uk'
__maintainer__ = "AJ McKenna"
__license__ = 'MIT'
__url__ = 'https://github.com/amckenna41/protPy'
__download_url__ = "https://github.com/amckenna41/protPy/archive/refs/heads/main.zip"
__status__ = "Production"
__keywords__ = ["bioinformatics", "protein engineering", "python", "pypi", "machine learning", \
                "aaindex", "protein descriptors", "physicochemical descriptors", "biochemical descriptors",
                "structural descriptors", "pySAR"]
__test_suite__ = "tests"

#list of all registered descriptor names in protPy (deprecated functions are not
#registered; see protpy.utils._DESCRIPTOR_REGISTRY)
from .utils import _DESCRIPTOR_REGISTRY as _registry
descriptor_names = list(_registry)
del _registry
