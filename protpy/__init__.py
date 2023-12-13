from .autocorrelation import *
from .composition import *
from .conjoint_triad import *
from .ctd import *
from .sequence_order import *

#protPy package metadata
__name__ = 'protpy'
__version__ = "1.2.1"
__description__ = "Python package for generating various biochemical, physicochemical and structural descriptors/features of protein sequences."
__author__ = 'AJ McKenna: https://github.com/amckenna41'
__authorEmail__ = 'amckenna41@qub.ac.uk'
__maintainer__ = "AJ McKenna"
__license__ = 'MIT'
__url__ = 'https://github.com/amckenna41/protPy'
__download_url__ = "https://github.com/amckenna41/protPy/archive/refs/heads/main.zip"
__status__ = "Production"
__keywords__ = ["bioinformatics", "protein engineering", "python", "pypi", "machine learning", \
                "aaindex", "protein descriptors", "physicochemical descriptors", "biochemical descriptors"
                "structural descriptors", "pySAR"]
__test_suite__ = "tests"

#list of all available descriptors in protPy
all_descriptors = ["amino_acid_composition", "dipeptide_composition", "tripeptide_composition",
    "moreaubroto_autocorrelation", "moran_autocorrelation", "geary_autocorrelation",
    "pseudo_amino_acid_composition", "amphiphilic_pseudo_amino_acid_composition", 
    "sequence_order_coupling_number", "conjoint_triad", "ctd_composition", 
    "ctd_transition", "ctd_distribution", "quasi_sequence_order"
    ]