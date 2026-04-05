from .autocorrelation import *
from .composition import *
from .conjoint_triad import *
from .ctd import *
from .sequence_order import *
from .utils import get_descriptor_info, list_descriptors, print_descriptor_info

#canonical amino acids available in the package
amino_acids = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P",
    "Q", "R", "S", "T", "V", "W", "Y"]

#protPy package metadata
__name__ = 'protpy'
__version__ = "1.3.0"
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

#list of all available descriptors in protPy
all_descriptors = ["amino_acid_composition", "dipeptide_composition", "tripeptide_composition",
    "gravy", "aromaticity", "instability_index", "isoelectric_point", "molecular_weight",
    "charge_distribution", "hydrophobic_polar_charged_composition", "secondary_structure_propensity",
    "kmer_composition", "reduced_alphabet_composition", "motif_composition",
    "amino_acid_pair_composition",
    "aliphatic_index", "extinction_coefficient", "boman_index",
    "aggregation_propensity", "hydrophobic_moment",
    "moreaubroto_autocorrelation", "moran_autocorrelation", "geary_autocorrelation",
    "pseudo_amino_acid_composition", "amphiphilic_pseudo_amino_acid_composition",
    "conjoint_triad", "ctd_composition", "ctd_transition", "ctd_distribution", "ctd_",
    "sequence_order_coupling_number_", "sequence_order_coupling_number",
    "sequence_order_coupling_number_all", "quasi_sequence_order", "quasi_sequence_order_all",
    "shannon_entropy"
    ]