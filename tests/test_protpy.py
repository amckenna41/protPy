###########################################################################
#################            protPy Module Tests          #################
###########################################################################

import os
import unittest
from Bio import SeqIO
unittest.TestLoader.sortTestMethodsUsing = None

import protpy as protpy

class ProtpyTests(unittest.TestCase):
    """
    Test suite for testing protpy package and its variables and metadata.

    Test Cases
    ----------
    test_protpy_metadata:
        testing pypi software package metadata.
    test_valid_descriptors:
        testing correct list of valid descriptors.
    """
    def test_protpy_metadata(self):
        """ Testing correct protpy version and metadata. """
        self.assertEqual(protpy.__version__, "1.4.2",
            "protpy version is not correct, expected version 1.4.2, got version {}".format(protpy.__version__))
        self.assertEqual(protpy.__name__, "protpy", 
            "protpy software name is not correct, expected protpy, got {}.".format(protpy.__name__))
        self.assertEqual(protpy.__author__, "AJ McKenna: https://github.com/amckenna41", 
            "protpy author is not correct, expected AJ McKenna, got {}.".format(protpy.__author__))
        self.assertEqual(protpy.__authorEmail__, "amckenna41@qub.ac.uk", 
            "protpy author email is not correct, expected amckenna41@qub.ac.uk, got {}.".format(protpy.__authorEmail__))
        self.assertEqual(protpy.__url__, "https://github.com/amckenna41/protPy", 
            "protpy repo URL is not correct, expected https://github.com/amckenna41/protPy, got {}.".format(protpy.__url__))
        self.assertEqual(protpy.__download_url__, "https://github.com/amckenna41/protPy/archive/refs/heads/main.zip", 
            "protpy repo download URL is not correct, expected https://github.com/amckenna41/protPy/archive/refs/heads/main.zip, got {}.".format(protpy.__download_url__))
        self.assertEqual(protpy.__status__, "Production", 
            "protpy status is not correct, expected Production got {}.".format(protpy.__status__))
        self.assertEqual(protpy.__license__, "MIT", 
            "protpy license type is not correct, expected MIT, got {}.".format(protpy.__license__))
        self.assertEqual(protpy.__maintainer__, "AJ McKenna", 
            "protpy maintainer is not correct, expected AJ McKenna, got {}.".format(protpy.__license__))
        self.assertEqual(protpy.__keywords__, ["bioinformatics", "protein engineering", 
            "python", "pypi", "machine learning", "aaindex", "protein descriptors", 
            "physicochemical descriptors", "biochemical descriptors", "structural descriptors", "pySAR"], 
            "protpy keywords is not correct, got\n{}.".format(protpy.__keywords__))

    def test_valid_descriptors(self):
        """ Testing correct list of available descriptors is included in package attribute. """
        self.assertEqual(len(protpy.descriptor_names), 36, 
            "Expected there to be 36 total descriptors, got {}.".format(len(protpy.descriptor_names)))
        self.assertIsInstance(protpy.descriptor_names, list, 
            "all_descriptors should be of type list, got {}.".format(type(protpy.descriptor_names)))
        # composition descriptors
        self.assertIn('amino_acid_composition', protpy.descriptor_names, 
            "amino_acid_composition should be in available descriptors list.")
        self.assertIn('dipeptide_composition', protpy.descriptor_names, 
            "dipeptide_composition should be in available descriptors list.")
        self.assertIn('tripeptide_composition', protpy.descriptor_names,
            "tripeptide_composition should be in available descriptors list.")
        self.assertIn('gravy', protpy.descriptor_names,
            "gravy should be in available descriptors list.")
        self.assertIn('aromaticity', protpy.descriptor_names,
            "aromaticity should be in available descriptors list.")
        self.assertIn('instability_index', protpy.descriptor_names,
            "instability_index should be in available descriptors list.")
        self.assertIn('isoelectric_point', protpy.descriptor_names,
            "isoelectric_point should be in available descriptors list.")
        self.assertIn('molecular_weight', protpy.descriptor_names,
            "molecular_weight should be in available descriptors list.")
        self.assertIn('charge_distribution', protpy.descriptor_names,
            "charge_distribution should be in available descriptors list.")
        self.assertIn('hydrophobic_polar_charged_composition', protpy.descriptor_names,
            "hydrophobic_polar_charged_composition should be in available descriptors list.")
        self.assertIn('secondary_structure_propensity', protpy.descriptor_names,
            "secondary_structure_propensity should be in available descriptors list.")
        self.assertIn('kmer_composition', protpy.descriptor_names,
            "kmer_composition should be in available descriptors list.")
        self.assertIn('reduced_alphabet_composition', protpy.descriptor_names,
            "reduced_alphabet_composition should be in available descriptors list.")
        self.assertIn('motif_composition', protpy.descriptor_names,
            "motif_composition should be in available descriptors list.")
        self.assertIn('amino_acid_pair_composition', protpy.descriptor_names,
            "amino_acid_pair_composition should be in available descriptors list.")
        self.assertIn('aliphatic_index', protpy.descriptor_names,
            "aliphatic_index should be in available descriptors list.")
        self.assertIn('extinction_coefficient', protpy.descriptor_names,
            "extinction_coefficient should be in available descriptors list.")
        self.assertIn('boman_index', protpy.descriptor_names,
            "boman_index should be in available descriptors list.")
        self.assertIn('aggregation_propensity', protpy.descriptor_names,
            "aggregation_propensity should be in available descriptors list.")
        self.assertIn('hydrophobic_moment', protpy.descriptor_names,
            "hydrophobic_moment should be in available descriptors list.")
        self.assertIn('pseudo_amino_acid_composition', protpy.descriptor_names,
            "pseudo_amino_acid_composition should be in available descriptors list.")
        self.assertIn('amphiphilic_pseudo_amino_acid_composition', protpy.descriptor_names, 
            "amphiphilic_pseudo_amino_acid_composition should be in available descriptors list.")
        # autocorrelation descriptors
        self.assertIn('moreaubroto_autocorrelation', protpy.descriptor_names,
            "moreaubroto_autocorrelation should be in available descriptors list.")
        self.assertIn('moran_autocorrelation', protpy.descriptor_names, 
            "moran_autocorrelation should be in available descriptors list.")
        self.assertIn('geary_autocorrelation', protpy.descriptor_names,
            "geary_autocorrelation should be in available descriptors list.")
        # conjoint triad descriptor
        self.assertIn('conjoint_triad', protpy.descriptor_names,
            "conjoint_triad should be in available descriptors list.")
        # CTD descriptors
        self.assertIn('ctd_composition', protpy.descriptor_names,
            "ctd_composition should be in available descriptors list.")
        self.assertIn('ctd_transition', protpy.descriptor_names,
            "ctd_transition should be in available descriptors list.")
        self.assertIn('ctd_distribution', protpy.descriptor_names,
            "ctd_distribution should be in available descriptors list.")
        self.assertIn('ctd_', protpy.descriptor_names,
            "ctd_ should be in available descriptors list.")
        # sequence order descriptors
        self.assertIn('sequence_order_coupling_number_', protpy.descriptor_names,
            "sequence_order_coupling_number_ should be in available descriptors list.")
        self.assertIn('sequence_order_coupling_number', protpy.descriptor_names,
            "sequence_order_coupling_number should be in available descriptors list.")
        self.assertIn('sequence_order_coupling_number_all', protpy.descriptor_names,
            "sequence_order_coupling_number_all should be in available descriptors list.")
        self.assertIn('quasi_sequence_order', protpy.descriptor_names,
            "quasi_sequence_order should be in available descriptors list.")
        self.assertIn('quasi_sequence_order_all', protpy.descriptor_names,
            "quasi_sequence_order_all should be in available descriptors list.")
        self.assertIn('shannon_entropy', protpy.descriptor_names,
            "shannon_entropy should be in available descriptors list.")