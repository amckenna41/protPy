###########################################################################
#################            protPy Module Tests          #################
###########################################################################

import os
import unittest
from Bio import SeqIO
unittest.TestLoader.sortTestMethodsUsing = None

import protpy as protpy

class ProtPyTests(unittest.TestCase):
    """
    Test suite for testing protpy package. 

    Test Cases
    ----------
    test_protpy_metadata:
        testing pypi software package metadata.
    test_valid_descriptors:
        testing correct list of valid descriptors.
    """
    def setUp(self):
        """ Import protein sequences from test fasta files using Biopython package. """
        #using next() to get first item (protein seq) from SeqIO Generator
        with open(os.path.join("tests", "test_fasta1.fasta")) as pro:
            self.protein_seq1 = str(next(SeqIO.parse(pro,'fasta')).seq)

        with open(os.path.join("tests", "test_fasta2.fasta")) as pro:
            self.protein_seq2 = str(next(SeqIO.parse(pro,'fasta')).seq)
        
        with open(os.path.join("tests", "test_fasta3.fasta")) as pro:
            self.protein_seq3 = str(next(SeqIO.parse(pro,'fasta')).seq)

        with open(os.path.join("tests", "test_fasta4.fasta")) as pro:
            self.protein_seq4 = str(next(SeqIO.parse(pro,'fasta')).seq)

        #list of canonical amino acids
        self.amino_acids = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", 
            "Q", "R", "S", "T", "V", "W", "Y"]

    def test_protpy_metadata(self):
        """ Testing correct protpy version and metadata. """
        self.assertEqual(protpy.__version__, "1.0.4", 
            "protpy version is not correct, got: {}".format(protpy.__version__))
        self.assertEqual(protpy.__name__, "protpy", 
            "protpy software name is not correct, got: {}".format(protpy.__name__))
        self.assertEqual(protpy.__author__, "AJ McKenna, https://github.com/amckenna41", 
            "protpy author is not correct, got: {}".format(protpy.__author__))
        self.assertEqual(protpy.__authorEmail__, "amckenna41@qub.ac.uk", 
            "protpy author email is not correct, got: {}".format(protpy.__authorEmail__))
        self.assertEqual(protpy.__url__, "https://github.com/amckenna41/protPy", 
            "protpy repo URL is not correct, got: {}".format(protpy.__url__))
        self.assertEqual(protpy.__download_url__, "https://github.com/amckenna41/protPy/archive/refs/heads/main.zip", 
            "protpy repo download URL is not correct, got: {}".format(protpy.__download_url__))
        self.assertEqual(protpy.__status__, "Production", 
            "protpy status is not correct, got: {}".format(protpy.__status__))
        self.assertEqual(protpy.__license__, "MIT", 
            "protpy license type is not correct, got: {}".format(protpy.__license__))
        self.assertEqual(protpy.__maintainer__, "AJ McKenna", 
            "protpy maintainer is not correct, got: {}".format(protpy.__license__))
        self.assertEqual(protpy.__keywords__, ["bioinformatics", "protein engineering", 
            "python", "pypi", "machine learning", "aaindex", "protein descriptors", 
            "physiochemical descriptors", "biochemical descriptors" "structural descriptors"], 
            "protpy keywords is not correct, got: {}".format(protpy.__keywords__))

    def test_valid_descriptors(self):
        """ Testing correct list of available descriptors is included in package attribute. """
        self.assertEqual(len(protpy.all_descriptors), 14, 
            "Expected there to be 14 total descriptors, got {}".format(len(protpy.all_descriptors)))
        self.assertIsInstance(protpy.all_descriptors, list, 
            "all_descriptors should be of type list, got {}.".format(type(protpy.all_descriptors)))
        self.assertIn('aa_composition', protpy.all_descriptors, 
            "aa_composition should be in available descriptors list.")
        self.assertIn('dipeptide_composition', protpy.all_descriptors, 
            "dipeptide_composition should be in available descriptors list.")
        self.assertIn('tripeptide_composition', protpy.all_descriptors,
            "tripeptide_composition should be in available descriptors list.")
        self.assertIn('pseudo_amino_acid_composition', protpy.all_descriptors,
            "pseudo_amino_acid_composition should be in available descriptors list.")
        self.assertIn('amphiphilic_pseudo_amino_acid_composition', protpy.all_descriptors, 
            "amphiphilic_pseudo_amino_acid_composition should be in available descriptors list.")
        self.assertIn('moreaubroto_autocorrelation', protpy.all_descriptors,
            "moreaubroto_autocorrelation should be in available descriptors list.")
        self.assertIn('moran_autocorrelation', protpy.all_descriptors, 
            "moran_autocorrelation should be in available descriptors list.")
        self.assertIn('geary_autocorrelation', protpy.all_descriptors,
            "geary_autocorrelation should be in available descriptors list.")
        self.assertIn('conjoint_triad', protpy.all_descriptors,
            "conjoint_triad should be in available descriptors list.")
        self.assertIn('ctd_composition', protpy.all_descriptors,
            "ctd_composition should be in available descriptors list.")
        self.assertIn('ctd_transition', protpy.all_descriptors,
            "ctd_transition should be in available descriptors list.")
        self.assertIn('ctd_distribution', protpy.all_descriptors,
            "ctd_distribution should be in available descriptors list.")
        self.assertIn('sequence_order_coupling_number', protpy.all_descriptors,
            "sequence_order_coupling_number should be in available descriptors list.")
        self.assertIn('quasi_sequence_order', protpy.all_descriptors,
            "quasi_sequence_order should be in available descriptors list.")
