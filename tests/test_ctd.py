###########################################################################
##################      protPy - CTD Module Tests      ####################
###########################################################################

import pandas as pd
import numpy as np
import os
import unittest
import re
from Bio import SeqIO
unittest.TestLoader.sortTestMethodsUsing = None
import protpy as protpy

class ProtPyCTDTests(unittest.TestCase):
    """
    Test suite for testing CTD (Composition, Transition, Distribution) 
    module and functionality in protpy package. 

    Test Cases
    ----------
    test_ctd:
        testing correct protpy CTD functionality.
    test_ctd_composition:
        testing correct protpy composition functionality.
    test_ctd_transition:
        testing correct protpy transition functionality.
    test_ctd_distribution:
        testing correct protpy distribution functionality.
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

    def test_ctd(self):
        """ Testing CTD descriptor attributes and functionality. """   
        properties = ["hydrophobicity", "normalized_vdwv", "polarity", "charge",
            "secondary_struct", "solvent_accessibility", "polarizability"]
 #1.)
        ctd_seq1 = protpy.ctd_(self.protein_seq1)

        self.assertEqual(ctd_seq1.shape, (1, 147), 'Descriptor not of correct shape.') 
        self.assertIsInstance(ctd_seq1, pd.DataFrame, 'Descriptor not of type DataFrame.')
        self.assertTrue(ctd_seq1.any().isnull().sum()==0, 'Descriptor should not contain any null values.')  
        self.assertTrue(all(col == np.float64 for col in list(ctd_seq1.dtypes)), "")

        #iterate over all columns, checking they follow naming convention using regex
        for col in list(ctd_seq1.columns):
            matching_col = False
            for prop in properties:
                if (col.startswith(prop)):
                    matching_col = True
                    self.assertTrue(((bool(re.search(prop + r"_CTD_[A-Z]{1}_[0-9][0-9]", col))) or 
                        (bool(re.search(prop + r"_CTD_[A-Z]{1}_[0-9][0-9]_[0-9][0-9]", col)))), 
                            "Column name does not follow expected format: {}.".format(col))
            self.assertTrue(matching_col, "Column name's property name not found and doesn't match format: {}.".format(matching_col))

    def test_ctd_composition(self):
        """ Testing CTD Composition descriptor attributes and functionality. """   
        properties = ["hydrophobicity", "normalized_vdwv", "polarity", "charge",
            "secondary_struct", "solvent_accessibility", "polarizability"]
#1.)
        for prop in properties:
            ctd_composition_seq1 = protpy.ctd_composition(self.protein_seq1, property=prop)

            self.assertEqual(ctd_composition_seq1.shape, (1, 3), 'Descriptor not of correct shape.') 
            self.assertIsInstance(ctd_composition_seq1, pd.DataFrame, 'Descriptor not of type DataFrame.')
            self.assertTrue(ctd_composition_seq1.any().isnull().sum()==0, 'Descriptor should not contain any null values.')  
            self.assertTrue(all(col == np.float64 for col in list(ctd_composition_seq1.dtypes)), "")
            
            #iterate over all columns, checking they follow naming convention using regex
            for col in list(ctd_composition_seq1.columns):
                matching_col = False
                if (col.startswith(prop)):
                    matching_col = True
                    self.assertTrue((bool(re.search(prop + r"_CTD_[A-Z]{1}_[0-9][0-9]", col))), 
                        "Column name does not follow expected format: {}.".format(col))
                self.assertTrue(matching_col, "Column name's property name not found and doesn't match format: {}.".format(matching_col))

    def test_ctd_distribution(self):
        """ Testing CTD Distribution descriptor attributes and functionality. """   
        properties = ["hydrophobicity", "normalized_vdwv", "polarity", "charge",
            "secondary_struct", "solvent_accessibility", "polarizability"]
#1.)
        for prop in properties:
            ctd_distribution_seq1 = protpy.ctd_distribution(self.protein_seq1, property=prop)

            self.assertEqual(ctd_distribution_seq1.shape, (1, 15), 'Descriptor not of correct shape.') 
            self.assertIsInstance(ctd_distribution_seq1, pd.DataFrame, 'Descriptor not of type DataFrame.')
            self.assertTrue(ctd_distribution_seq1.any().isnull().sum()==0, 'Descriptor should not contain any null values.')  
            self.assertTrue(all(col == np.float64 for col in list(ctd_distribution_seq1.dtypes)), "")
            
            #iterate over all columns, checking they follow naming convention using regex
            for col in list(ctd_distribution_seq1.columns):
                matching_col = False
                if (col.startswith(prop)):
                    matching_col = True
                    self.assertTrue((bool(re.search(prop + r"_CTD_[A-Z]{1}_[0-9][0-9]", col))), 
                        "Column name does not follow expected format: {}.".format(col))
                self.assertTrue(matching_col, "Column name's property name not found and doesn't match format: {}.".format(matching_col))

    def test_ctd_transition(self):
        """ Testing CTD Transition descriptor attributes and functionality. """   
        properties = ["hydrophobicity", "normalized_vdwv", "polarity", "charge",
            "secondary_struct", "solvent_accessibility", "polarizability"]
#1.)    
        for prop in properties:
            ctd_transition_seq1 = protpy.ctd_transition(self.protein_seq1, property=prop)

            self.assertEqual(ctd_transition_seq1.shape, (1, 3), 'Descriptor not of correct shape.') 
            self.assertIsInstance(ctd_transition_seq1, pd.DataFrame, 'Descriptor not of type DataFrame.')
            self.assertTrue(ctd_transition_seq1.any().isnull().sum()==0, 'Descriptor should not contain any null values.')  
            self.assertTrue(all(col == np.float64 for col in list(ctd_transition_seq1.dtypes)), "")
            
            #iterate over all columns and check its name follows expected format
            for col in list(ctd_transition_seq1.columns):
                matching_col = False
                if (col.startswith(prop)):
                    matching_col = True
                    self.assertTrue(((bool(re.search(prop + r"_CTD_[A-Z]{1}_[0-9][0-9]", col))) or 
                        (bool(re.search(prop + r"_CTD_[A-Z]{1}_[0-9][0-9]_[0-9][0-9]", col)))), 
                            "Column name does not follow expected format: {}.".format(col))
                self.assertTrue(matching_col, "Column name's property name not found and doesn't match format: {}.".format(matching_col))