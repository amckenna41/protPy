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

class ProtpyCTDTests(unittest.TestCase):
    """
    Test suite for testing CTD (Composition, Transition, Distribution) 
    module and functionality in protpy package, including the CTD 
    descriptors.

    Test Cases
    ----------
    test_ctd:
        testing correct protpy overall CTD functionality.
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
        with open(os.path.join("tests", "test_data", "test_fasta1.fasta")) as pro:
            self.protein_seq1 = str(next(SeqIO.parse(pro,'fasta')).seq)

        with open(os.path.join("tests", "test_data", "test_fasta2.fasta")) as pro:
            self.protein_seq2 = str(next(SeqIO.parse(pro,'fasta')).seq)
        
        with open(os.path.join("tests", "test_data", "test_fasta3.fasta")) as pro:
            self.protein_seq3 = str(next(SeqIO.parse(pro,'fasta')).seq)

        with open(os.path.join("tests", "test_data", "test_fasta4.fasta")) as pro:
            self.protein_seq4 = str(next(SeqIO.parse(pro,'fasta')).seq)

        self.all_protein_seqs = [self.protein_seq1, self.protein_seq2, self.protein_seq3, self.protein_seq4]

    def test_ctd(self):
        """ Testing CTD descriptor attributes and functionality. """   
        properties = ["hydrophobicity", "normalized_vdwv", "polarity", "charge",
            "secondary_struct", "solvent_accessibility", "polarizability"]
 #1.)   
        for seq in self.all_protein_seqs:
            ctd_seq1 = protpy.ctd_(seq, all_ctd=True) #using all properties

            self.assertIsInstance(ctd_seq1, pd.DataFrame,
                'Expected output to be a DataFrame, got {}.'.format(type(ctd_seq1)))
            self.assertEqual(ctd_seq1.shape, (1, 147), 
                'Expected output to be of shape 1 x 147, got {}.'.format(ctd_seq1.shape)) 
            self.assertTrue(ctd_seq1.any().isnull().sum()==0,
                'Expected output to contain no null values.')        
            self.assertTrue(all(col == np.float64 for col in list(ctd_seq1.dtypes)),
                "Expected output values to be of datatype np.float64, got {}.".format(list(ctd_seq1.dtypes)))

            #iterate over all columns, checking they follow naming convention using regex
            for col in list(ctd_seq1.columns):
                matching_col = False
                for prop in properties:
                    if (col.endswith(prop)):
                        matching_col = True
                        self.assertTrue(((bool(re.search(r"CTD_[A-Z]_[0-9]{2}_" + prop, col))) or 
                            (bool(re.search(r"CTD_[A-Z]_[0-9]{2}_[0-9]{2}_" + prop, col))) or 
                            (bool(re.search(r"CTD_[A-Z]_[0-9]{2}_[0-9]{3}_" + prop, col)))), 
                                "Column name does not follow expected regex format: {}.".format(col))
                self.assertTrue(matching_col, 
                    "Column name {} not found in list of available properties:\n{}.".format(col, properties))
#2.)    
            for prop in properties:
                ctd_seq2 = protpy.ctd_(seq, property=prop, all_ctd=False) #using one property at a time

                self.assertIsInstance(ctd_seq2, pd.DataFrame,
                    'Expected output to be a DataFrame, got {}.'.format(type(ctd_seq2)))
                self.assertEqual(ctd_seq2.shape, (1, 21), 
                    'Expected output to be of shape 1 x 21, got {}.'.format(ctd_seq2.shape))
                self.assertTrue(ctd_seq2.any().isnull().sum()==0,
                    'Expected output to contain no null values.')        
                self.assertTrue(all(col == np.float64 for col in list(ctd_seq2.dtypes)),
                    "Expected output values to be of datatype np.float64, got {}.".format(list(ctd_seq2.dtypes)))

                #iterate over all columns, checking they follow naming convention using regex
                for col in list(ctd_seq2.columns):
                    matching_col = False
                    for prop in properties:
                        if (col.endswith(prop)):
                            matching_col = True
                            self.assertTrue(((bool(re.search(r"CTD_[A-Z]_[0-9]{2}_" + prop, col))) or 
                                (bool(re.search(r"CTD_[A-Z]_[0-9]{2}_[0-9]{2}_" + prop, col))) or 
                                (bool(re.search(r"CTD_[A-Z]_[0-9]{2}_[0-9]{3}_" + prop, col)))), 
                                    "Column name does not follow expected regex format: {}.".format(col))
                    self.assertTrue(matching_col, 
                        "Column name {} not found in list of available properties:\n{}.".format(col, properties))

    def test_ctd_composition(self):
        """ Testing CTD Composition descriptor attributes and functionality. """   
        properties = ["hydrophobicity", "normalized_vdwv", "polarity", "charge",
            "secondary_struct", "solvent_accessibility", "polarizability"]
#1.)
        for seq in self.all_protein_seqs:
            for prop in properties:
                ctd_composition_seq1 = protpy.ctd_composition(seq, property=prop)

                self.assertIsInstance(ctd_composition_seq1, pd.DataFrame,
                    'Expected output to be a DataFrame, got {}.'.format(type(ctd_composition_seq1)))
                self.assertEqual(ctd_composition_seq1.shape, (1, 3), 
                    'Expected output to be of shape 1 x 3, got {}.'.format(ctd_composition_seq1.shape)) 
                self.assertTrue(ctd_composition_seq1.any().isnull().sum()==0,
                    'Expected output to contain no null values.')        
                self.assertTrue(all(col == np.float64 for col in list(ctd_composition_seq1.dtypes)),
                    "Expected output values to be of datatype np.float64, got {}.".format(list(ctd_composition_seq1.dtypes)))

                #iterate over all columns, checking they follow naming convention using regex
                for col in list(ctd_composition_seq1.columns):
                    matching_col = False
                    for prop in properties:
                        if (col.endswith(prop)):
                            matching_col = True
                            self.assertTrue(((bool(re.search(r"CTD_C_[0-9]{2}_" + prop, col))) or 
                                (bool(re.search(r"CTD_C_[0-9]{2}_[0-9]{2}_" + prop, col))) or 
                                (bool(re.search(r"CTD_C_[0-9]{2}_[0-9]{3}_" + prop, col)))), 
                                    "Column name does not follow expected regex format: {}.".format(col))
            self.assertTrue(matching_col, 
                "Column name {} not found in list of available properties:\n{}.".format(col, properties))

    def test_ctd_transition(self):
        """ Testing CTD Transition descriptor attributes and functionality. """   
        properties = ["hydrophobicity", "normalized_vdwv", "polarity", "charge",
            "secondary_struct", "solvent_accessibility", "polarizability"]
#1.)    
        for seq in self.all_protein_seqs:
            for prop in properties:
                ctd_transition_seq1 = protpy.ctd_transition(seq, property=prop)

                self.assertIsInstance(ctd_transition_seq1, pd.DataFrame, 
                    'Expected output to be a DataFrame, got {}.'.format(type(ctd_transition_seq1)))
                self.assertEqual(ctd_transition_seq1.shape, (1, 3), 
                    'Expected output to be of shape 1 x 3, got {}.'.format(ctd_transition_seq1.shape))
                self.assertTrue(ctd_transition_seq1.any().isnull().sum()==0,
                    'Expected output to contain no null values.')        
                self.assertTrue(all(col == np.float64 for col in list(ctd_transition_seq1.dtypes)), 
                    "Expected output values to be of datatype np.float64, got {}.".format(list(ctd_transition_seq1.dtypes)))
#2.)
                #iterate over all columns and check its name follows expected format
                for col in list(ctd_transition_seq1.columns):
                    matching_col = False
                    if (col.endswith(prop)):
                        matching_col = True
                        self.assertTrue(((bool(re.search(r"CTD_T_[0-9]{2}_" + prop, col))) or 
                            (bool(re.search(r"CTD_T_[0-9]{2}_[0-9]{2}_" + prop, col)))), 
                                "Column name does not follow expected regex format: {}.".format(col))
            self.assertTrue(matching_col, 
                "Column name {} not found in list of available properties:\n{}.".format(col, properties))

    def test_ctd_distribution(self):
        """ Testing CTD Distribution descriptor attributes and functionality. """   
        properties = ["hydrophobicity", "normalized_vdwv", "polarity", "charge",
            "secondary_struct", "solvent_accessibility", "polarizability"]
#1.)
        for seq in self.all_protein_seqs:
            for prop in properties:
                ctd_distribution_seq1 = protpy.ctd_distribution(seq, property=prop)

                self.assertIsInstance(ctd_distribution_seq1, pd.DataFrame, 
                    'Expected output to be a DataFrame, got {}.'.format(type(ctd_distribution_seq1)))
                self.assertEqual(ctd_distribution_seq1.shape, (1, 15), 
                    'Expected output to be of shape {}, got {}.'.format((1, 15), ctd_distribution_seq1.shape)) 
                self.assertTrue(ctd_distribution_seq1.any().isnull().sum()==0,
                    'Expected output to contain no null values.')        
                self.assertTrue(all(col == np.float64 for col in list(ctd_distribution_seq1.dtypes)), 
                    "Expected output values to be of datatype np.float64, got {}.".format(list(ctd_distribution_seq1.dtypes)))
                
                #iterate over all columns and check its name follows expected format
                for col in list(ctd_distribution_seq1.columns):
                    matching_col = False
                    if (col.endswith(prop)):
                        matching_col = True
                        self.assertTrue(((bool(re.search(r"CTD_D_[0-9][0-9]_" + prop, col))) or (bool(re.search(r"CTD_D_[0-9]{2}_[0-9]{3}_" + prop, col)))), 
                            "Column name does not follow expected format: {}.".format(col))
            self.assertTrue(matching_col, 
                "Column name {} not found in list of available properties:\n{}.".format(col, properties))