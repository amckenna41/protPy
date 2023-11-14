###########################################################################
#############      protPy - Conjoint Triad Module Tests      ##############
###########################################################################

import pandas as pd
import numpy as np
import os
import unittest
import re
from Bio import SeqIO
unittest.TestLoader.sortTestMethodsUsing = None

import protpy as protpy

class ProtpyConjointTriadTests(unittest.TestCase):
    """
    Test suite for testing conjoint triad module and functionality 
    in protpy package, including the Conjoint Triad descriptor.

    Test Cases
    ----------
    test_conjoint_triad:
        testing correct protpy conjoint triad functionality.
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

    def test_conjoint_triad(self):
        """ Testing conjoint triad descriptor attributes and functionality. """
#1.)
        conjoint_triad_seq1 = protpy.conjoint_triad(self.protein_seq1)
        #testing first 10 columns 
        conjoint_triad_seq1_expected_values = pd.DataFrame(columns=["conj_triad_111", "conj_triad_112", "conj_triad_113", "conj_triad_114", "conj_triad_115", "conj_triad_116", "conj_triad_117", "conj_triad_121", "conj_triad_122", "conj_triad_123"])
        conjoint_triad_seq1_expected_values.loc[0] = [7, 17, 11, 3, 6, 6, 2, 16, 11, 19]

        self.assertIsInstance(conjoint_triad_seq1, pd.DataFrame, 
            'Expected output to be a DataFrame, got {}.'.format(type(conjoint_triad_seq1)))
        self.assertEqual(conjoint_triad_seq1.shape, (1, 343),
            'Expected output to be of shape 1 x 343, got {}.'.format(conjoint_triad_seq1.shape)) 
        self.assertTrue(conjoint_triad_seq1.iloc[:, : 10].equals(conjoint_triad_seq1_expected_values), 
            "Output did not match expected values, got\n{}.".format(conjoint_triad_seq1.iloc[:, : 10]))
        self.assertTrue(conjoint_triad_seq1.any().isnull().sum()==0,
            'Expected output to contain no null values.')        
        self.assertTrue(all(col == np.int64 for col in list(conjoint_triad_seq1.dtypes)),
            "Expected output values to be of datatype np.int64, got\n{}.".format(list(conjoint_triad_seq1.dtypes)))
        #iterate over all columns, checking they follow naming convention using regex
        for col in list(conjoint_triad_seq1.columns):
            self.assertTrue(bool(re.match(r"conj_triad_[0-9]{3}", col)), 
                "Column name doesn't match expected regex pattern: {}.".format(col))   
#2.)
        conjoint_triad_seq2 = protpy.conjoint_triad(self.protein_seq2)
        #testing first 10 columns 
        conjoint_triad_seq2_expected_values = pd.DataFrame(columns=["conj_triad_111", "conj_triad_112", "conj_triad_113", "conj_triad_114", "conj_triad_115", "conj_triad_116", "conj_triad_117", "conj_triad_121", "conj_triad_122", "conj_triad_123"])
        conjoint_triad_seq2_expected_values.loc[0] = [1, 4, 4, 1, 2, 4, 0, 3, 5, 4]

        self.assertIsInstance(conjoint_triad_seq2, pd.DataFrame, 
            'Expected output to be a DataFrame, got {}.'.format(type(conjoint_triad_seq2)))
        self.assertEqual(conjoint_triad_seq2.shape, (1, 343),
            'Expected output to be of shape 1 x 343, got {}.'.format(conjoint_triad_seq2.shape)) 
        self.assertTrue(conjoint_triad_seq2.iloc[:, : 10].equals(conjoint_triad_seq2_expected_values), 
            "Output did not match expected values, got\n{}.".format(conjoint_triad_seq2.iloc[:, : 10]))
        self.assertTrue(conjoint_triad_seq2.any().isnull().sum()==0,
            'Expected output to contain no null values.')        
        self.assertTrue(all(col == np.int64 for col in list(conjoint_triad_seq2.dtypes)),
            "Expected output values to be of datatype np.int64, got {}.".format(list(conjoint_triad_seq2.dtypes)))
        #iterate over all columns, checking they follow naming convention using regex
        for col in list(conjoint_triad_seq2.columns):
            self.assertTrue(bool(re.match(r"conj_triad_[0-9]{3}", col)), 
                "Column name doesn't match expected regex pattern: {}.".format(col))   
#3.)
        conjoint_triad_seq3 = protpy.conjoint_triad(self.protein_seq3)
        #testing first 10 columns 
        conjoint_triad_seq3_expected_values = pd.DataFrame(columns=["conj_triad_111", "conj_triad_112", "conj_triad_113", "conj_triad_114", "conj_triad_115", "conj_triad_116", "conj_triad_117", "conj_triad_121", "conj_triad_122", "conj_triad_123"])
        conjoint_triad_seq3_expected_values.loc[0] = [0, 2, 0, 0, 0, 0, 0, 1, 3, 0]

        self.assertIsInstance(conjoint_triad_seq3, pd.DataFrame, 
            'Expected output to be a DataFrame, got {}.'.format(type(conjoint_triad_seq3)))
        self.assertEqual(conjoint_triad_seq3.shape, (1, 343),
            'Expected output to be of shape 1 x 343, got {}.'.format(conjoint_triad_seq3.shape)) 
        self.assertTrue(conjoint_triad_seq3.iloc[:, : 10].equals(conjoint_triad_seq3_expected_values), 
            "Output did not match expected values, got\n{}.".format(conjoint_triad_seq3.iloc[:, : 10]))
        self.assertTrue(conjoint_triad_seq3.any().isnull().sum()==0,
            'Expected output to contain no null values.')        
        self.assertTrue(all(col == np.int64 for col in list(conjoint_triad_seq3.dtypes)),
            "Expected output values to be of datatype np.int64, got {}.".format(list(conjoint_triad_seq3.dtypes)))
        #iterate over all columns, checking they follow naming convention using regex
        for col in list(conjoint_triad_seq3.columns):
            self.assertTrue(bool(re.match(r"conj_triad_[0-9]{3}", col)), 
                "Column name doesn't match expected regex pattern: {}.".format(col))    
#4.)
        conjoint_triad_seq4 = protpy.conjoint_triad(self.protein_seq4)
        #testing first 10 columns 
        conjoint_triad_seq4_expected_values = pd.DataFrame(columns=["conj_triad_111", "conj_triad_112", "conj_triad_113", "conj_triad_114", "conj_triad_115", "conj_triad_116", "conj_triad_117", "conj_triad_121", "conj_triad_122", "conj_triad_123"])
        conjoint_triad_seq4_expected_values.loc[0] = [0, 7, 2, 2, 1, 4, 0, 3, 8, 4]

        self.assertIsInstance(conjoint_triad_seq4, pd.DataFrame, 
            'Expected output to be a DataFrame, got {}.'.format(type(conjoint_triad_seq4)))
        self.assertEqual(conjoint_triad_seq4.shape, (1, 343),
            'Expected output to be of shape 1 x 343, got {}.'.format(conjoint_triad_seq4.shape)) 
        self.assertTrue(conjoint_triad_seq4.iloc[:, : 10].equals(conjoint_triad_seq4_expected_values), 
            "Output did not match expected values, got\n{}.".format(conjoint_triad_seq4.iloc[:, : 10]))
        self.assertTrue(conjoint_triad_seq4.any().isnull().sum()==0,
            'Expected output to contain no null values.')        
        self.assertTrue(all(col == np.int64 for col in list(conjoint_triad_seq4.dtypes)),
            "Expected output values to be of datatype np.int64, got {}.".format(list(conjoint_triad_seq4.dtypes)))
        #iterate over all columns, checking they follow naming convention using regex
        for col in list(conjoint_triad_seq4.columns):
            self.assertTrue(bool(re.match(r"conj_triad_[0-9]{3}", col)), 
                "Column name doesn't match expected regex pattern: {}.".format(col))   
#5.)
        invalid_seq5 = "ABCDEF"
        invalid_seq6 = ""
        with (self.assertRaises(ValueError)):
            conjoint_triad_seq5 = protpy.conjoint_triad(invalid_seq5)
            conjoint_triad_seq6 = protpy.conjoint_triad(invalid_seq6)
#6.)
        invalid_seq7 = 12345
        invalid_seq8 = False
        with (self.assertRaises(TypeError)):
            conjoint_triad_seq7 = protpy.conjoint_triad(invalid_seq7)
            conjoint_triad_seq8 = protpy.conjoint_triad(invalid_seq8)