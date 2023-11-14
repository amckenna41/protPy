###########################################################################
###############      protPy - Composition Module Tests      ###############
###########################################################################

import pandas as pd
import numpy as np
import os
import unittest
import re
from Bio import SeqIO
unittest.TestLoader.sortTestMethodsUsing = None

import protpy as protpy

class ProtpyCompositionTests(unittest.TestCase):
    """
    Test suite for testing composition module and functionality in protpy package, 
    including the Amino Acid, Dipeptide and Tripeptide Composition descriptors as 
    well as the Pseudo and Amphipillic Pseudo Amino Acid descriptors.

    Test Cases
    ----------
    test_amino_acid_composition:
        testing correct protpy amino acid composition functionality.
    test_dipeptide_composition:
        testing correct protpy dipeptide amino acid composition functionality.
    test_tripeptide_composition:
        testing correct protpy tripeptide amino acid composition functionality.
    test_pseudo_amino_acid_composition:
        testing correct protpy pseudo amino acid composition functionality.
    test_amphiphilic_pseudo_amino_acid_composition:
        testing correct protpy amphophilic amino acid composition functionality.
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

        #list of canonical amino acids
        self.amino_acids = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", 
                            "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]

    def test_amino_acid_composition(self): 
        """ Testing Amino Acid Composition protein descriptor attributes and methods. """
#1.)
        amino_acid_composition_seq1 = protpy.amino_acid_composition(self.protein_seq1)
        amino_acid_composition_seq1_expected_values = pd.DataFrame(columns=["A", "C", "D", "E", "F", "G", "H", 
            "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"])
        amino_acid_composition_seq1_expected_values.loc[0] = [6.693, 3.108, 5.817, 3.347, 6.614, 6.295, 1.195, 
            6.215, 4.781, 7.888, 1.594, 6.454, 4.542, 4.382, 3.108, 7.649, 7.888, 7.251, 0.876, 4.303]

        self.assertIsInstance(amino_acid_composition_seq1, pd.DataFrame,
            'Expected output to be of type DataFrame, got {}.'.format(type(amino_acid_composition_seq1)))
        self.assertEqual(amino_acid_composition_seq1.shape, (1, 20),
            'Expected output to be of shape {}, got {}.'.format((1, 20), amino_acid_composition_seq1.shape)) 
        self.assertEqual(self.amino_acids, list(amino_acid_composition_seq1.columns), 
            'Incorrect column values found in output:\n{}.'.format(list(amino_acid_composition_seq1.columns)))
        self.assertTrue(amino_acid_composition_seq1.any().isnull().sum()==0,
            'Expected output to contain no null values.')        
        self.assertTrue(amino_acid_composition_seq1.equals(amino_acid_composition_seq1_expected_values), 
            "Output did not match expected values, got\n{}.".format(amino_acid_composition_seq1))
        self.assertTrue(all(col == np.float64 for col in list(amino_acid_composition_seq1.dtypes)),
                "Expected output values to be of datatype np.float64, got {}.".format(list(amino_acid_composition_seq1.dtypes)))
#2.)
        amino_acid_composition_seq2 = protpy.amino_acid_composition(self.protein_seq2)
        amino_acid_composition_seq2_expected_values = pd.DataFrame(columns=["A", "C", "D", "E", "F", "G", "H", 
            "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"])
        amino_acid_composition_seq2_expected_values.loc[0] = [8.057, 0.0, 5.213, 3.318, 3.081, 10.664, 1.185, 
            2.607, 6.872, 6.161, 1.659, 5.924, 7.346, 8.057, 7.346, 8.294, 7.82, 2.607, 1.185, 2.607]

        self.assertIsInstance(amino_acid_composition_seq2, pd.DataFrame,
            'Expected output to be of type DataFrame, got {}.'.format(type(amino_acid_composition_seq2)))
        self.assertEqual(amino_acid_composition_seq2.shape, (1, 20),
            'Expected output to be of shape {}, got {}.'.format((1, 20), amino_acid_composition_seq2.shape)) 
        self.assertEqual(self.amino_acids, list(amino_acid_composition_seq2.columns), 
            'Incorrect column values found in output:\n{}'.format(list(amino_acid_composition_seq2.columns)))
        self.assertTrue(amino_acid_composition_seq2.any().isnull().sum()==0,
            'Expected output to contain no null values.')        
        self.assertTrue(amino_acid_composition_seq2.equals(amino_acid_composition_seq2_expected_values), 
            "Output did not match expected values, got\n{}.".format(amino_acid_composition_seq2))
        self.assertTrue(all(col == np.float64 for col in list(amino_acid_composition_seq2.dtypes)),
                "Expected output values to be of datatype np.float64, got {}.".format(list(amino_acid_composition_seq2.dtypes)))
#3.)
        amino_acid_composition_seq3 = protpy.amino_acid_composition(self.protein_seq3)
        amino_acid_composition_seq3_expected_values = pd.DataFrame(columns=["A", "C", "D", "E", "F", "G", "H", 
            "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"])
        amino_acid_composition_seq3_expected_values.loc[0] = [5.263, 3.947, 1.316, 3.947, 5.263, 2.632, 0.0, 
            3.947, 2.632, 18.421, 1.316, 6.579, 2.632, 0.0, 2.632, 9.211, 6.579, 18.421, 0.0, 5.263]

        self.assertIsInstance(amino_acid_composition_seq3, pd.DataFrame,
            'Expected output to be of type DataFrame, got {}.'.format(type(amino_acid_composition_seq3)))
        self.assertEqual(amino_acid_composition_seq3.shape, (1, 20),
            'Expected output to be of shape {}, got {}.'.format((1, 20), amino_acid_composition_seq3.shape)) 
        self.assertEqual(self.amino_acids, list(amino_acid_composition_seq3.columns), 
            'Incorrect column values found in output:\n{}'.format(list(amino_acid_composition_seq3.columns)))
        self.assertTrue(amino_acid_composition_seq3.any().isnull().sum()==0,
            'Expected output to contain no null values.')        
        self.assertTrue(amino_acid_composition_seq3.equals(amino_acid_composition_seq3_expected_values), 
            "Output did not match expected values, got\n{}.".format(amino_acid_composition_seq3))
        self.assertTrue(all(col == np.float64 for col in list(amino_acid_composition_seq3.dtypes)),
                "Expected output values to be of datatype np.float64, got {}.".format(list(amino_acid_composition_seq3.dtypes)))
#4.)
        amino_acid_composition_seq4 = protpy.amino_acid_composition(self.protein_seq4)
        amino_acid_composition_seq4_expected_values = pd.DataFrame(columns=["A", "C", "D", "E", "F", "G", "H", 
            "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"])
        amino_acid_composition_seq4_expected_values.loc[0] = [8.831, 0.0, 5.728, 2.864, 3.103, 10.263, 0.955, 
            3.341, 7.399, 6.444, 1.671, 5.251, 6.683, 8.353, 6.921, 8.831, 7.637, 1.909, 1.193, 2.625]

        self.assertIsInstance(amino_acid_composition_seq4, pd.DataFrame,
            'Expected output to be of type DataFrame, got {}.'.format(type(amino_acid_composition_seq4)))
        self.assertEqual(amino_acid_composition_seq4.shape, (1, 20),
            'Expected output to be of shape {}, got {}.'.format((1, 20), amino_acid_composition_seq4.shape)) 
        self.assertEqual(self.amino_acids, list(amino_acid_composition_seq4.columns), 
            'Incorrect column values found in output:\n{}'.format(list(amino_acid_composition_seq4.columns)))
        self.assertTrue(amino_acid_composition_seq4.any().isnull().sum()==0,
            'Expected output to contain no null values.')        
        self.assertTrue(amino_acid_composition_seq4.equals(amino_acid_composition_seq4_expected_values), 
            "Output did not match expected values, got\n{}.".format(amino_acid_composition_seq4))
        self.assertTrue(all(col == np.float64 for col in list(amino_acid_composition_seq4.dtypes)),
                "Expected output values to be of datatype np.float64, got {}.".format(list(amino_acid_composition_seq4.dtypes)))
#5.)
        invalid_seq5 = "ABCDEF"
        invalid_seq6 = "OOOOO"
        with (self.assertRaises(ValueError)):
            amino_acid_composition_seq5 = protpy.amino_acid_composition(invalid_seq5)
            amino_acid_composition_seq6 = protpy.amino_acid_composition(invalid_seq6)
#6.)
        invalid_seq7 = 12345
        invalid_seq8 = True
        with (self.assertRaises(TypeError)):
            amino_acid_composition_seq7 = protpy.amino_acid_composition(invalid_seq7)
            amino_acid_composition_seq8 = protpy.amino_acid_composition(invalid_seq8)

    def test_dipeptide_composition(self): 
        """ Testing Dipeptide Composition protein descriptor attributes and functionality. """
#1.)
        dipeptide_comp_seq1 = protpy.dipeptide_composition(self.protein_seq1)
        #testing first 10 columns 
        dipeptide_comp_seq1_expected_values = pd.DataFrame(columns=["AA", "AC", "AD", "AE", "AF", "AG", "AH", "AI", "AK", "AL"])
        dipeptide_comp_seq1_expected_values.loc[0] = [0.72, 0.16, 0.48, 0.4, 0.24, 0.48, 0.0, 0.56, 0.16, 0.48] 

        self.assertIsInstance(dipeptide_comp_seq1, pd.DataFrame,
            'Expected output to be of type DataFrame, got {}.'.format(type(dipeptide_comp_seq1)))
        self.assertEqual(dipeptide_comp_seq1.shape, (1, 400),
            'Expected output to be of shape {}, got {}.'.format((1, 400), dipeptide_comp_seq1.shape)) 
        self.assertTrue(dipeptide_comp_seq1.iloc[:, : 10].equals(dipeptide_comp_seq1_expected_values), 
            "Output did not match expected values, got\n{}.".format(dipeptide_comp_seq1.iloc[:, : 10]))
        for col in list(dipeptide_comp_seq1.columns):
            #check all columns follow pattern of XY where x & y are amino acids 
            self.assertTrue(bool(re.match(r'^[A-Z]{2}$', col)), 
                "Column doesn't follow correct regex naming convention: {}.".format(col))
            self.assertIn(col[0], self.amino_acids, 
                "Amino acid in column {} not found in list of amino acids:\n{}.".format(col[0], self.amino_acids))
            self.assertIn(col[1], self.amino_acids, 
                "Amino acid in column {} not found in list of amino acids:\n{}.".format(col[1], self.amino_acids))
        self.assertTrue(dipeptide_comp_seq1.any().isnull().sum()==0, 
            'Expected output to contain no null values.')        
        self.assertTrue(all(col == np.float64 for col in list(dipeptide_comp_seq1.dtypes)), 
            "Expected output values to be of datatype np.float64, got {}.".format(list(dipeptide_comp_seq1.dtypes)))
#2.)
        dipeptide_comp_seq2 = protpy.dipeptide_composition(self.protein_seq2)
        #testing first 10 columns 
        dipeptide_comp_seq2_expected_values = pd.DataFrame(columns=["AA", "AC", "AD", "AE", "AF", "AG", "AH", "AI", "AK", "AL"])
        dipeptide_comp_seq2_expected_values.loc[0] = [0.71, 0.0, 0.48, 0.48, 0.48, 0.0, 0.0, 0.24, 0.0, 0.95]
        
        self.assertIsInstance(dipeptide_comp_seq2, pd.DataFrame,
            'Expected output to be of type DataFrame, got {}.'.format(type(dipeptide_comp_seq2)))
        self.assertEqual(dipeptide_comp_seq2.shape, (1, 400),
            'Expected output to be of shape {}, got {}.'.format((1, 400), dipeptide_comp_seq2.shape)) 
        self.assertTrue(dipeptide_comp_seq2.iloc[:, : 10].equals(dipeptide_comp_seq2_expected_values), 
            "Output did not match expected values, got\n{}.".format(dipeptide_comp_seq2.iloc[:, : 10]))
        for col in list(dipeptide_comp_seq2.columns):
            #check all columns follow pattern of XY where x & y are amino acids 
            self.assertTrue(bool(re.match(r'^[A-Z]{2}$', col)), 
                "Column doesn't follow correct regex naming convention: {}.".format(col))
            self.assertIn(col[0], self.amino_acids, 
                "Amino acid in column {} not found in list of amino acids:\n{}.".format(col[0], self.amino_acids))
            self.assertIn(col[1], self.amino_acids, 
                "Amino acid column {} not found in list of amino acids:\n{}.".format(col[1], self.amino_acids))
        self.assertTrue(dipeptide_comp_seq2.any().isnull().sum()==0, 
            'Expected output to contain no null values.')        
        self.assertTrue(all(col == np.float64 for col in list(dipeptide_comp_seq2.dtypes)), 
            "Expected output values to be of datatype np.float64, got {}.".format(list(dipeptide_comp_seq2.dtypes)))
#3.)
        dipeptide_comp_seq3 = protpy.dipeptide_composition(self.protein_seq3)
        #testing first 10 columns 
        dipeptide_comp_seq3_expected_values = pd.DataFrame(columns=["AA", "AC", "AD", "AE", "AF", "AG", "AH", "AI", "AK", "AL"])
        dipeptide_comp_seq3_expected_values.loc[0] = [0.0, 0.0, 0.0, 0.0, 1.33, 0.0, 0.0, 1.33, 0.0, 1.33]

        self.assertIsInstance(dipeptide_comp_seq3, pd.DataFrame,
            'Expected output to be of type DataFrame, got {}.'.format(type(dipeptide_comp_seq3)))
        self.assertEqual(dipeptide_comp_seq3.shape, (1, 400),
            'Expected output to be of shape {}, got {}.'.format((1, 400), dipeptide_comp_seq3.shape)) 
        self.assertTrue(dipeptide_comp_seq3.iloc[:, : 10].equals(dipeptide_comp_seq3_expected_values), 
            "Output did not match expected values, got\n{}.".format(dipeptide_comp_seq3.iloc[:, : 10]))
        for col in list(dipeptide_comp_seq3.columns):
            #check all columns follow pattern of XY where x & y are amino acids 
            self.assertTrue(bool(re.match(r'^[A-Z]{2}$', col)), 
                "Column doesn't follow correct regex naming convention: {}.".format(col))
            self.assertIn(col[0], self.amino_acids, 
                "Amino acid column {} not found in list of amino acids:\n{}.".format(col[0], self.amino_acids))
            self.assertIn(col[1], self.amino_acids, 
                "Amino acid column {} not found in list of amino acids:\n{}.".format(col[1], self.amino_acids))
        self.assertTrue(dipeptide_comp_seq3.any().isnull().sum()==0, 
            'Expected output to contain no null values.')        
        self.assertTrue(all(col == np.float64 for col in list(dipeptide_comp_seq3.dtypes)), 
            "Expected output values to be of datatype np.float64, got {}.".format(list(dipeptide_comp_seq3.dtypes)))
#4.)
        dipeptide_comp_seq4 = protpy.dipeptide_composition(self.protein_seq4)
        #testing first 10 columns 
        dipeptide_comp_seq4_expected_values = pd.DataFrame(columns=["AA", "AC", "AD", "AE", "AF", "AG", "AH", "AI", "AK", "AL"])
        dipeptide_comp_seq4_expected_values.loc[0] = [0.96, 0.0, 0.72, 0.48, 0.48, 0.48, 0.0, 0.48, 0.0, 1.2]

        self.assertIsInstance(dipeptide_comp_seq4, pd.DataFrame,
            'Expected output to be of type DataFrame, got {}.'.format(type(dipeptide_comp_seq4)))
        self.assertEqual(dipeptide_comp_seq3.shape, (1, 400),
            'Expected output to be of shape {}, got {}.'.format((1, 400), dipeptide_comp_seq4.shape)) 
        self.assertTrue(dipeptide_comp_seq4.iloc[:, : 10].equals(dipeptide_comp_seq4_expected_values), 
            "Output did not match expected values, got\n{}.".format(dipeptide_comp_seq4.iloc[:, : 10]))
        for col in list(dipeptide_comp_seq4.columns):
            #check all columns follow pattern of XY where x & y are amino acids 
            self.assertTrue(bool(re.match(r'^[A-Z]{2}$', col)), 
                "Column doesn't follow correct regex naming convention: {}.".format(col))
            self.assertIn(col[0], self.amino_acids, 
                "Amino acid column {} not found in list of amino acids:\n{}.".format(col[0], self.amino_acids))
            self.assertIn(col[1], self.amino_acids, 
                "Amino acid column {} not found in list of amino acids:\n{}.".format(col[1], self.amino_acids))
        self.assertTrue(dipeptide_comp_seq4.any().isnull().sum()==0, 
            'Expected output to contain no null values.')        
        self.assertTrue(all(col == np.float64 for col in list(dipeptide_comp_seq4.dtypes)), 
            "Expected output values to be of datatype np.float64, **got {}.".format(list(dipeptide_comp_seq4.dtypes)))
#5.)
        invalid_seq5 = "ABCDEF"
        invalid_seq6 = "OOOOO"
        with (self.assertRaises(ValueError)):
            dipeptide_composition_seq5 = protpy.dipeptide_composition(invalid_seq5)
            dipeptide_composition_seq6 = protpy.dipeptide_composition(invalid_seq6)
#6.)
        invalid_seq7 = 12345
        invalid_seq8 = True
        with (self.assertRaises(TypeError)):
            dipeptide_composition_seq7 = protpy.dipeptide_composition(invalid_seq7)
            dipeptide_composition_seq8 = protpy.dipeptide_composition(invalid_seq8)

    # @unittest.skip("Descriptor can take quite a bit of time to calculate therefore skipping.")
    def test_tripeptide_composition(self):
        """ Testing Tripeptide Composition protein descriptor attributes and functionality. """
#1.)
        tripeptide_comp_seq1 = protpy.tripeptide_composition(self.protein_seq1)
        #testing first 10 columns 
        tripeptide_comp_seq1_expected_values = pd.DataFrame(columns=["AAA", "AAC", "AAD", "AAE", "AAF", "AAG", "AAH", "AAI", "AAK", "AAL"])
        tripeptide_comp_seq1_expected_values.loc[0] = [1, 0, 0, 2, 0, 0, 0, 0, 0, 2]

        self.assertIsInstance(tripeptide_comp_seq1, pd.DataFrame,
            'Expected output to be of type DataFrame, got {}.'.format(type(tripeptide_comp_seq1)))
        self.assertEqual(tripeptide_comp_seq1.shape, (1, 8000),
            'Expected output to be of shape {}, got {}.'.format((1, 8000), tripeptide_comp_seq1.shape)) 
        self.assertTrue(tripeptide_comp_seq1.iloc[:, : 10].equals(tripeptide_comp_seq1_expected_values), 
            "Output did not match expected values, got\n{}.".format(tripeptide_comp_seq1.iloc[:, : 10]))
        for col in list(tripeptide_comp_seq1.columns):
            #check all columns follow pattern of XY where x & y are amino acids 
            self.assertTrue(bool(re.match(r'^[A-Z]{3}$', col)), 
                "Column doesn't follow correct regex naming convention: {}.".format(col))      
            self.assertIn(col[0], self.amino_acids,
                "Amino acid {} not found in list of amino acids:\n{}.".format(col[0], self.amino_acids))
            self.assertIn(col[1], self.amino_acids,
                "Amino acid {} not found in list of amino acids:\n{}.".format(col[1], self.amino_acids))
            self.assertIn(col[2], self.amino_acids,
                "Amino acid {} not found in list of amino acids:\n{}.".format(col[2], self.amino_acids))
        self.assertTrue(tripeptide_comp_seq1.any().isnull().sum()==0,
            'Expected output to contain no null values.')        
        self.assertTrue(all(col == np.int64 for col in list(tripeptide_comp_seq1.dtypes)),
            "Expected output values to be of datatype np.float64, got {}.".format(list(tripeptide_comp_seq1.dtypes)))
#2.)
        tripeptide_comp_seq2 = protpy.tripeptide_composition(self.protein_seq2)
        #testing first 10 columns 
        tripeptide_comp_seq2_expected_values = pd.DataFrame(columns=["AAA", "AAC", "AAD", "AAE", "AAF", "AAG", "AAH", "AAI", "AAK", "AAL"])
        tripeptide_comp_seq2_expected_values.loc[0] = [0, 0, 1, 1, 0, 0, 0, 0, 0, 0]

        self.assertIsInstance(tripeptide_comp_seq2, pd.DataFrame,
            'Expected output to be of type DataFrame, got {}.'.format(type(tripeptide_comp_seq2)))
        self.assertEqual(tripeptide_comp_seq2.shape, (1, 8000),
            'Expected output to be of shape {}, got {}.'.format((1, 8000), tripeptide_comp_seq2.shape)) 
        self.assertTrue(tripeptide_comp_seq2.iloc[:, : 10].equals(tripeptide_comp_seq2_expected_values), 
            "Output did not match expected values, got\n{}.".format(tripeptide_comp_seq2.iloc[:, : 10]))
        for col in list(tripeptide_comp_seq2.columns):
            #check all columns follow pattern of XY where x & y are amino acids 
            self.assertTrue(bool(re.match(r'^[A-Z]{3}$', col)), 
                "Column doesn't follow correct regex naming convention: {}.".format(col))      
            self.assertIn(col[0], self.amino_acids,
                "Amino acid {} not found in list of amino acids:\n{}.".format(col[0], self.amino_acids))
            self.assertIn(col[1], self.amino_acids,
                "Amino acid {} not found in list of amino acids:\n{}.".format(col[1], self.amino_acids))
            self.assertIn(col[2], self.amino_acids,
                "Amino acid {} not found in list of amino acids:\n{}.".format(col[2], self.amino_acids))
        self.assertTrue(tripeptide_comp_seq2.any().isnull().sum()==0,
            'Expected output to contain no null values.')        
        self.assertTrue(all(col == np.int64 for col in list(tripeptide_comp_seq2.dtypes)),
            "Expected output values to be of datatype np.float64, got {}.".format(list(tripeptide_comp_seq2.dtypes)))
#3.)
        tripeptide_comp_seq3 = protpy.tripeptide_composition(self.protein_seq3)
        #testing first 10 columns 
        tripeptide_comp_seq3_expected_values = pd.DataFrame(columns=["AAA", "AAC", "AAD", "AAE", "AAF", "AAG", "AAH", "AAI", "AAK", "AAL"])
        tripeptide_comp_seq3_expected_values.loc[0] = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

        self.assertIsInstance(tripeptide_comp_seq3, pd.DataFrame,
            'Expected output to be of type DataFrame, got {}.'.format(type(tripeptide_comp_seq3)))
        self.assertEqual(tripeptide_comp_seq3.shape, (1, 8000),
            'Expected output to be of shape {}, got {}.'.format((1, 8000), tripeptide_comp_seq3.shape)) 
        self.assertTrue(tripeptide_comp_seq3.iloc[:, : 10].equals(tripeptide_comp_seq3_expected_values), 
            "Output did not match expected values, got\n{}.".format(tripeptide_comp_seq3.iloc[:, : 10]))
        for col in list(tripeptide_comp_seq3.columns):
            #check all columns follow pattern of XY where x & y are amino acids 
            self.assertTrue(bool(re.match(r'^[A-Z]{3}$', col)), 
                "Column doesn't follow correct regex naming convention: {}.".format(col))      
            self.assertIn(col[0], self.amino_acids,
                "Amino acid {} not found in list of amino acids:\n{}.".format(col[0], self.amino_acids))
            self.assertIn(col[1], self.amino_acids,
                "Amino acid {} not found in list of amino acids:\n{}.".format(col[1], self.amino_acids))
            self.assertIn(col[2], self.amino_acids,
                "Amino acid {} not found in list of amino acids:\n{}.".format(col[2], self.amino_acids))
        self.assertTrue(tripeptide_comp_seq3.any().isnull().sum()==0,
            'Expected output to contain no null values.')        
        self.assertTrue(all(col == np.int64 for col in list(tripeptide_comp_seq3.dtypes)),
            "Expected output values to be of datatype np.float64, got {}.".format(list(tripeptide_comp_seq3.dtypes)))
#4.)
        tripeptide_comp_seq4 = protpy.tripeptide_composition(self.protein_seq4)
        #testing first 10 columns 
        tripeptide_comp_seq4_expected_values = pd.DataFrame(columns=["AAA", "AAC", "AAD", "AAE", "AAF", "AAG", "AAH", "AAI", "AAK", "AAL"])
        tripeptide_comp_seq4_expected_values.loc[0] = [0, 0, 1, 1, 0, 0, 0, 1, 0, 1]

        self.assertIsInstance(tripeptide_comp_seq4, pd.DataFrame,
            'Expected output to be of type DataFrame, got {}.'.format(type(tripeptide_comp_seq4)))
        self.assertEqual(tripeptide_comp_seq4.shape, (1, 8000),
            'Expected output to be of shape {}, got {}.'.format((1, 8000), tripeptide_comp_seq4.shape)) 
        self.assertTrue(tripeptide_comp_seq4.iloc[:, : 10].equals(tripeptide_comp_seq4_expected_values), 
            "Output did not match expected values, got\n{}.".format(tripeptide_comp_seq4.iloc[:, : 10]))
        for col in list(tripeptide_comp_seq4.columns):
            #check all columns follow pattern of XY where x & y are amino acids 
            self.assertTrue(bool(re.match(r'^[A-Z]{3}$', col)), 
                "Column doesn't follow correct regex naming convention: {}.".format(col))      
            self.assertIn(col[0], self.amino_acids,
                "Amino acid {} not found in list of amino acids:\n{}.".format(col[0], self.amino_acids))
            self.assertIn(col[1], self.amino_acids,
                "Amino acid {} not found in list of amino acids:\n{}.".format(col[1], self.amino_acids))
            self.assertIn(col[2], self.amino_acids,
                "Amino acid {} not found in list of amino acids:\n{}.".format(col[2], self.amino_acids))
        self.assertTrue(tripeptide_comp_seq4.any().isnull().sum()==0,
            'Expected output to contain no null values.')        
        self.assertTrue(all(col == np.int64 for col in list(tripeptide_comp_seq4.dtypes)),
            "Expected output values to be of datatype np.float64, got {}.".format(list(tripeptide_comp_seq4.dtypes)))
#5.)
        invalid_seq5 = "ABCDEF"
        invalid_seq6 = "OOOOO"
        with (self.assertRaises(ValueError)):
            tripeptide_composition_seq5 = protpy.tripeptide_composition(invalid_seq5)
            tripeptide_composition_seq6 = protpy.tripeptide_composition(invalid_seq6)
#6.)
        invalid_seq7 = 12345
        invalid_seq8 = True
        with (self.assertRaises(TypeError)):
            tripeptide_composition_seq7 = protpy.tripeptide_composition(invalid_seq7)
            tripeptide_composition_seq8 = protpy.tripeptide_composition(invalid_seq8)

    # @unittest.skip("Descriptor can take quite a bit of time to calculate therefore skipping.") #**
    def test_pseudo_amino_acid_composition(self): 
        """ Testing pseudo amino acid composition protein descriptor attributes and functionality. """
        lamda = 30
        weight = 0.05
#1.)
        pseudo_amino_acid_composition_seq1 = protpy.pseudo_amino_acid_composition(self.protein_seq1, weight=weight, lamda=lamda)
        #testing first 10 columns 
        pseudo_amino_acid_composition_seq1_expected_values = pd.DataFrame(columns=["PAAC_1", "PAAC_2", "PAAC_3", "PAAC_4", "PAAC_5", "PAAC_6", "PAAC_7", "PAAC_8", "PAAC_9", "PAAC_10"])
        pseudo_amino_acid_composition_seq1_expected_values.loc[0] = [0.127, 0.059, 0.111, 0.064, 0.126, 0.12, 0.023, 0.118, 0.091, 0.15]
     
        self.assertIsInstance(pseudo_amino_acid_composition_seq1, pd.DataFrame, 
            'Expected output to be of type DataFrame, got {}.'.format(type(pseudo_amino_acid_composition_seq1)))
        self.assertEqual(pseudo_amino_acid_composition_seq1.shape, (1, 20+lamda), 
            'Expected output to be of shape {}, got {}.'.format((1, 20+lamda), pseudo_amino_acid_composition_seq1.shape)) 
        self.assertTrue(pseudo_amino_acid_composition_seq1.iloc[:, : 10].equals(pseudo_amino_acid_composition_seq1_expected_values), 
            "Output did not match expected values, got\n{}.".format(pseudo_amino_acid_composition_seq1.iloc[:, : 10]))
        for col in list(pseudo_amino_acid_composition_seq1.columns):
            #check all columns follow correct naming convention
            self.assertTrue(bool(re.match(r"PAAC_[0-9]", col)), 
                "Column doesn't follow correct regex naming convention: {}.".format(col))      
        self.assertTrue(pseudo_amino_acid_composition_seq1.any().isnull().sum()==0, 
            'Expected output to contain no null values.')        
        self.assertTrue(all(col == np.float64 for col in list(pseudo_amino_acid_composition_seq1.dtypes)),
            "Expected output values to be of datatype np.float64, got {}.".format(list(pseudo_amino_acid_composition_seq1.dtypes)))
#2.)
        pseudo_amino_acid_composition_seq2 = protpy.pseudo_amino_acid_composition(self.protein_seq2, weight=weight, lamda=lamda)
        #testing first 10 columns 
        pseudo_amino_acid_composition_seq2_expected_values = pd.DataFrame(columns=["PAAC_1", "PAAC_2", "PAAC_3", "PAAC_4", "PAAC_5", "PAAC_6", "PAAC_7", "PAAC_8", "PAAC_9", "PAAC_10"])
        pseudo_amino_acid_composition_seq2_expected_values.loc[0] = [0.137, 0.0, 0.088, 0.056, 0.052, 0.181, 0.02, 0.044, 0.117, 0.105]

        self.assertIsInstance(pseudo_amino_acid_composition_seq2, pd.DataFrame, 
            'Expected output to be of type DataFrame, got {}.'.format(type(pseudo_amino_acid_composition_seq2)))
        self.assertEqual(pseudo_amino_acid_composition_seq2.shape, (1, 20+lamda), 
            'Expected output to be of shape {}, got {}.'.format((1, 20+lamda), pseudo_amino_acid_composition_seq2.shape)) 
        self.assertTrue(pseudo_amino_acid_composition_seq2.iloc[:, : 10].equals(pseudo_amino_acid_composition_seq2_expected_values), 
            "Output did not match expected values, got\n{}.".format(pseudo_amino_acid_composition_seq2.iloc[:, : 10]))
        for col in list(pseudo_amino_acid_composition_seq2.columns):
            #check all columns follow correct naming convention
            self.assertTrue(bool(re.match(r"PAAC_[0-9]", col)), 
                "Column doesn't follow correct regex naming convention: {}.".format(col))      
        self.assertTrue(pseudo_amino_acid_composition_seq2.any().isnull().sum()==0, 
            'Expected output to contain no null values.')        
        self.assertTrue(all(col == np.float64 for col in list(pseudo_amino_acid_composition_seq2.dtypes)),
            "Expected output values to be of datatype np.float64, got {}.".format(list(pseudo_amino_acid_composition_seq2.dtypes)))
#3.)
        pseudo_amino_acid_composition_seq3 = protpy.pseudo_amino_acid_composition(self.protein_seq3, weight=weight, lamda=lamda)
        #testing first 10 columns 
        pseudo_amino_acid_composition_seq3_expected_values = pd.DataFrame(columns=["PAAC_1", "PAAC_2", "PAAC_3", "PAAC_4", "PAAC_5", "PAAC_6", "PAAC_7", "PAAC_8", "PAAC_9", "PAAC_10"])
        pseudo_amino_acid_composition_seq3_expected_values.loc[0] = [0.126, 0.095, 0.032, 0.095, 0.126, 0.063, 0.0, 0.095, 0.063, 0.441]

        self.assertIsInstance(pseudo_amino_acid_composition_seq3, pd.DataFrame, 
            'Expected output to be of type DataFrame, got {}.'.format(type(pseudo_amino_acid_composition_seq3)))
        self.assertEqual(pseudo_amino_acid_composition_seq3.shape, (1, 20+lamda), 
            'Expected output to be of shape {}, got {}.'.format((1, 20+lamda), pseudo_amino_acid_composition_seq3.shape)) 
        self.assertTrue(pseudo_amino_acid_composition_seq3.iloc[:, : 10].equals(pseudo_amino_acid_composition_seq3_expected_values), 
            "Output did not match expected values, got\n{}.".format(pseudo_amino_acid_composition_seq3.iloc[:, : 10]))
        for col in list(pseudo_amino_acid_composition_seq3.columns):
            #check all columns follow correct naming convention
            self.assertTrue(bool(re.match(r"PAAC_[0-9]", col)), 
                "Column doesn't follow correct regex naming convention: {}.".format(col))      
        self.assertTrue(pseudo_amino_acid_composition_seq3.any().isnull().sum()==0, 
            'Expected output to contain no null values.')        
        self.assertTrue(all(col == np.float64 for col in list(pseudo_amino_acid_composition_seq3.dtypes)),
            "Expected output values to be of datatype np.float64, got {}.".format(list(pseudo_amino_acid_composition_seq3.dtypes)))
#4.)
        pseudo_amino_acid_composition_seq4 = protpy.pseudo_amino_acid_composition(self.protein_seq4, weight=weight, lamda=lamda)
        #testing first 10 columns 
        pseudo_amino_acid_composition_seq4_expected_values = pd.DataFrame(columns=["PAAC_1", "PAAC_2", "PAAC_3", "PAAC_4", "PAAC_5", "PAAC_6", "PAAC_7", "PAAC_8", "PAAC_9", "PAAC_10"])
        pseudo_amino_acid_composition_seq4_expected_values.loc[0] = [0.15, 0.0, 0.097, 0.048, 0.053, 0.174, 0.016, 0.057, 0.125, 0.109]

        self.assertIsInstance(pseudo_amino_acid_composition_seq4, pd.DataFrame, 
            'Expected output to be of type DataFrame, got {}.'.format(type(pseudo_amino_acid_composition_seq4)))
        self.assertEqual(pseudo_amino_acid_composition_seq4.shape, (1, 20+lamda), 
            'Expected output to be of shape {}, got {}.'.format((1, 20+lamda), pseudo_amino_acid_composition_seq4.shape)) 
        self.assertTrue(pseudo_amino_acid_composition_seq4.iloc[:, : 10].equals(pseudo_amino_acid_composition_seq4_expected_values), 
            "Output did not match expected values, got\n{}.".format(pseudo_amino_acid_composition_seq4.iloc[:, : 10]))
        for col in list(pseudo_amino_acid_composition_seq4.columns):
            #check all columns follow correct naming convention
            self.assertTrue(bool(re.match(r"PAAC_[0-9]", col)), 
                "Column doesn't follow correct regex naming convention: {}.".format(col))      
        self.assertTrue(pseudo_amino_acid_composition_seq4.any().isnull().sum()==0, 
            'Expected output to contain no null values.')        
        self.assertTrue(all(col == np.float64 for col in list(pseudo_amino_acid_composition_seq4.dtypes)),
            "Expected output values to be of datatype np.float64, got {}.".format(list(pseudo_amino_acid_composition_seq4.dtypes)))
#5.)
        invalid_seq5 = "ABCDEF"
        invalid_seq6 = "OOOOO"
        with (self.assertRaises(ValueError)):
            pseudo_amino_acid_composition_seq5 = protpy.pseudo_amino_acid_composition(invalid_seq5)
            pseudo_amino_acid_composition_seq6 = protpy.pseudo_amino_acid_composition(invalid_seq6)
#6.)
        invalid_seq7 = 12345
        invalid_seq8 = False
        with (self.assertRaises(TypeError)):
            pseudo_amino_acid_composition_seq7 = protpy.pseudo_amino_acid_composition(invalid_seq7)
            pseudo_amino_acid_composition_seq8 = protpy.pseudo_amino_acid_composition(invalid_seq8)

    # @unittest.skip("Descriptor can take quite a bit of time to calculate therefore skipping.")
    def test_amphiphilic_pseudo_amino_acid_composition(self):
        """ Testing amphipillic pseudo composition protein descriptor attributes and functionality. """
        lamda = 30
        weight = 0.05
#1.)
        amp_pseudo_amino_acid_composition_seq1 = protpy.amphiphilic_pseudo_amino_acid_composition(self.protein_seq1, weight=weight, lamda=lamda)
        #testing first 10 columns 
        amp_pseudo_amino_acid_composition_seq1_expected_values = pd.DataFrame(columns=["APAAC_1", "APAAC_2", "APAAC_3", "APAAC_4", "APAAC_5", "APAAC_6", "APAAC_7", "APAAC_8", "APAAC_9", "APAAC_10"])
        amp_pseudo_amino_acid_composition_seq1_expected_values.loc[0] = [6.06, 2.814, 5.267, 3.03, 5.988, 5.699, 1.082, 5.627, 4.329, 7.142]

        self.assertIsInstance(amp_pseudo_amino_acid_composition_seq1, pd.DataFrame, 
            'Expected output to be of type DataFrame, got {}.'.format(type(amp_pseudo_amino_acid_composition_seq1)))
        self.assertEqual(amp_pseudo_amino_acid_composition_seq1.shape, (1, 20+(2*lamda)),
            'Expected output to be of shape {}, got {}.'.format((1, 20+(2*lamda)), amp_pseudo_amino_acid_composition_seq1.shape)) 
        self.assertTrue(amp_pseudo_amino_acid_composition_seq1.iloc[:, : 10].equals(amp_pseudo_amino_acid_composition_seq1_expected_values), 
            "Output did not match expected values, got\n{}.".format(amp_pseudo_amino_acid_composition_seq1.iloc[:, : 10]))
        for col in list(amp_pseudo_amino_acid_composition_seq1.columns):
            #check all columns follow correct naming convention
            self.assertTrue(bool(re.match(r"APAAC_[0-9]", col)),
                "Column doesn't follow correct regex naming convention: {}.".format(col))      
        self.assertTrue(amp_pseudo_amino_acid_composition_seq1.any().isnull().sum()==0,
            'Expected output to contain no null values.')        
        self.assertTrue(all(col == np.float64 for col in list(amp_pseudo_amino_acid_composition_seq1.dtypes)),
            "Expected output values to be of datatype np.float64, got {}.".format(list(amp_pseudo_amino_acid_composition_seq1.dtypes)))
#2.)
        amp_pseudo_amino_acid_composition_seq2 = protpy.amphiphilic_pseudo_amino_acid_composition(self.protein_seq2, weight=weight, lamda=lamda)
        #testing first 10 columns 
        amp_pseudo_amino_acid_composition_seq2_expected_values = pd.DataFrame(columns=["APAAC_1", "APAAC_2", "APAAC_3", "APAAC_4", "APAAC_5", "APAAC_6", "APAAC_7", "APAAC_8", "APAAC_9", "APAAC_10"])
        amp_pseudo_amino_acid_composition_seq2_expected_values.loc[0] = [2.867, 0.0, 1.855, 1.181, 1.096, 3.795, 0.422, 0.928, 2.446, 2.193]

        self.assertIsInstance(amp_pseudo_amino_acid_composition_seq2, pd.DataFrame, 
            'Expected output to be of type DataFrame, got {}.'.format(type(amp_pseudo_amino_acid_composition_seq2)))
        self.assertEqual(amp_pseudo_amino_acid_composition_seq2.shape, (1, 20+(2*lamda)),
            'Expected output to be of shape {}, got {}.'.format((1, 20+(2*lamda)), amp_pseudo_amino_acid_composition_seq2.shape)) 
        self.assertTrue(amp_pseudo_amino_acid_composition_seq2.iloc[:, : 10].equals(amp_pseudo_amino_acid_composition_seq2_expected_values), 
            "Output did not match expected values, got\n{}.".format(amp_pseudo_amino_acid_composition_seq2.iloc[:, : 10]))
        for col in list(amp_pseudo_amino_acid_composition_seq2.columns):
            #check all columns follow correct naming convention
            self.assertTrue(bool(re.match(r"APAAC_[0-9]", col)),
                "Column doesn't follow correct regex naming convention: {}.".format(col))      
        self.assertTrue(amp_pseudo_amino_acid_composition_seq2.any().isnull().sum()==0,
            'Expected output to contain no null values.')        
        self.assertTrue(all(col == np.float64 for col in list(amp_pseudo_amino_acid_composition_seq2.dtypes)),
            "Expected output values to be of datatype np.float64, got {}.".format(list(amp_pseudo_amino_acid_composition_seq2.dtypes)))
#3.)
        amp_pseudo_amino_acid_composition_seq3 = protpy.amphiphilic_pseudo_amino_acid_composition(self.protein_seq3, weight=weight, lamda=lamda)
        #testing first 10 columns 
        amp_pseudo_amino_acid_composition_seq3_expected_values = pd.DataFrame(columns=["APAAC_1", "APAAC_2", "APAAC_3", "APAAC_4", "APAAC_5", "APAAC_6", "APAAC_7", "APAAC_8", "APAAC_9", "APAAC_10"])
        amp_pseudo_amino_acid_composition_seq3_expected_values.loc[0] = [1.071, 0.803, 0.268, 0.803, 1.071, 0.536, 0.0, 0.803, 0.536, 3.749]

        self.assertIsInstance(amp_pseudo_amino_acid_composition_seq3, pd.DataFrame, 
            'Expected output to be of type DataFrame, got {}.'.format(type(amp_pseudo_amino_acid_composition_seq3)))
        self.assertEqual(amp_pseudo_amino_acid_composition_seq3.shape, (1, 20+(2*lamda)),
            'Expected output to be of shape {}, got {}.'.format((1, 20+(2*lamda)), amp_pseudo_amino_acid_composition_seq3.shape)) 
        self.assertTrue(amp_pseudo_amino_acid_composition_seq3.iloc[:, : 10].equals(amp_pseudo_amino_acid_composition_seq3_expected_values), 
            "Output did not match expected values, got\n{}.".format(amp_pseudo_amino_acid_composition_seq3.iloc[:, : 10]))
        for col in list(amp_pseudo_amino_acid_composition_seq3.columns):
            #check all columns follow correct naming convention
            self.assertTrue(bool(re.match(r"APAAC_[0-9]", col)),
                "Column doesn't follow correct regex naming convention: {}.".format(col))      
        self.assertTrue(amp_pseudo_amino_acid_composition_seq3.any().isnull().sum()==0,
            'Expected output to contain no null values.')        
        self.assertTrue(all(col == np.float64 for col in list(amp_pseudo_amino_acid_composition_seq3.dtypes)),
            "Expected output values to be of datatype np.float64, got {}.".format(list(amp_pseudo_amino_acid_composition_seq3.dtypes)))
#4.)
        amp_pseudo_amino_acid_composition_seq4 = protpy.amphiphilic_pseudo_amino_acid_composition(self.protein_seq4, weight=weight, lamda=lamda)
        #testing first 10 columns 
        amp_pseudo_amino_acid_composition_seq4_expected_values = pd.DataFrame(columns=["APAAC_1", "APAAC_2", "APAAC_3", "APAAC_4", "APAAC_5", "APAAC_6", "APAAC_7", "APAAC_8", "APAAC_9", "APAAC_10"])
        amp_pseudo_amino_acid_composition_seq4_expected_values.loc[0] = [ 3.209, 0.0, 2.081, 1.041, 1.128, 3.729, 0.347, 1.214, 2.689, 2.342]

        self.assertIsInstance(amp_pseudo_amino_acid_composition_seq4, pd.DataFrame, 
            'Expected output to be of type DataFrame, got {}.'.format(type(amp_pseudo_amino_acid_composition_seq4)))
        self.assertEqual(amp_pseudo_amino_acid_composition_seq4.shape, (1, 20+(2*lamda)),
            'Expected output to be of shape {}, got {}.'.format((1, 20+(2*lamda)), amp_pseudo_amino_acid_composition_seq4.shape)) 
        self.assertTrue(amp_pseudo_amino_acid_composition_seq4.iloc[:, : 10].equals(amp_pseudo_amino_acid_composition_seq4_expected_values), 
            "Output did not match expected values, got\n{}.".format(amp_pseudo_amino_acid_composition_seq4.iloc[:, : 10]))
        for col in list(amp_pseudo_amino_acid_composition_seq4.columns):
            #check all columns follow correct naming convention
            self.assertTrue(bool(re.match(r"APAAC_[0-9]", col)),
                "Column doesn't follow correct regex naming convention: {}.".format(col))      
        self.assertTrue(amp_pseudo_amino_acid_composition_seq4.any().isnull().sum()==0,
            'Expected output to contain no null values.')        
        self.assertTrue(all(col == np.float64 for col in list(amp_pseudo_amino_acid_composition_seq4.dtypes)),
            "Expected output values to be of datatype np.float64, got {}.".format(list(amp_pseudo_amino_acid_composition_seq4.dtypes)))
#5.)
        invalid_seq5 = "ABCDEF"
        invalid_seq6 = "OOOOO"
        with (self.assertRaises(ValueError)):
            amp_pseudo_amino_acid_composition_seq5 = protpy.amphiphilic_pseudo_amino_acid_composition(invalid_seq5)
            amp_pseudo_amino_acid_composition_seq6 = protpy.amphiphilic_pseudo_amino_acid_composition(invalid_seq6)
#6.)
        invalid_seq7 = 12345
        invalid_seq8 = True
        with (self.assertRaises(TypeError)):
            amp_pseudo_amino_acid_composition_seq7 = protpy.amphiphilic_pseudo_amino_acid_composition(invalid_seq7)
            amp_pseudo_amino_acid_composition_seq8 = protpy.amphiphilic_pseudo_amino_acid_composition(invalid_seq8)