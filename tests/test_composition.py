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
    well as the Pseudo and Amphipillic Pseudo Amino Acid descriptors, GRAVY, and 
    the new physicochemical descriptors.

    Test Cases
    ----------
    test_amino_acid_composition:
        testing correct protpy amino acid composition functionality.
    test_gravy:
        testing correct protpy GRAVY (Grand Average of Hydropathy) functionality.
    test_aromaticity:
        testing correct protpy Aromaticity functionality.
    test_instability_index:
        testing correct protpy Instability Index functionality.
    test_isoelectric_point:
        testing correct protpy Isoelectric Point (pI) functionality.
    test_molecular_weight:
        testing correct protpy Molecular Weight functionality.
    test_charge_distribution:
        testing correct protpy Charge Distribution functionality.
    test_hydrophobic_polar_charged_composition:
        testing correct protpy Hydrophobic/Polar/Charged Composition functionality.
    test_secondary_structure_propensity:
        testing correct protpy Secondary Structure Propensity functionality.
    test_kmer_composition:
        testing correct protpy k-mer Composition functionality.
    test_reduced_alphabet_composition:
        testing correct protpy Reduced Alphabet Composition functionality.
    test_motif_composition:
        testing correct protpy Motif Composition functionality.
    test_amino_acid_pair_composition:
        testing correct protpy Amino Acid Pair Composition functionality.
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
        with self.assertRaises(ValueError):
            amino_acid_composition_seq5 = protpy.amino_acid_composition(invalid_seq5)
        with self.assertRaises(ValueError):
            amino_acid_composition_seq6 = protpy.amino_acid_composition(invalid_seq6)
#6.)
        invalid_seq7 = 12345
        invalid_seq8 = True
        with self.assertRaises(TypeError):
            amino_acid_composition_seq7 = protpy.amino_acid_composition(invalid_seq7)
        with self.assertRaises(TypeError):
            amino_acid_composition_seq8 = protpy.amino_acid_composition(invalid_seq8)

    def test_gravy(self):
        """ Testing GRAVY (Grand Average of Hydropathy) descriptor attributes and functionality. """
#1.)
        gravy_seq1 = protpy.gravy(self.protein_seq1)

        self.assertIsInstance(gravy_seq1, pd.DataFrame,
            'Expected output to be of type DataFrame, got {}.'.format(type(gravy_seq1)))
        self.assertEqual(gravy_seq1.shape, (1, 1),
            'Expected output to be of shape {}, got {}.'.format((1, 1), gravy_seq1.shape))
        self.assertEqual(list(gravy_seq1.columns), ["GRAVY"],
            'Incorrect column values found in output:\n{}.'.format(list(gravy_seq1.columns)))
        self.assertTrue(gravy_seq1.any().isnull().sum() == 0,
            'Expected output to contain no null values.')
        self.assertEqual(gravy_seq1["GRAVY"].values[0], -0.045,
            'Expected GRAVY value of -0.045 for seq1, got {}.'.format(gravy_seq1["GRAVY"].values[0]))
#2.)
        gravy_seq2 = protpy.gravy(self.protein_seq2)

        self.assertIsInstance(gravy_seq2, pd.DataFrame,
            'Expected output to be of type DataFrame, got {}.'.format(type(gravy_seq2)))
        self.assertEqual(gravy_seq2.shape, (1, 1),
            'Expected output to be of shape {}, got {}.'.format((1, 1), gravy_seq2.shape))
        self.assertEqual(gravy_seq2["GRAVY"].values[0], -1.027,
            'Expected GRAVY value of -1.027 for seq2, got {}.'.format(gravy_seq2["GRAVY"].values[0]))
#3.)
        gravy_seq3 = protpy.gravy(self.protein_seq3)

        self.assertIsInstance(gravy_seq3, pd.DataFrame,
            'Expected output to be of type DataFrame, got {}.'.format(type(gravy_seq3)))
        self.assertEqual(gravy_seq3.shape, (1, 1),
            'Expected output to be of shape {}, got {}.'.format((1, 1), gravy_seq3.shape))
        self.assertEqual(gravy_seq3["GRAVY"].values[0], 1.141,
            'Expected GRAVY value of 1.141 for seq3, got {}.'.format(gravy_seq3["GRAVY"].values[0]))
#4.)
        gravy_seq4 = protpy.gravy(self.protein_seq4)

        self.assertIsInstance(gravy_seq4, pd.DataFrame,
            'Expected output to be of type DataFrame, got {}.'.format(type(gravy_seq4)))
        self.assertEqual(gravy_seq4.shape, (1, 1),
            'Expected output to be of shape {}, got {}.'.format((1, 1), gravy_seq4.shape))
        self.assertEqual(gravy_seq4["GRAVY"].values[0], -0.971,
            'Expected GRAVY value of -0.971 for seq4, got {}.'.format(gravy_seq4["GRAVY"].values[0]))
#5.)
        invalid_seq5 = "ABCDEF"
        invalid_seq6 = "OOOOO"
        with self.assertRaises(ValueError):
            protpy.gravy(invalid_seq5)
        with self.assertRaises(ValueError):
            protpy.gravy(invalid_seq6)
#6.)
        invalid_seq7 = 12345
        invalid_seq8 = True
        with self.assertRaises(TypeError):
            protpy.gravy(invalid_seq7)
        with self.assertRaises(TypeError):
            protpy.gravy(invalid_seq8)

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
        self.assertEqual(dipeptide_comp_seq4.shape, (1, 400),
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
        with self.assertRaises(ValueError):
            dipeptide_composition_seq5 = protpy.dipeptide_composition(invalid_seq5)
        with self.assertRaises(ValueError):
            dipeptide_composition_seq6 = protpy.dipeptide_composition(invalid_seq6)
#6.)
        invalid_seq7 = 12345
        invalid_seq8 = True
        with self.assertRaises(TypeError):
            dipeptide_composition_seq7 = protpy.dipeptide_composition(invalid_seq7)
        with self.assertRaises(TypeError):
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
        with self.assertRaises(ValueError):
            tripeptide_composition_seq5 = protpy.tripeptide_composition(invalid_seq5)
        with self.assertRaises(ValueError):
            tripeptide_composition_seq6 = protpy.tripeptide_composition(invalid_seq6)
#6.)
        invalid_seq7 = 12345
        invalid_seq8 = True
        with self.assertRaises(TypeError):
            tripeptide_composition_seq7 = protpy.tripeptide_composition(invalid_seq7)
        with self.assertRaises(TypeError):
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
        with self.assertRaises(ValueError):
            pseudo_amino_acid_composition_seq5 = protpy.pseudo_amino_acid_composition(invalid_seq5)
        with self.assertRaises(ValueError):
            pseudo_amino_acid_composition_seq6 = protpy.pseudo_amino_acid_composition(invalid_seq6)
#6.)
        invalid_seq7 = 12345
        invalid_seq8 = False
        with self.assertRaises(TypeError):
            pseudo_amino_acid_composition_seq7 = protpy.pseudo_amino_acid_composition(invalid_seq7)
        with self.assertRaises(TypeError):
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
        amp_pseudo_amino_acid_composition_seq1_expected_values.loc[0] = [6.624, 3.076, 5.757, 3.312, 6.546, 6.23, 1.183, 6.151, 4.732, 7.806]

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
        amp_pseudo_amino_acid_composition_seq2_expected_values.loc[0] = [6.822, 0.0, 4.414, 2.809, 2.609, 9.03, 1.003, 2.207, 5.819, 5.217]

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
        amp_pseudo_amino_acid_composition_seq3_expected_values.loc[0] = [3.783, 2.837, 0.946, 2.837, 3.783, 1.892, 0.0, 2.837, 1.892, 13.24]

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
        amp_pseudo_amino_acid_composition_seq4_expected_values.loc[0] = [7.514, 0.0, 4.874, 2.437, 2.64, 8.733, 0.813, 2.843, 6.296, 5.483]

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
        with self.assertRaises(ValueError):
            amp_pseudo_amino_acid_composition_seq5 = protpy.amphiphilic_pseudo_amino_acid_composition(invalid_seq5)
        with self.assertRaises(ValueError):
            amp_pseudo_amino_acid_composition_seq6 = protpy.amphiphilic_pseudo_amino_acid_composition(invalid_seq6)
#6.)
        invalid_seq7 = 12345
        invalid_seq8 = True
        with self.assertRaises(TypeError):
            amp_pseudo_amino_acid_composition_seq7 = protpy.amphiphilic_pseudo_amino_acid_composition(invalid_seq7)
        with self.assertRaises(TypeError):
            amp_pseudo_amino_acid_composition_seq8 = protpy.amphiphilic_pseudo_amino_acid_composition(invalid_seq8)

    def test_aromaticity(self):
        """ Testing Aromaticity protein descriptor attributes and functionality. """
#1.)
        aromaticity_seq1 = protpy.aromaticity(self.protein_seq1)

        self.assertIsInstance(aromaticity_seq1, pd.DataFrame,
            'Expected output to be of type DataFrame, got {}.'.format(type(aromaticity_seq1)))
        self.assertEqual(aromaticity_seq1.shape, (1, 1),
            'Expected output shape {}, got {}.'.format((1, 1), aromaticity_seq1.shape))
        self.assertEqual(list(aromaticity_seq1.columns), ["Aromaticity"],
            'Incorrect columns: {}.'.format(list(aromaticity_seq1.columns)))
        self.assertTrue(aromaticity_seq1.any().isnull().sum() == 0,
            'Expected no null values.')
        self.assertEqual(aromaticity_seq1["Aromaticity"].values[0], 0.118,
            'Expected 0.118 for seq1, got {}.'.format(aromaticity_seq1["Aromaticity"].values[0]))
#2.)
        aromaticity_seq2 = protpy.aromaticity(self.protein_seq2)
        self.assertEqual(aromaticity_seq2["Aromaticity"].values[0], 0.069,
            'Expected 0.069 for seq2, got {}.'.format(aromaticity_seq2["Aromaticity"].values[0]))
#3.)
        aromaticity_seq3 = protpy.aromaticity(self.protein_seq3)
        self.assertEqual(aromaticity_seq3["Aromaticity"].values[0], 0.105,
            'Expected 0.105 for seq3, got {}.'.format(aromaticity_seq3["Aromaticity"].values[0]))
#4.)
        aromaticity_seq4 = protpy.aromaticity(self.protein_seq4)
        self.assertEqual(aromaticity_seq4["Aromaticity"].values[0], 0.069,
            'Expected 0.069 for seq4, got {}.'.format(aromaticity_seq4["Aromaticity"].values[0]))
#5.)
        with self.assertRaises(ValueError):
            protpy.aromaticity("ABCDEF")
        with self.assertRaises(ValueError):
            protpy.aromaticity("OOOOO")
#6.)
        with self.assertRaises(TypeError):
            protpy.aromaticity(12345)
        with self.assertRaises(TypeError):
            protpy.aromaticity(True)

    def test_instability_index(self):
        """ Testing Instability Index protein descriptor attributes and functionality. """
#1.)
        ii_seq1 = protpy.instability_index(self.protein_seq1)

        self.assertIsInstance(ii_seq1, pd.DataFrame,
            'Expected output to be of type DataFrame, got {}.'.format(type(ii_seq1)))
        self.assertEqual(ii_seq1.shape, (1, 1),
            'Expected output shape {}, got {}.'.format((1, 1), ii_seq1.shape))
        self.assertEqual(list(ii_seq1.columns), ["InstabilityIndex"],
            'Incorrect columns: {}.'.format(list(ii_seq1.columns)))
        self.assertTrue(ii_seq1.any().isnull().sum() == 0,
            'Expected no null values.')
        self.assertEqual(ii_seq1["InstabilityIndex"].values[0], 31.836,
            'Expected 31.836 for seq1, got {}.'.format(ii_seq1["InstabilityIndex"].values[0]))
#2.)
        ii_seq2 = protpy.instability_index(self.protein_seq2)
        self.assertEqual(ii_seq2["InstabilityIndex"].values[0], 57.591,
            'Expected 57.591 for seq2, got {}.'.format(ii_seq2["InstabilityIndex"].values[0]))
#3.)
        ii_seq3 = protpy.instability_index(self.protein_seq3)
        self.assertEqual(ii_seq3["InstabilityIndex"].values[0], 24.357,
            'Expected 24.357 for seq3, got {}.'.format(ii_seq3["InstabilityIndex"].values[0]))
#4.)
        ii_seq4 = protpy.instability_index(self.protein_seq4)
        self.assertEqual(ii_seq4["InstabilityIndex"].values[0], 57.345,
            'Expected 57.345 for seq4, got {}.'.format(ii_seq4["InstabilityIndex"].values[0]))
#5.)
        with self.assertRaises(ValueError):
            protpy.instability_index("ABCDEF")
        with self.assertRaises(ValueError):
            protpy.instability_index("OOOOO")
#6.)
        with self.assertRaises(TypeError):
            protpy.instability_index(12345)
        with self.assertRaises(TypeError):
            protpy.instability_index(True)

    def test_isoelectric_point(self):
        """ Testing Isoelectric Point protein descriptor attributes and functionality. """
#1.)
        pi_seq1 = protpy.isoelectric_point(self.protein_seq1)

        self.assertIsInstance(pi_seq1, pd.DataFrame,
            'Expected output to be of type DataFrame, got {}.'.format(type(pi_seq1)))
        self.assertEqual(pi_seq1.shape, (1, 1),
            'Expected output shape {}, got {}.'.format((1, 1), pi_seq1.shape))
        self.assertEqual(list(pi_seq1.columns), ["IsoelectricPoint"],
            'Incorrect columns: {}.'.format(list(pi_seq1.columns)))
        self.assertTrue(pi_seq1.any().isnull().sum() == 0,
            'Expected no null values.')
        self.assertEqual(pi_seq1["IsoelectricPoint"].values[0], 5.412,
            'Expected 5.412 for seq1, got {}.'.format(pi_seq1["IsoelectricPoint"].values[0]))
#2.)
        pi_seq2 = protpy.isoelectric_point(self.protein_seq2)
        self.assertEqual(pi_seq2["IsoelectricPoint"].values[0], 10.514,
            'Expected 10.514 for seq2, got {}.'.format(pi_seq2["IsoelectricPoint"].values[0]))
#3.)
        pi_seq3 = protpy.isoelectric_point(self.protein_seq3)
        self.assertEqual(pi_seq3["IsoelectricPoint"].values[0], 6.134,
            'Expected 6.134 for seq3, got {}.'.format(pi_seq3["IsoelectricPoint"].values[0]))
#4.)
        pi_seq4 = protpy.isoelectric_point(self.protein_seq4)
        self.assertEqual(pi_seq4["IsoelectricPoint"].values[0], 10.472,
            'Expected 10.472 for seq4, got {}.'.format(pi_seq4["IsoelectricPoint"].values[0]))
#5.)
        with self.assertRaises(ValueError):
            protpy.isoelectric_point("ABCDEF")
        with self.assertRaises(ValueError):
            protpy.isoelectric_point("OOOOO")
#6.)
        with self.assertRaises(TypeError):
            protpy.isoelectric_point(12345)
        with self.assertRaises(TypeError):
            protpy.isoelectric_point(True)

    def test_molecular_weight(self):
        """ Testing Molecular Weight protein descriptor attributes and functionality. """
#1.)
        mw_seq1 = protpy.molecular_weight(self.protein_seq1)

        self.assertIsInstance(mw_seq1, pd.DataFrame,
            'Expected output to be of type DataFrame, got {}.'.format(type(mw_seq1)))
        self.assertEqual(mw_seq1.shape, (1, 1),
            'Expected output shape {}, got {}.'.format((1, 1), mw_seq1.shape))
        self.assertEqual(list(mw_seq1.columns), ["MolecularWeight"],
            'Incorrect columns: {}.'.format(list(mw_seq1.columns)))
        self.assertTrue(mw_seq1.any().isnull().sum() == 0,
            'Expected no null values.')
        self.assertEqual(mw_seq1["MolecularWeight"].values[0], 139122.355,
            'Expected 139122.355 for seq1, got {}.'.format(mw_seq1["MolecularWeight"].values[0]))
#2.)
        mw_seq2 = protpy.molecular_weight(self.protein_seq2)
        self.assertEqual(mw_seq2["MolecularWeight"].values[0], 46023.455,
            'Expected 46023.455 for seq2, got {}.'.format(mw_seq2["MolecularWeight"].values[0]))
#3.)
        mw_seq3 = protpy.molecular_weight(self.protein_seq3)
        self.assertEqual(mw_seq3["MolecularWeight"].values[0], 8360.93,
            'Expected 8360.93 for seq3, got {}.'.format(mw_seq3["MolecularWeight"].values[0]))
#4.)
        mw_seq4 = protpy.molecular_weight(self.protein_seq4)
        self.assertEqual(mw_seq4["MolecularWeight"].values[0], 45624.194,
            'Expected 45624.194 for seq4, got {}.'.format(mw_seq4["MolecularWeight"].values[0]))
#5.)
        with self.assertRaises(ValueError):
            protpy.molecular_weight("ABCDEF")
        with self.assertRaises(ValueError):
            protpy.molecular_weight("OOOOO")
#6.)
        with self.assertRaises(TypeError):
            protpy.molecular_weight(12345)
        with self.assertRaises(TypeError):
            protpy.molecular_weight(True)

    def test_charge_distribution(self):
        """ Testing Charge Distribution protein descriptor attributes and functionality. """
#1.)
        chg_seq1 = protpy.charge_distribution(self.protein_seq1)

        self.assertIsInstance(chg_seq1, pd.DataFrame,
            'Expected output to be of type DataFrame, got {}.'.format(type(chg_seq1)))
        self.assertEqual(chg_seq1.shape, (1, 3),
            'Expected output shape {}, got {}.'.format((1, 3), chg_seq1.shape))
        self.assertEqual(list(chg_seq1.columns), ["PositiveCharge", "NegativeCharge", "NetCharge"],
            'Incorrect columns: {}.'.format(list(chg_seq1.columns)))
        self.assertTrue(chg_seq1.any().isnull().sum() == 0,
            'Expected no null values.')
        self.assertEqual(chg_seq1["PositiveCharge"].values[0], 99.526,
            'Expected PositiveCharge 99.526 for seq1, got {}.'.format(chg_seq1["PositiveCharge"].values[0]))
        self.assertEqual(chg_seq1["NegativeCharge"].values[0], 114.956,
            'Expected NegativeCharge 114.956 for seq1, got {}.'.format(chg_seq1["NegativeCharge"].values[0]))
        self.assertEqual(chg_seq1["NetCharge"].values[0], -15.43,
            'Expected NetCharge -15.43 for seq1, got {}.'.format(chg_seq1["NetCharge"].values[0]))
#2.)
        chg_seq2 = protpy.charge_distribution(self.protein_seq2)
        self.assertEqual(chg_seq2["PositiveCharge"].values[0], 60.168)
        self.assertEqual(chg_seq2["NegativeCharge"].values[0], 35.986)
        self.assertEqual(chg_seq2["NetCharge"].values[0], 24.182)
#3.)
        chg_seq3 = protpy.charge_distribution(self.protein_seq3)
        self.assertEqual(chg_seq3["PositiveCharge"].values[0], 3.998)
        self.assertEqual(chg_seq3["NegativeCharge"].values[0], 3.998)
        self.assertEqual(chg_seq3["NetCharge"].values[0], 0.0)
#4.)
        chg_seq4 = protpy.charge_distribution(self.protein_seq4)
        self.assertEqual(chg_seq4["NetCharge"].values[0], 24.142)
#5.)
        with self.assertRaises(ValueError):
            protpy.charge_distribution("ABCDEF")
        with self.assertRaises(ValueError):
            protpy.charge_distribution("OOOOO")
#6.)
        with self.assertRaises(TypeError):
            protpy.charge_distribution(12345)
        with self.assertRaises(TypeError):
            protpy.charge_distribution(True)

    def test_hydrophobic_polar_charged_composition(self):
        """ Testing Hydrophobic/Polar/Charged Composition descriptor attributes and functionality. """
#1.)
        hpc_seq1 = protpy.hydrophobic_polar_charged_composition(self.protein_seq1)

        self.assertIsInstance(hpc_seq1, pd.DataFrame,
            'Expected output to be of type DataFrame, got {}.'.format(type(hpc_seq1)))
        self.assertEqual(hpc_seq1.shape, (1, 3),
            'Expected output shape {}, got {}.'.format((1, 3), hpc_seq1.shape))
        self.assertEqual(list(hpc_seq1.columns), ["Hydrophobic", "Polar", "Charged"],
            'Incorrect columns: {}.'.format(list(hpc_seq1.columns)))
        self.assertTrue(hpc_seq1.any().isnull().sum() == 0,
            'Expected no null values.')
        self.assertEqual(hpc_seq1["Hydrophobic"].values[0], 44.542)
        self.assertEqual(hpc_seq1["Polar"].values[0], 32.669)
        self.assertEqual(hpc_seq1["Charged"].values[0], 18.247)
#2.)
        hpc_seq2 = protpy.hydrophobic_polar_charged_composition(self.protein_seq2)
        self.assertEqual(hpc_seq2["Hydrophobic"].values[0], 27.962)
        self.assertEqual(hpc_seq2["Polar"].values[0], 40.758)
        self.assertEqual(hpc_seq2["Charged"].values[0], 23.934)
#3.)
        hpc_seq3 = protpy.hydrophobic_polar_charged_composition(self.protein_seq3)
        self.assertEqual(hpc_seq3["Hydrophobic"].values[0], 61.842)
        self.assertEqual(hpc_seq3["Polar"].values[0], 25.0)
        self.assertEqual(hpc_seq3["Charged"].values[0], 10.526)
#4.)
        hpc_seq4 = protpy.hydrophobic_polar_charged_composition(self.protein_seq4)
        self.assertEqual(hpc_seq4["Hydrophobic"].values[0], 29.117)
        self.assertEqual(hpc_seq4["Polar"].values[0], 40.334)
        self.assertEqual(hpc_seq4["Charged"].values[0], 23.866)
#5.)
        with self.assertRaises(ValueError):
            protpy.hydrophobic_polar_charged_composition("ABCDEF")
        with self.assertRaises(ValueError):
            protpy.hydrophobic_polar_charged_composition("OOOOO")
#6.)
        with self.assertRaises(TypeError):
            protpy.hydrophobic_polar_charged_composition(12345)
        with self.assertRaises(TypeError):
            protpy.hydrophobic_polar_charged_composition(True)

    def test_secondary_structure_propensity(self):
        """ Testing Secondary Structure Propensity descriptor attributes and functionality. """
#1.)
        ssp_seq1 = protpy.secondary_structure_propensity(self.protein_seq1)

        self.assertIsInstance(ssp_seq1, pd.DataFrame,
            'Expected output to be of type DataFrame, got {}.'.format(type(ssp_seq1)))
        self.assertEqual(ssp_seq1.shape, (1, 3),
            'Expected output shape {}, got {}.'.format((1, 3), ssp_seq1.shape))
        self.assertEqual(list(ssp_seq1.columns), ["Helix", "Sheet", "Coil"],
            'Incorrect columns: {}.'.format(list(ssp_seq1.columns)))
        self.assertTrue(ssp_seq1.any().isnull().sum() == 0,
            'Expected no null values.')
        self.assertEqual(ssp_seq1["Helix"].values[0], 0.983)
        self.assertEqual(ssp_seq1["Sheet"].values[0], 1.05)
        self.assertEqual(ssp_seq1["Coil"].values[0], 1.043)
#2.)
        ssp_seq2 = protpy.secondary_structure_propensity(self.protein_seq2)
        self.assertEqual(ssp_seq2["Helix"].values[0], 0.954)
        self.assertEqual(ssp_seq2["Sheet"].values[0], 0.952)
        self.assertEqual(ssp_seq2["Coil"].values[0], 1.151)
#3.)
        ssp_seq3 = protpy.secondary_structure_propensity(self.protein_seq3)
        self.assertEqual(ssp_seq3["Helix"].values[0], 1.036)
        self.assertEqual(ssp_seq3["Sheet"].values[0], 1.141)
        self.assertEqual(ssp_seq3["Coil"].values[0], 0.912)
#4.)
        ssp_seq4 = protpy.secondary_structure_propensity(self.protein_seq4)
        self.assertEqual(ssp_seq4["Helix"].values[0], 0.962)
        self.assertEqual(ssp_seq4["Sheet"].values[0], 0.959)
        self.assertEqual(ssp_seq4["Coil"].values[0], 1.138)
#5.)
        with self.assertRaises(ValueError):
            protpy.secondary_structure_propensity("ABCDEF")
        with self.assertRaises(ValueError):
            protpy.secondary_structure_propensity("OOOOO")
#6.)
        with self.assertRaises(TypeError):
            protpy.secondary_structure_propensity(12345)
        with self.assertRaises(TypeError):
            protpy.secondary_structure_propensity(True)

    def test_kmer_composition(self):
        """ Testing k-mer Composition protein descriptor attributes and functionality. """
#1.)
        kmer_seq1 = protpy.kmer_composition(self.protein_seq1, k=2)

        self.assertIsInstance(kmer_seq1, pd.DataFrame,
            'Expected output to be of type DataFrame, got {}.'.format(type(kmer_seq1)))
        self.assertEqual(kmer_seq1.shape, (1, 400),
            'Expected output shape {}, got {}.'.format((1, 400), kmer_seq1.shape))
        self.assertTrue(kmer_seq1.any().isnull().sum() == 0,
            'Expected no null values.')
        self.assertEqual(kmer_seq1["AA"].values[0], 0.797)
        self.assertEqual(kmer_seq1["AC"].values[0], 0.159)
#2.)
        kmer_seq2 = protpy.kmer_composition(self.protein_seq2, k=2)
        self.assertEqual(kmer_seq2.shape, (1, 400))
        self.assertEqual(kmer_seq2["AA"].values[0], 0.713)
#3.)
        kmer_seq3 = protpy.kmer_composition(self.protein_seq3, k=2)
        self.assertEqual(kmer_seq3.shape, (1, 400))
        self.assertEqual(kmer_seq3["AA"].values[0], 0.0)
#4.)
        kmer_seq4 = protpy.kmer_composition(self.protein_seq4, k=2)
        self.assertEqual(kmer_seq4.shape, (1, 400))
        self.assertEqual(kmer_seq4["AA"].values[0], 0.957)
#5.) k=1 should give same shape as amino_acid_composition (1 x 20)
        kmer_k1 = protpy.kmer_composition(self.protein_seq1, k=1)
        self.assertEqual(kmer_k1.shape, (1, 20))
#6.)
        with self.assertRaises(ValueError):
            protpy.kmer_composition("ABCDEF")
        with self.assertRaises(ValueError):
            protpy.kmer_composition("OOOOO")
#7.)
        with self.assertRaises(TypeError):
            protpy.kmer_composition(12345)
        with self.assertRaises(TypeError):
            protpy.kmer_composition(True)

    def test_reduced_alphabet_composition(self):
        """ Testing Reduced Alphabet Composition descriptor attributes and functionality. """
#1.)
        ra_seq1 = protpy.reduced_alphabet_composition(self.protein_seq1)

        self.assertIsInstance(ra_seq1, pd.DataFrame,
            'Expected output to be of type DataFrame, got {}.'.format(type(ra_seq1)))
        self.assertEqual(ra_seq1.shape, (1, 6),
            'Expected output shape {}, got {}.'.format((1, 6), ra_seq1.shape))
        self.assertTrue(ra_seq1.any().isnull().sum() == 0,
            'Expected no null values.')
        self.assertEqual(ra_seq1["ReducedAlphabet_1"].values[0], 25.339)
        self.assertEqual(ra_seq1["ReducedAlphabet_2"].values[0], 34.741)
#2.)
        ra_seq2 = protpy.reduced_alphabet_composition(self.protein_seq2)
        self.assertEqual(ra_seq2.shape, (1, 6))
        self.assertEqual(ra_seq2["ReducedAlphabet_1"].values[0], 24.171)
#3.)
        ra_seq3 = protpy.reduced_alphabet_composition(self.protein_seq3)
        self.assertEqual(ra_seq3.shape, (1, 6))
        self.assertEqual(ra_seq3["ReducedAlphabet_1"].values[0], 25.0)
#4.)
        ra_seq4 = protpy.reduced_alphabet_composition(self.protein_seq4)
        self.assertEqual(ra_seq4.shape, (1, 6))
        self.assertEqual(ra_seq4["ReducedAlphabet_1"].values[0], 25.298)
#5.) alphabet_size=2 should yield 2 columns
        ra_2 = protpy.reduced_alphabet_composition(self.protein_seq1, alphabet_size=2)
        self.assertEqual(ra_2.shape, (1, 2))
#6.)
        with self.assertRaises(ValueError):
            protpy.reduced_alphabet_composition("ABCDEF")
        with self.assertRaises(ValueError):
            protpy.reduced_alphabet_composition("OOOOO")
#7.)
        with self.assertRaises(TypeError):
            protpy.reduced_alphabet_composition(12345)
        with self.assertRaises(TypeError):
            protpy.reduced_alphabet_composition(True)

    def test_motif_composition(self):
        """ Testing Motif Composition descriptor attributes and functionality. """
#1.)
        mc_seq1 = protpy.motif_composition(self.protein_seq1)

        self.assertIsInstance(mc_seq1, pd.DataFrame,
            'Expected output to be of type DataFrame, got {}.'.format(type(mc_seq1)))
        self.assertEqual(mc_seq1.shape, (1, 8),
            'Expected output shape {}, got {}.'.format((1, 8), mc_seq1.shape))
        self.assertTrue(mc_seq1.any().isnull().sum() == 0,
            'Expected no null values.')
        self.assertEqual(mc_seq1["NxS/T_glycosylation"].values[0], 23)
        self.assertEqual(mc_seq1["RGD_integrin"].values[0], 0)
#2.)
        mc_seq2 = protpy.motif_composition(self.protein_seq2)
        self.assertEqual(mc_seq2.shape, (1, 8))
        self.assertEqual(mc_seq2["NxS/T_glycosylation"].values[0], 3)
#3.)
        mc_seq3 = protpy.motif_composition(self.protein_seq3)
        self.assertEqual(mc_seq3.shape, (1, 8))
        self.assertEqual(mc_seq3["NxS/T_glycosylation"].values[0], 2)
        self.assertEqual(mc_seq3["CxxC_zinc_finger"].values[0], 1)
#4.)
        mc_seq4 = protpy.motif_composition(self.protein_seq4)
        self.assertEqual(mc_seq4.shape, (1, 8))
        self.assertEqual(mc_seq4["NxS/T_glycosylation"].values[0], 5)
#5.) custom motifs dict
        custom = protpy.motif_composition(self.protein_seq1, motifs={"RGD": r"RGD", "RXXR": r"R..R"})
        self.assertEqual(custom.shape[1], 2)
        self.assertIn("RGD", custom.columns)
#6.)
        with self.assertRaises(ValueError):
            protpy.motif_composition("ABCDEF")
        with self.assertRaises(ValueError):
            protpy.motif_composition("OOOOO")
#7.)
        with self.assertRaises(TypeError):
            protpy.motif_composition(12345)
        with self.assertRaises(TypeError):
            protpy.motif_composition(True)

    def test_amino_acid_pair_composition(self):
        """ Testing Amino Acid Pair Composition descriptor attributes and functionality. """
#1.)
        pair_seq1 = protpy.amino_acid_pair_composition(self.protein_seq1)

        self.assertIsInstance(pair_seq1, pd.DataFrame,
            'Expected output to be of type DataFrame, got {}.'.format(type(pair_seq1)))
        self.assertEqual(pair_seq1.shape, (1, 400),
            'Expected output shape {}, got {}.'.format((1, 400), pair_seq1.shape))
        self.assertTrue(pair_seq1.any().isnull().sum() == 0,
            'Expected no null values.')
        #check first column label contains physicochemical class annotation
        self.assertIn("_", pair_seq1.columns[0],
            'Expected underscore separator in column names, got {}.'.format(pair_seq1.columns[0]))
        self.assertEqual(pair_seq1["AA_Hydrophobic-Hydrophobic"].values[0], 0.797)
#2.)
        pair_seq2 = protpy.amino_acid_pair_composition(self.protein_seq2)
        self.assertEqual(pair_seq2.shape, (1, 400))
        self.assertEqual(pair_seq2["AA_Hydrophobic-Hydrophobic"].values[0], 0.713)
#3.)
        pair_seq3 = protpy.amino_acid_pair_composition(self.protein_seq3)
        self.assertEqual(pair_seq3.shape, (1, 400))
        self.assertEqual(pair_seq3["AA_Hydrophobic-Hydrophobic"].values[0], 0.0)
#4.)
        pair_seq4 = protpy.amino_acid_pair_composition(self.protein_seq4)
        self.assertEqual(pair_seq4.shape, (1, 400))
        self.assertEqual(pair_seq4["AA_Hydrophobic-Hydrophobic"].values[0], 0.957)
#5.)
        with self.assertRaises(ValueError):
            protpy.amino_acid_pair_composition("ABCDEF")
        with self.assertRaises(ValueError):
            protpy.amino_acid_pair_composition("OOOOO")
#6.)
        with self.assertRaises(TypeError):
            protpy.amino_acid_pair_composition(12345)
        with self.assertRaises(TypeError):
            protpy.amino_acid_pair_composition(True)

    def test_aliphatic_index(self):
        """ Testing Aliphatic Index protein descriptor attributes and functionality. """
#1.)
        ai_seq1 = protpy.aliphatic_index(self.protein_seq1)
        self.assertIsInstance(ai_seq1, pd.DataFrame,
            'Expected output to be of type DataFrame, got {}.'.format(type(ai_seq1)))
        self.assertEqual(ai_seq1.shape, (1, 1),
            'Expected output shape {}, got {}.'.format((1, 1), ai_seq1.shape))
        self.assertEqual(list(ai_seq1.columns), ["AliphaticIndex"],
            'Incorrect columns: {}.'.format(list(ai_seq1.columns)))
        self.assertTrue(ai_seq1.any().isnull().sum() == 0, 'Expected no null values.')
        self.assertEqual(ai_seq1["AliphaticIndex"].values[0], 82.725)
#2.)
        self.assertEqual(protpy.aliphatic_index(self.protein_seq2)["AliphaticIndex"].values[0], 49.81)
#3.)
        self.assertEqual(protpy.aliphatic_index(self.protein_seq3)["AliphaticIndex"].values[0], 145.921)
#4.)
        self.assertEqual(protpy.aliphatic_index(self.protein_seq4)["AliphaticIndex"].values[0], 52.53)
#5.)
        with self.assertRaises(ValueError): protpy.aliphatic_index("ABCDEF")
        with self.assertRaises(ValueError): protpy.aliphatic_index("OOOOO")
#6.)
        with self.assertRaises(TypeError): protpy.aliphatic_index(12345)
        with self.assertRaises(TypeError): protpy.aliphatic_index(True)

    def test_extinction_coefficient(self):
        """ Testing Extinction Coefficient protein descriptor attributes and functionality. """
#1.)
        ec_seq1 = protpy.extinction_coefficient(self.protein_seq1)
        self.assertIsInstance(ec_seq1, pd.DataFrame,
            'Expected output to be of type DataFrame, got {}.'.format(type(ec_seq1)))
        self.assertEqual(ec_seq1.shape, (1, 2),
            'Expected output shape {}, got {}.'.format((1, 2), ec_seq1.shape))
        self.assertEqual(list(ec_seq1.columns), ["ExtCoeff_Reduced", "ExtCoeff_Oxidized"],
            'Incorrect columns: {}.'.format(list(ec_seq1.columns)))
        self.assertTrue(ec_seq1.any().isnull().sum() == 0, 'Expected no null values.')
        self.assertEqual(ec_seq1["ExtCoeff_Reduced"].values[0], 140960)
        self.assertEqual(ec_seq1["ExtCoeff_Oxidized"].values[0], 143335)
#2.)
        ec_seq2 = protpy.extinction_coefficient(self.protein_seq2)
        self.assertEqual(ec_seq2["ExtCoeff_Reduced"].values[0], 43890)
        self.assertEqual(ec_seq2["ExtCoeff_Oxidized"].values[0], 43890)
#3.)
        ec_seq3 = protpy.extinction_coefficient(self.protein_seq3)
        self.assertEqual(ec_seq3["ExtCoeff_Reduced"].values[0], 5960)
        self.assertEqual(ec_seq3["ExtCoeff_Oxidized"].values[0], 6085)
#4.)
        ec_seq4 = protpy.extinction_coefficient(self.protein_seq4)
        self.assertEqual(ec_seq4["ExtCoeff_Reduced"].values[0], 43890)
        self.assertEqual(ec_seq4["ExtCoeff_Oxidized"].values[0], 43890)
#5.)
        with self.assertRaises(ValueError): protpy.extinction_coefficient("ABCDEF")
        with self.assertRaises(ValueError): protpy.extinction_coefficient("OOOOO")
#6.)
        with self.assertRaises(TypeError): protpy.extinction_coefficient(12345)
        with self.assertRaises(TypeError): protpy.extinction_coefficient(True)

    def test_boman_index(self):
        """ Testing Boman Index protein descriptor attributes and functionality. """
#1.)
        bi_seq1 = protpy.boman_index(self.protein_seq1)
        self.assertIsInstance(bi_seq1, pd.DataFrame,
            'Expected output to be of type DataFrame, got {}.'.format(type(bi_seq1)))
        self.assertEqual(bi_seq1.shape, (1, 1),
            'Expected output shape {}, got {}.'.format((1, 1), bi_seq1.shape))
        self.assertEqual(list(bi_seq1.columns), ["BomanIndex"],
            'Incorrect columns: {}.'.format(list(bi_seq1.columns)))
        self.assertTrue(bi_seq1.any().isnull().sum() == 0, 'Expected no null values.')
        self.assertEqual(bi_seq1["BomanIndex"].values[0], 0.119)
#2.)
        self.assertEqual(protpy.boman_index(self.protein_seq2)["BomanIndex"].values[0], -0.498)
#3.)
        self.assertEqual(protpy.boman_index(self.protein_seq3)["BomanIndex"].values[0], 1.237)
#4.)
        self.assertEqual(protpy.boman_index(self.protein_seq4)["BomanIndex"].values[0], -0.455)
#5.)
        with self.assertRaises(ValueError): protpy.boman_index("ABCDEF")
        with self.assertRaises(ValueError): protpy.boman_index("OOOOO")
#6.)
        with self.assertRaises(TypeError): protpy.boman_index(12345)
        with self.assertRaises(TypeError): protpy.boman_index(True)

    def test_aggregation_propensity(self):
        """ Testing Aggregation Propensity protein descriptor attributes and functionality. """
#1.)
        ap_seq1 = protpy.aggregation_propensity(self.protein_seq1)
        self.assertIsInstance(ap_seq1, pd.DataFrame,
            'Expected output to be of type DataFrame, got {}.'.format(type(ap_seq1)))
        self.assertEqual(ap_seq1.shape, (1, 2),
            'Expected output shape {}, got {}.'.format((1, 2), ap_seq1.shape))
        self.assertEqual(list(ap_seq1.columns), ["AggregProneRegions", "AggregProneFraction"],
            'Incorrect columns: {}.'.format(list(ap_seq1.columns)))
        self.assertTrue(ap_seq1.any().isnull().sum() == 0, 'Expected no null values.')
        self.assertEqual(ap_seq1["AggregProneRegions"].values[0], 58)
        self.assertEqual(ap_seq1["AggregProneFraction"].values[0], 11.793)
#2.)
        ap_seq2 = protpy.aggregation_propensity(self.protein_seq2)
        self.assertEqual(ap_seq2["AggregProneRegions"].values[0], 10)
        self.assertEqual(ap_seq2["AggregProneFraction"].values[0], 7.109)
#3.)
        ap_seq3 = protpy.aggregation_propensity(self.protein_seq3)
        self.assertEqual(ap_seq3["AggregProneRegions"].values[0], 20)
        self.assertEqual(ap_seq3["AggregProneFraction"].values[0], 42.105)
#4.)
        ap_seq4 = protpy.aggregation_propensity(self.protein_seq4)
        self.assertEqual(ap_seq4["AggregProneRegions"].values[0], 12)
        self.assertEqual(ap_seq4["AggregProneFraction"].values[0], 7.637)
#5.)
        with self.assertRaises(ValueError): protpy.aggregation_propensity("ABCDEF")
        with self.assertRaises(ValueError): protpy.aggregation_propensity("OOOOO")
#6.)
        with self.assertRaises(TypeError): protpy.aggregation_propensity(12345)
        with self.assertRaises(TypeError): protpy.aggregation_propensity(True)

    def test_hydrophobic_moment(self):
        """ Testing Hydrophobic Moment protein descriptor attributes and functionality. """
#1.)
        hm_seq1 = protpy.hydrophobic_moment(self.protein_seq1)
        self.assertIsInstance(hm_seq1, pd.DataFrame,
            'Expected output to be of type DataFrame, got {}.'.format(type(hm_seq1)))
        self.assertEqual(hm_seq1.shape, (1, 2),
            'Expected output shape {}, got {}.'.format((1, 2), hm_seq1.shape))
        self.assertEqual(list(hm_seq1.columns), ["HydrophobicMoment_Mean", "HydrophobicMoment_Max"],
            'Incorrect columns: {}.'.format(list(hm_seq1.columns)))
        self.assertTrue(hm_seq1.any().isnull().sum() == 0, 'Expected no null values.')
        self.assertEqual(hm_seq1["HydrophobicMoment_Mean"].values[0], 0.272)
        self.assertEqual(hm_seq1["HydrophobicMoment_Max"].values[0], 0.813)
#2.)
        hm_seq2 = protpy.hydrophobic_moment(self.protein_seq2)
        self.assertEqual(hm_seq2["HydrophobicMoment_Mean"].values[0], 0.296)
        self.assertEqual(hm_seq2["HydrophobicMoment_Max"].values[0], 0.760)
#3.)
        hm_seq3 = protpy.hydrophobic_moment(self.protein_seq3)
        self.assertEqual(hm_seq3["HydrophobicMoment_Mean"].values[0], 0.228)
        self.assertEqual(hm_seq3["HydrophobicMoment_Max"].values[0], 0.491)
#4.)
        hm_seq4 = protpy.hydrophobic_moment(self.protein_seq4)
        self.assertEqual(hm_seq4["HydrophobicMoment_Mean"].values[0], 0.294)
        self.assertEqual(hm_seq4["HydrophobicMoment_Max"].values[0], 0.786)
#5.) custom window and angle
        hm_custom = protpy.hydrophobic_moment(self.protein_seq1, window=7, angle=160)
        self.assertEqual(hm_custom.shape, (1, 2))
#6.)
        with self.assertRaises(ValueError): protpy.hydrophobic_moment("ABCDEF")
        with self.assertRaises(ValueError): protpy.hydrophobic_moment("OOOOO")
#7.)
        with self.assertRaises(TypeError): protpy.hydrophobic_moment(12345)
        with self.assertRaises(TypeError): protpy.hydrophobic_moment(True)

    def test_shannon_entropy(self):
        """ Testing Shannon Entropy protein descriptor attributes and functionality. """
#1.)
        se_seq1 = protpy.shannon_entropy(self.protein_seq1)

        self.assertIsInstance(se_seq1, pd.DataFrame,
            'Expected output to be of type DataFrame, got {}.'.format(type(se_seq1)))
        self.assertEqual(se_seq1.shape, (1, 1),
            'Expected output shape {}, got {}.'.format((1, 1), se_seq1.shape))
        self.assertEqual(list(se_seq1.columns), ["ShannonEntropy"],
            'Incorrect columns: {}.'.format(list(se_seq1.columns)))
        self.assertTrue(se_seq1.any().isnull().sum() == 0, 'Expected no null values.')
        self.assertEqual(se_seq1["ShannonEntropy"].values[0], 4.163,
            'Expected ShannonEntropy of 4.163 for seq1, got {}.'.format(se_seq1["ShannonEntropy"].values[0]))
#2.)
        se_seq2 = protpy.shannon_entropy(self.protein_seq2)
        self.assertIsInstance(se_seq2, pd.DataFrame)
        self.assertEqual(se_seq2.shape, (1, 1))
        self.assertEqual(se_seq2["ShannonEntropy"].values[0], 4.024,
            'Expected ShannonEntropy of 4.024 for seq2, got {}.'.format(se_seq2["ShannonEntropy"].values[0]))
#3.)
        se_seq3 = protpy.shannon_entropy(self.protein_seq3)
        self.assertEqual(se_seq3.shape, (1, 1))
        self.assertEqual(se_seq3["ShannonEntropy"].values[0], 3.672,
            'Expected ShannonEntropy of 3.672 for seq3, got {}.'.format(se_seq3["ShannonEntropy"].values[0]))
#4.)
        se_seq4 = protpy.shannon_entropy(self.protein_seq4)
        self.assertEqual(se_seq4.shape, (1, 1))
        self.assertEqual(se_seq4["ShannonEntropy"].values[0], 4.01,
            'Expected ShannonEntropy of 4.01 for seq4, got {}.'.format(se_seq4["ShannonEntropy"].values[0]))
#5.) low-complexity (single amino acid) sequence should have entropy 0
        se_low = protpy.shannon_entropy("AAAAAAAAAA")
        self.assertEqual(se_low["ShannonEntropy"].values[0], 0.0,
            'Expected ShannonEntropy of 0.0 for homopolymer sequence.')
#6.)
        with self.assertRaises(ValueError): protpy.shannon_entropy("ABCDEF")
        with self.assertRaises(ValueError): protpy.shannon_entropy("OOOOO")
#7.)
        with self.assertRaises(TypeError): protpy.shannon_entropy(12345)
        with self.assertRaises(TypeError): protpy.shannon_entropy(True)