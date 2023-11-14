###########################################################################
#############      protPy - Sequence Order Module Tests      ##############
###########################################################################

import pandas as pd
import numpy as np
import os
import unittest
import re
from decimal import *
from Bio import SeqIO
unittest.TestLoader.sortTestMethodsUsing = None

import protpy as protpy

class ProtpySequenceOrderTests(unittest.TestCase):
    """
    Test suite for testing sequence order module and functionality in protpy 
    package, including the sequence order coupling number and quasi sequence
    order descriptors.

    Test Cases
    ----------
    test_sequence_order_coupling_number_:
        testing correct SOCN calculation protpy functionality.
    test_sequence_order_coupling_number:
        testing correct protpy sequence order coupling number functionality
        when generating values for specified distance matrix.
    test_sequence_order_coupling_number_all:
        testing correct protpy sequence order coupling number functionality
        when generating values for both distance matrices.
    test_quasi_sequence_order:
        testing correct protpy quasi sequence order functionality when
        generating values for specified distance matrix.
    test_quasi_sequence_order_all:
        testing correct protpy quasi sequence order functionality when 
        generating values for both distance matrices.
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

        self.distance_matrix_schneider_wrede = "schneider-wrede"
        self.distance_matrix_grantham = "grantham"

    def test_sequence_order_coupling_number_(self):
        """ Testing SOCN function that calculates the individial SOCN value. """
        #expected outputs from SOCN function for 4 test sequences, using both matrices
        socn1_expected_grantham = [401.387, 409.243, 376.946, 393.042, 396.196, 384.772, 383.711, 
            391.951, 386.081, 386.139, 389.241, 393.244, 404.283, 385.307, 379.764]
        socn2_expected_grantham = [140.587, 142.97, 137.385, 139.629, 150.663, 145.954, 148.194, 
            147.217, 153.02, 153.583, 148.249, 151.87, 147.762, 141.598, 148.125]
        socn3_expected_grantham = [17.252, 16.808, 17.018, 16.892, 17.204, 15.861, 16.668, 12.764, 
            16.279, 13.025, 16.401, 14.546, 15.033, 14.816, 12.769]
        socn4_expected_grantham = [138.258, 142.37, 137.803, 139.443, 148.422, 146.01, 146.63, 
            149.337, 150.478, 148.856, 149.214, 154.144, 146.718, 139.508, 149.819]

        socn1_expected_sw = [401.387, 409.243, 376.946, 393.042, 396.196, 384.772, 383.711, 391.951, 
            386.081, 386.139, 389.241, 393.244, 404.283, 385.307, 379.764]
        socn2_expected_sw = [140.587, 142.97, 137.385, 139.629, 150.663, 145.954, 148.194, 147.217, 
            153.02, 153.583, 148.249, 151.87, 147.762, 141.598, 148.125]
        socn3_expected_sw = [17.252, 16.808, 17.018, 16.892, 17.204, 15.861, 16.668, 12.764, 16.279, 
            13.025, 16.401, 14.546, 15.033, 14.816, 12.769]
        socn4_expected_sw = [138.258, 142.37, 137.803, 139.443, 148.422, 146.01, 146.63, 149.337, 
            150.478, 148.856, 149.214, 154.144, 146.718, 139.508, 149.819]

        lag = 15
        for i in range(lag):
#1.)        #testing using schneider wrede distance matrix
            socn1_grantham = protpy.sequence_order_coupling_number_(self.protein_seq1, i+1, self.distance_matrix_schneider_wrede)
            self.assertIsInstance(socn1_grantham, float, 
                'Expected descriptor to be of type float, got {}.'.format(type(socn1_grantham)))
            self.assertEqual(socn1_grantham, socn1_expected_grantham[i], 
                "Output from function doesnt match expected:\n{}.".format(socn1_grantham))
#2.)
            socn2_grantham = protpy.sequence_order_coupling_number_(self.protein_seq2, i+1, self.distance_matrix_schneider_wrede)
            self.assertIsInstance(socn2_grantham, float, 
                'Expected descriptor to be of type float, got {}.'.format(type(socn2_grantham)))
            self.assertEqual(socn2_grantham, socn2_expected_grantham[i], 
                "Output from function doesnt match expected:\n{}.".format(socn2_grantham))
#3.)
            socn3_grantham = protpy.sequence_order_coupling_number_(self.protein_seq3, i+1, self.distance_matrix_schneider_wrede)
            self.assertIsInstance(socn3_grantham, float, 
                'Expected descriptor to be of type float, got {}.'.format(type(socn3_grantham)))
            self.assertEqual(socn3_grantham, socn3_expected_grantham[i], 
                "Output from function doesnt match expected:\n{}.".format(socn3_grantham))
#4.)
            socn4_grantham = protpy.sequence_order_coupling_number_(self.protein_seq4, i+1, self.distance_matrix_schneider_wrede)
            self.assertIsInstance(socn4_grantham, float, 
                'Expected descriptor to be of type float, got {}.'.format(type(socn4_grantham)))
            self.assertEqual(socn4_grantham, socn4_expected_grantham[i], 
                "Output from function doesnt match expected:\n{}.".format(socn4_grantham))
#5.)        #testing using grantham distance matrix
            socn1_sw = protpy.sequence_order_coupling_number_(self.protein_seq1, i+1, self.distance_matrix_grantham)
            self.assertIsInstance(socn1_sw, float, 
                'Expected descriptor to be of type float, got {}.'.format(type(socn1_sw)))
            self.assertEqual(socn1_sw, socn1_expected_sw[i], 
                "Output from function doesnt match expected:\n{}.".format(socn1_sw))
#6.)
            socn2_sw = protpy.sequence_order_coupling_number_(self.protein_seq2, i+1, self.distance_matrix_grantham)
            self.assertIsInstance(socn2_sw, float, 
                'Expected descriptor to be of type float, got {}.'.format(type(socn2_sw)))
            self.assertEqual(socn2_sw, socn2_expected_sw[i], 
                "Output from function doesnt match expected:\n{}.".format(socn2_sw))
#7.)
            socn3_sw = protpy.sequence_order_coupling_number_(self.protein_seq3, i+1, self.distance_matrix_grantham)
            self.assertIsInstance(socn3_sw, float, 
                'Expected descriptor to be of type float, got {}.'.format(type(socn3_sw)))
            self.assertEqual(socn3_sw, socn3_expected_sw[i], 
                "Output from function doesnt match expected:\n{}.".format(socn3_sw))
#8.)
            socn4_sw = protpy.sequence_order_coupling_number_(self.protein_seq4, i+1, self.distance_matrix_grantham)
            self.assertIsInstance(socn4_sw, float, 
                'Expected descriptor to be of type float, got {}.'.format(type(socn4_sw)))
            self.assertEqual(socn4_sw, socn4_expected_sw[i], 
                "Output from function doesnt match expected:\n{}.".format(socn4_sw))
#9.)
            invalid_seq5 = "ABCDEF"
            invalid_seq6 = "OOOOO"
            with (self.assertRaises(ValueError)):
                socn_seq5 = protpy.sequence_order_coupling_number_(invalid_seq5, i+1, 
                    distance_matrix=self.distance_matrix_schneider_wrede)
                socn_seq6 = protpy.sequence_order_coupling_number_(invalid_seq6, i+1, 
                    distance_matrix=self.distance_matrix_grantham)
#10.)
            invalid_seq7 = 1234
            invalid_seq8 = False
            with (self.assertRaises(TypeError)):
                socn_seq7 = protpy.sequence_order_coupling_number_(invalid_seq7, i+1, 
                    distance_matrix=self.distance_matrix_grantham)
                socn_seq8 = protpy.sequence_order_coupling_number_(invalid_seq8, i+1, 
                    distance_matrix=self.distance_matrix_grantham)

    def test_sequence_order_coupling_number(self):
        """ Testing sequence order coupling number descriptor attributes and functionality. """
#1.)
        lag=5
        socn_seq1 = protpy.sequence_order_coupling_number(self.protein_seq1, lag=lag, 
            distance_matrix=self.distance_matrix_schneider_wrede)
        #testing first 10 columns 
        socn_seq1_expected_values = pd.DataFrame(columns=["SOCN_SW1", "SOCN_SW2", "SOCN_SW3", "SOCN_SW4", "SOCN_SW5"])
        socn_seq1_expected_values.loc[0] = [401.387, 409.243, 376.946, 393.042, 396.196]

        self.assertIsInstance(socn_seq1, pd.DataFrame, 
            'Expected output to be a DataFrame, got {}.'.format(type(socn_seq1)))
        self.assertEqual(socn_seq1.shape, (1, lag), 
            'Expected output to be of shape 1 x 5, got {}.'.format(socn_seq1.shape)) 
        self.assertTrue(socn_seq1.iloc[:, : 5].equals(socn_seq1_expected_values), 
            "Output did not match expected values, got\n{}.".format(socn_seq1.iloc[:, : 10]))
        for col in list(socn_seq1.columns):
            #check all columns follow pattern of SOCNX or SOCNXY where x & y integers between 0 and 9
            self.assertTrue((bool(re.match(r'SOCN_SW[0-9]', col)) or bool(re.match(r'SOCN_SW[0-9][0-9]', col))), 
                "Column name doesn't match expected regex pattern: {}.".format(col))     
        self.assertTrue(socn_seq1.any().isnull().sum()==0, 
            'Expected output to contain no null values.')        
        self.assertTrue(all(col == np.float64 for col in list(socn_seq1.dtypes)), 
            "Descriptor dataframe dtypes expected to be np.float64, got {}.".format(list(socn_seq1.dtypes)))
#2.)
        lag=10
        socn_seq2 = protpy.sequence_order_coupling_number(self.protein_seq2, lag=lag, 
            distance_matrix=self.distance_matrix_schneider_wrede)
        #testing first 10 columns 
        socn_seq2_expected_values = pd.DataFrame(columns=["SOCN_SW1", "SOCN_SW2", "SOCN_SW3", "SOCN_SW4", "SOCN_SW5", "SOCN_SW6", "SOCN_SW7", "SOCN_SW8", "SOCN_SW9", "SOCN_SW10"])
        socn_seq2_expected_values.loc[0] = [140.587, 142.97, 137.385, 139.629, 150.663, 145.954, 148.194, 147.217, 153.02, 153.583]

        self.assertIsInstance(socn_seq2, pd.DataFrame, 
            'Expected output to be a DataFrame, got {}.'.format(type(socn_seq2)))
        self.assertEqual(socn_seq2.shape, (1, lag), 
            'Expected output to be of shape 1 x 10, got {}.'.format(socn_seq2.shape)) 
        self.assertTrue(socn_seq2.iloc[:, : 10].equals(socn_seq2_expected_values), 
            "Output did not match expected values, got\n{}.".format(socn_seq2.iloc[:, : 10]))
        for col in list(socn_seq2.columns):
            #check all columns follow pattern of SOCNX or SOCNXY where x & y integers between 0 and 9
            self.assertTrue((bool(re.match(r'SOCN_SW[0-9]', col)) or bool(re.match(r'SOCN_SW[0-9][0-9]', col))), 
                "Column name doesn't match expected regex pattern: {}.".format(col))     
        self.assertTrue(socn_seq2.any().isnull().sum()==0, 
            'Expected output to contain no null values.')        
        self.assertTrue(all(col == np.float64 for col in list(socn_seq2.dtypes)), 
            "Descriptor dataframe dtypes expected to be np.float64, got {}.".format(list(socn_seq2.dtypes)))
#3.)
        lag=15
        socn_seq3 = protpy.sequence_order_coupling_number(self.protein_seq3, lag=lag, 
            distance_matrix=self.distance_matrix_schneider_wrede)
        #testing first 10 columns 
        socn_seq3_expected_values = pd.DataFrame(columns=["SOCN_SW1", "SOCN_SW2", "SOCN_SW3", "SOCN_SW4", "SOCN_SW5", "SOCN_SW6", "SOCN_SW7", "SOCN_SW8", "SOCN_SW9", "SOCN_SW10"])
        socn_seq3_expected_values.loc[0] = [17.252, 16.808, 17.018, 16.892, 17.204, 15.861, 16.668, 12.764, 16.279, 13.025]

        self.assertIsInstance(socn_seq3, pd.DataFrame, 
            'Expected output to be a DataFrame, got {}.'.format(type(socn_seq3)))
        self.assertEqual(socn_seq3.shape, (1, lag), 
            'Expected output to be of shape 1 x 15, got {}.'.format(socn_seq3.shape)) 
        self.assertTrue(socn_seq3.iloc[:, : 10].equals(socn_seq3_expected_values), 
            "Output did not match expected values, got\n{}.".format(socn_seq3.iloc[:, : 10]))
        for col in list(socn_seq3.columns):
            #check all columns follow pattern of SOCNX or SOCNXY where x & y integers between 0 and 9
            self.assertTrue((bool(re.match(r'SOCN_SW[0-9]', col)) or bool(re.match(r'SOCN_SW[0-9][0-9]', col))), 
                "Column name doesn't match expected regex pattern: {}.".format(col))     
        self.assertTrue(socn_seq3.any().isnull().sum()==0, 
            'Expected output to contain no null values.')        
        self.assertTrue(all(col == np.float64 for col in list(socn_seq3.dtypes)), 
            "Descriptor dataframe dtypes expected to be np.float64, got {}.".format(list(socn_seq3.dtypes)))
#4.)
        invalid_lag = 500 #greater than len of protein sequence
        socn_seq4 = protpy.sequence_order_coupling_number(self.protein_seq4, lag=invalid_lag, 
            distance_matrix=self.distance_matrix_schneider_wrede)
        #testing first 10 columns 
        socn_seq4_expected_values = pd.DataFrame(columns=["SOCN_SW1", "SOCN_SW2", "SOCN_SW3", "SOCN_SW4", "SOCN_SW5", "SOCN_SW6", "SOCN_SW7", "SOCN_SW8", "SOCN_SW9", "SOCN_SW10"])
        socn_seq4_expected_values.loc[0] = [138.258, 142.37, 137.803, 139.443, 148.422, 146.01, 146.63, 149.337, 150.478, 148.856]

        self.assertIsInstance(socn_seq4, pd.DataFrame, 
            'Expected output to be a DataFrame, got {}.'.format(type(socn_seq4)))
        self.assertEqual(socn_seq4.shape, (1, 30), 
            'Expected output to be of shape 1 x 30, got {}.'.format(socn_seq4.shape)) 
        self.assertTrue(socn_seq4.iloc[:, : 10].equals(socn_seq4_expected_values), 
            "Output did not match expected values, got\n{}.".format(socn_seq4.iloc[:, : 10]))
        for col in list(socn_seq4.columns):
            #check all columns follow pattern of SOCNX or SOCNXY where x & y integers between 0 and 9
            self.assertTrue((bool(re.match(r'SOCN_SW[0-9]', col)) or bool(re.match(r'SOCN_SW[0-9][0-9]', col))), 
                "Column name doesn't match expected regex pattern: {}.".format(col))     
        self.assertTrue(socn_seq4.any().isnull().sum()==0, 
            'Expected output to contain no null values.')        
        self.assertTrue(all(col == np.float64 for col in list(socn_seq4.dtypes)), 
            "Descriptor dataframe dtypes expected to be np.float64, got {}.".format(list(socn_seq4.dtypes)))
#5.)
        invalid_seq5 = "ABCDEF"
        invalid_seq6 = "OOOOO"
        invalid_lag = "BLAHBLAH"
        with (self.assertRaises(ValueError)):
            socn_seq5 = protpy.sequence_order_coupling_number(invalid_seq5, lag=10, 
                distance_matrix=self.distance_matrix_schneider_wrede)
            socn_seq6 = protpy.sequence_order_coupling_number(invalid_seq6, lag=10, 
                distance_matrix=self.distance_matrix_schneider_wrede)
            socn_seq7 = protpy.sequence_order_coupling_number(self.protein_seq4, lag=invalid_lag, 
                distance_matrix=self.distance_matrix_schneider_wrede)
#6.)
        invalid_seq7 = 12345
        invalid_seq8 = False
        with (self.assertRaises(TypeError)):
            socn_seq7 = protpy.sequence_order_coupling_number(invalid_seq7, lag=10, 
                distance_matrix=self.distance_matrix_schneider_wrede)
            socn_seq8 = protpy.sequence_order_coupling_number(invalid_seq8, lag=10, 
                distance_matrix=self.distance_matrix_schneider_wrede)

    def test_sequence_order_coupling_number_all(self):
        """ Testing function that calculates SOCN using both matrices for a sequence and lag. """
#1.)
        lag=30
        socn_all_seq1 = protpy.sequence_order_coupling_number_all(self.protein_seq1, lag=lag)
        #testing first 10 columns 
        socn_all_seq1_expected_values = pd.DataFrame(columns=["SOCN_SW1", "SOCN_SW2", "SOCN_SW3", "SOCN_SW4", "SOCN_SW5", "SOCN_SW6", "SOCN_SW7", "SOCN_SW8", "SOCN_SW9", "SOCN_SW10"])
        socn_all_seq1_expected_values.loc[0] = [401.387, 409.243, 376.946, 393.042, 396.196, 384.772, 383.711, 391.951, 386.081, 386.139]

        self.assertIsInstance(socn_all_seq1, pd.DataFrame, 
            'Expected output to be a DataFrame, got {}.'.format(type(socn_all_seq1)))
        self.assertEqual(socn_all_seq1.shape, (1, lag*2), 
            'Expected output to be of shape {}, got {}.'.format((1, lag*2), socn_all_seq1.shape))
        self.assertTrue(socn_all_seq1.iloc[:, : 10].equals(socn_all_seq1_expected_values), 
            "Output did not match expected values, got\n{}.".format(socn_all_seq1.iloc[:, : 10]))
        for col in list(socn_all_seq1.columns):
            #check all columns follow pattern of SOCNX or SOCNXY where x & y integers between 0 and 9
            self.assertTrue((bool(re.match(r'SOCN_SW[0-9]', col)) or bool(re.match(r'SOCN_SW[0-9][0-9]', col)) or
                bool(re.match(r'SOCN_Grant[0-9]', col)) or bool(re.match(r'SOCN_Grant[0-9][0-9]', col))), 
                "Column name doesn't match expected regex pattern: {}.".format(col))     
        self.assertTrue(socn_all_seq1.any().isnull().sum()==0, 
            'Expected output to contain no null values.')        
        self.assertTrue(all(col == np.float64 for col in list(socn_all_seq1.dtypes)), 
            "Descriptor dataframe dtypes expected to be np.float64, got {}.".format(list(socn_all_seq1.dtypes)))
#2.)
        lag=15
        socn_all_seq2 = protpy.sequence_order_coupling_number_all(self.protein_seq2, lag=lag)
        #testing first 10 columns 
        socn_all_seq2_expected_values = pd.DataFrame(columns=["SOCN_SW1", "SOCN_SW2", "SOCN_SW3", "SOCN_SW4", "SOCN_SW5", "SOCN_SW6", "SOCN_SW7", "SOCN_SW8", "SOCN_SW9", "SOCN_SW10"])
        socn_all_seq2_expected_values.loc[0] = [140.587, 142.97, 137.385, 139.629, 150.663, 145.954, 148.194, 147.217, 153.02, 153.583]

        self.assertIsInstance(socn_all_seq2, pd.DataFrame, 
            'Expected output to be a DataFrame, got {}.'.format(type(socn_all_seq2)))
        self.assertEqual(socn_all_seq2.shape, (1, lag*2), 
            'Expected output to be of shape {}, got {}.'.format((1, lag*2), socn_all_seq2.shape))
        self.assertTrue(socn_all_seq2.iloc[:, : 10].equals(socn_all_seq2_expected_values), 
            "Output did not match expected values, got\n{}.".format(socn_all_seq2.iloc[:, : 10]))
        for col in list(socn_all_seq2.columns):
            #check all columns follow pattern of SOCNX or SOCNXY where x & y integers between 0 and 9
            self.assertTrue((bool(re.match(r'SOCN_SW[0-9]', col)) or bool(re.match(r'SOCN_SW[0-9][0-9]', col)) or
                bool(re.match(r'SOCN_Grant[0-9]', col)) or bool(re.match(r'SOCN_Grant[0-9][0-9]', col))), 
                "Column name doesn't match expected regex pattern: {}.".format(col))     
        self.assertTrue(socn_all_seq2.any().isnull().sum()==0, 
            'Expected output to contain no null values.')        
        self.assertTrue(all(col == np.float64 for col in list(socn_all_seq2.dtypes)), 
            "Descriptor dataframe dtypes expected to be np.float64, got {}.".format(list(socn_all_seq2.dtypes)))
#3.)
        lag=10
        socn_all_seq3 = protpy.sequence_order_coupling_number_all(self.protein_seq3, lag=lag)
        #testing first 10 columns 
        socn_all_seq3_expected_values = pd.DataFrame(columns=["SOCN_SW1", "SOCN_SW2", "SOCN_SW3", "SOCN_SW4", "SOCN_SW5", "SOCN_SW6", "SOCN_SW7", "SOCN_SW8", "SOCN_SW9", "SOCN_SW10"])
        socn_all_seq3_expected_values.loc[0] = [17.252, 16.808, 17.018, 16.892, 17.204, 15.861, 16.668, 12.764, 16.279, 13.025]

        self.assertIsInstance(socn_all_seq3, pd.DataFrame, 
            'Expected output to be a DataFrame, got {}.'.format(type(socn_all_seq3)))
        self.assertEqual(socn_all_seq3.shape, (1, lag*2), 
            'Expected output to be of shape {}, got {}.'.format((1, lag*2), socn_all_seq3.shape)) 
        self.assertTrue(socn_all_seq3.iloc[:, : 10].equals(socn_all_seq3_expected_values), 
            "Output did not match expected values, got\n{}.".format(socn_all_seq3.iloc[:, : 10]))
        for col in list(socn_all_seq3.columns):
            #check all columns follow pattern of SOCNX or SOCNXY where x & y integers between 0 and 9
            self.assertTrue((bool(re.match(r'SOCN_SW[0-9]', col)) or bool(re.match(r'SOCN_SW[0-9][0-9]', col)) or
                bool(re.match(r'SOCN_Grant[0-9]', col)) or bool(re.match(r'SOCN_Grant[0-9][0-9]', col))), 
                "Column name doesn't match expected regex pattern: {}.".format(col))     
        self.assertTrue(socn_all_seq3.any().isnull().sum()==0, 
            'Expected output to contain no null values.')        
        self.assertTrue(all(col == np.float64 for col in list(socn_all_seq3.dtypes)), 
            "Descriptor dataframe dtypes expected to be np.float64, got {}.".format(list(socn_all_seq3.dtypes)))
#4.)
        lag=5
        socn_all_seq4 = protpy.sequence_order_coupling_number_all(self.protein_seq4, lag=lag)
        #testing first 10 columns 
        socn_all_seq4_expected_values = pd.DataFrame(columns=["SOCN_SW1", "SOCN_SW2", "SOCN_SW3", "SOCN_SW4", "SOCN_SW5", "SOCN_Grant1", "SOCN_Grant2", "SOCN_Grant3", "SOCN_Grant4", "SOCN_Grant5"])
        socn_all_seq4_expected_values.loc[0] = [138.258, 142.37, 137.803, 139.443, 148.422, 138.258, 142.37, 137.803, 139.443, 148.422]

        self.assertIsInstance(socn_all_seq4, pd.DataFrame, 
            'Expected output to be a DataFrame, got {}.'.format(type(socn_all_seq4)))
        self.assertEqual(socn_all_seq4.shape, (1, lag*2), 
            'Expected output to be of shape {}, got {}.'.format((1, lag*2), socn_all_seq4.shape))
        self.assertTrue(socn_all_seq4.iloc[:, : 10].equals(socn_all_seq4_expected_values), 
            "Output did not match expected values, got\n{}.".format(socn_all_seq4.iloc[:, : 10]))
        for col in list(socn_all_seq4.columns):
            #check all columns follow pattern of SOCNX or SOCNXY where x & y integers between 0 and 9
            self.assertTrue((bool(re.match(r'SOCN_SW[0-9]', col)) or bool(re.match(r'SOCN_SW[0-9][0-9]', col)) or
                bool(re.match(r'SOCN_Grant[0-9]', col)) or bool(re.match(r'SOCN_Grant[0-9][0-9]', col))), 
                "Column name doesn't match expected regex pattern: {}.".format(col))     
        self.assertTrue(socn_all_seq4.any().isnull().sum()==0, 
            'Expected output to contain no null values.')        
        self.assertTrue(all(col == np.float64 for col in list(socn_all_seq4.dtypes)), 
            "Descriptor dataframe dtypes expected to be np.float64, got {}.".format(list(socn_all_seq4.dtypes)))
#5.)
        invalid_seq5 = "ABCDEF"
        invalid_seq6 = "OOOOO"
        with (self.assertRaises(ValueError)):
            socn_all_seq5 = protpy.sequence_order_coupling_number_all(invalid_seq5, lag=lag)
            socn_all_seq6 = protpy.sequence_order_coupling_number_all(invalid_seq6, lag=lag)
#6.)
        invalid_seq7 = 1234
        invalid_lag = "BLAHBLAH"
        with (self.assertRaises(TypeError)):
            socn_all_seq7 = protpy.sequence_order_coupling_number_all(invalid_seq7, lag=lag)
            socn_all_seq8 = protpy.sequence_order_coupling_number_all(self.protein_seq4, lag=invalid_lag)

    def test_quasi_sequence_order(self):
        """ Testing quasi sequence order descriptor attributes and functionality. """
#1.)    #testing using distance_matrix_schneider_wrede distance matrix
        lag = 30
        quasi_sequence_order_seq1 = protpy.quasi_sequence_order(self.protein_seq1, lag=lag, 
            distance_matrix=self.distance_matrix_schneider_wrede)
        #testing first 10 columns 
        quasi_sequence_order_seq1_expected_values = pd.DataFrame(columns=["QSO_SW1", "QSO_SW2", "QSO_SW3", "QSO_SW4", "QSO_SW5", "QSO_SW6", "QSO_SW7", "QSO_SW8", "QSO_SW9", "QSO_SW10"])
        quasi_sequence_order_seq1_expected_values.loc[0] = [0.005692, 0.002643, 0.004947, 0.002846, 0.005625, 0.005354, 0.001016, 0.005286, 0.004066, 0.006708]

        self.assertIsInstance(quasi_sequence_order_seq1, pd.DataFrame,
            'Expected output to be a DataFrame, got {}.'.format(type(quasi_sequence_order_seq1)))
        self.assertEqual(quasi_sequence_order_seq1.shape, (1, lag+20), 
            'Expected output to be of shape 1 x 50, got {}.'.format(quasi_sequence_order_seq1.shape))
        self.assertTrue(quasi_sequence_order_seq1.iloc[:, : 10].equals(quasi_sequence_order_seq1_expected_values), 
            "Output did not match expected values, got\n{}.".format(quasi_sequence_order_seq1.iloc[:, : 10]))
        for col in list(quasi_sequence_order_seq1.columns):
            #check all columns follow pattern of SOCNX or SOCNXY where x & y integers between 0 and 9
            self.assertTrue((bool(re.match(r'QSO_SW[0-9]', col)) or bool(re.match(r'QSO_SW[0-9][0-9]', col))), 
                "Column name doesn't match expected regex pattern: {}.".format(col))     
        self.assertTrue(quasi_sequence_order_seq1.any().isnull().sum()==0, 
            'Expected output to contain no null values.')        
        self.assertTrue(all(col == np.float64 for col in list(quasi_sequence_order_seq1.dtypes)), 
            "Descriptor dataframe dtypes expected to be np.float64, got {}.".format(list(quasi_sequence_order_seq1.dtypes)))
#2.)
        lag = 20
        quasi_sequence_order_seq2 = protpy.quasi_sequence_order(self.protein_seq2, lag=lag, 
            distance_matrix=self.distance_matrix_schneider_wrede)
        #testing first 10 columns 
        quasi_sequence_order_seq2_expected_values = pd.DataFrame(columns=["QSO_SW1", "QSO_SW2", "QSO_SW3", "QSO_SW4", "QSO_SW5", "QSO_SW6", "QSO_SW7", "QSO_SW8", "QSO_SW9", "QSO_SW10"])
        quasi_sequence_order_seq2_expected_values.loc[0] = [0.027474, 0.0, 0.017776, 0.011314, 0.010506, 0.036364, 0.004041, 0.00889, 0.023434, 0.021009]

        self.assertIsInstance(quasi_sequence_order_seq2, pd.DataFrame,
            'Expected output to be a DataFrame, got {}.'.format(type(quasi_sequence_order_seq2)))
        self.assertEqual(quasi_sequence_order_seq2.shape, (1, lag+20), 
            'Expected output to be of shape 1 x 40, got {}.'.format(quasi_sequence_order_seq2.shape)) 
        self.assertTrue(quasi_sequence_order_seq2.iloc[:, : 10].equals(quasi_sequence_order_seq2_expected_values), 
            "Output did not match expected values, got\n{}.".format(quasi_sequence_order_seq2.iloc[:, : 10]))
        for col in list(quasi_sequence_order_seq2.columns):
            #check all columns follow pattern of SOCNX or SOCNXY where x & y integers between 0 and 9
            self.assertTrue((bool(re.match(r'QSO_SW[0-9]', col)) or bool(re.match(r'QSO_SW[0-9][0-9]', col))), 
                "Column name doesn't match expected regex pattern: {}.".format(col))     
        self.assertTrue(quasi_sequence_order_seq2.any().isnull().sum()==0, 
            'Expected output to contain no null values.')        
        self.assertTrue(all(col == np.float64 for col in list(quasi_sequence_order_seq2.dtypes)), 
            "Descriptor dataframe dtypes expected to be np.float64, got {}.".format(list(quasi_sequence_order_seq2.dtypes)))
#3.)
        lag = 10
        quasi_sequence_order_seq3 = protpy.quasi_sequence_order(self.protein_seq3, lag=lag, 
            distance_matrix=self.distance_matrix_schneider_wrede)
        #testing first 10 columns 
        quasi_sequence_order_seq3_expected_values = pd.DataFrame(columns=["QSO_SW1", "QSO_SW2", "QSO_SW3", "QSO_SW4", "QSO_SW5", "QSO_SW6", "QSO_SW7", "QSO_SW8", "QSO_SW9", "QSO_SW10"])
        quasi_sequence_order_seq3_expected_values.loc[0] = [0.310006, 0.23249, 0.077516, 0.23249, 0.310006, 0.155032, 0.0, 0.23249, 0.155032, 1.08505]

        self.assertIsInstance(quasi_sequence_order_seq3, pd.DataFrame,
            'Expected output to be a DataFrame, got {}.'.format(type(quasi_sequence_order_seq3)))
        self.assertEqual(quasi_sequence_order_seq3.shape, (1, lag+20), 
            'Expected output to be of shape 1 x 30, got {}.'.format(quasi_sequence_order_seq3.shape)) 
        self.assertTrue(quasi_sequence_order_seq3.iloc[:, : 10].equals(quasi_sequence_order_seq3_expected_values), 
            "Output did not match expected values, got\n{}.".format(quasi_sequence_order_seq3.iloc[:, : 10]))
        for col in list(quasi_sequence_order_seq3.columns):
            #check all columns follow pattern of SOCNX or SOCNXY where x & y integers between 0 and 9
            self.assertTrue((bool(re.match(r'QSO_SW[0-9]', col)) or bool(re.match(r'QSO_SW[0-9][0-9]', col))), 
                "Column name doesn't match expected regex pattern: {}.".format(col))     
        self.assertTrue(quasi_sequence_order_seq3.any().isnull().sum()==0, 
            'Expected output to contain no null values.')        
        self.assertTrue(all(col == np.float64 for col in list(quasi_sequence_order_seq3.dtypes)), 
            "Descriptor dataframe dtypes expected to be np.float64, got {}.".format(list(quasi_sequence_order_seq3.dtypes)))
#4.)
        lag = 5
        quasi_sequence_order_seq4 = protpy.quasi_sequence_order(self.protein_seq4, lag=lag, 
            distance_matrix=self.distance_matrix_schneider_wrede)
        #testing first 10 columns 
        quasi_sequence_order_seq4_expected_values = pd.DataFrame(columns=["QSO_SW1", "QSO_SW2", "QSO_SW3", "QSO_SW4", "QSO_SW5", "QSO_SW6", "QSO_SW7", "QSO_SW8", "QSO_SW9", "QSO_SW10"])
        quasi_sequence_order_seq4_expected_values.loc[0] = [0.123287, 0.0, 0.079967, 0.039983, 0.04332, 0.143279, 0.013332, 0.046643, 0.103295, 0.089963]

        self.assertIsInstance(quasi_sequence_order_seq4, pd.DataFrame,
            'Expected output to be a DataFrame, got {}.'.format(type(quasi_sequence_order_seq4)))
        self.assertEqual(quasi_sequence_order_seq4.shape, (1, lag+20), 
            'Expected output to be of shape 1 x 25, got {}.'.format(quasi_sequence_order_seq4.shape)) 
        self.assertTrue(quasi_sequence_order_seq4.iloc[:, : 10].equals(quasi_sequence_order_seq4_expected_values), 
            "Output did not match expected values, got\n{}.".format(quasi_sequence_order_seq4.iloc[:, : 10]))
        for col in list(quasi_sequence_order_seq4.columns):
            #check all columns follow pattern of SOCNX or SOCNXY where x & y integers between 0 and 9
            self.assertTrue((bool(re.match(r'QSO_SW[0-9]', col)) or bool(re.match(r'QSO_SW[0-9][0-9]', col))), 
                "Column name doesn't match expected regex pattern: {}.".format(col))     
        self.assertTrue(quasi_sequence_order_seq4.any().isnull().sum()==0, 
            'Expected output to contain no null values.')        
        self.assertTrue(all(col == np.float64 for col in list(quasi_sequence_order_seq4.dtypes)), 
            "Descriptor dataframe dtypes expected to be np.float64, got {}.".format(list(quasi_sequence_order_seq4.dtypes)))
#5.)
        invalid_seq5 = "ABCDEF"
        invalid_seq6 = "OOOOO"
        invalid_lag = "BLAHBLAH"
        with (self.assertRaises(ValueError)):
            quasi_sequence_order_seq5 = protpy.quasi_sequence_order(invalid_seq5, lag=30, 
                distance_matrix=self.distance_matrix_schneider_wrede)
            quasi_sequence_order_seq6 = protpy.quasi_sequence_order(invalid_seq5, lag=10, 
                distance_matrix=self.distance_matrix_schneider_wrede)
            quasi_sequence_order_seq7 = protpy.quasi_sequence_order(self.protein_seq4, lag=invalid_lag, 
                distance_matrix=self.distance_matrix_schneider_wrede)
#6.)
        invalid_seq7 = 1234
        invalid_seq8 = False
        with (self.assertRaises(TypeError)):
            quasi_sequence_order_seq8 = protpy.quasi_sequence_order(invalid_seq7, lag=20, 
                distance_matrix=self.distance_matrix_schneider_wrede)
            quasi_sequence_order_seq9 = protpy.quasi_sequence_order(invalid_seq8, lag=20, 
                distance_matrix=self.distance_matrix_schneider_wrede)
            
    def test_quasi_sequence_order_all(self):
        """ Testing function that calculates QSO using both matrices for a sequence and lag. """
#1.)
        lag = 30
        qso_all_seq1 = protpy.quasi_sequence_order_all(self.protein_seq1, lag=lag)
        #testing first 10 columns 
        qso_all_seq1_expected_values = pd.DataFrame(columns=["QSO_SW1", "QSO_SW2", "QSO_SW3", "QSO_SW4", "QSO_SW5", "QSO_SW6", "QSO_SW7", "QSO_SW8", "QSO_SW9", "QSO_SW10"])
        qso_all_seq1_expected_values.loc[0] = [0.005692, 0.002643, 0.004947, 0.002846, 0.005625, 0.005354, 0.001016, 0.005286, 0.004066, 0.006708]

        self.assertIsInstance(qso_all_seq1, pd.DataFrame, 
            'Expected output to be a DataFrame, got {}.'.format(type(qso_all_seq1)))
        self.assertEqual(qso_all_seq1.shape, (1, (lag+20)*2),
            'Expected output to be of shape 1 x 100, got {}.'.format((1, (lag+20)*2), qso_all_seq1.shape)) 
        self.assertTrue(qso_all_seq1.iloc[:, : 10].equals(qso_all_seq1_expected_values), 
            "Output did not match expected values, got\n{}.".format(qso_all_seq1.iloc[:, : 10]))
        for col in list(qso_all_seq1.columns):
            #check all columns follow pattern of SOCNX or SOCNXY where x & y integers between 0 and 9
            self.assertTrue((bool(re.match(r'QSO_SW[0-9]', col)) or bool(re.match(r'QSO_SW[0-9][0-9]', col)) or
                bool(re.match(r'QSO_Grant[0-9]', col)) or bool(re.match(r'QSO_Grant[0-9][0-9]', col))), 
                "Column name doesn't match expected regex pattern: {}.".format(col))     
        self.assertTrue(qso_all_seq1.any().isnull().sum()==0,
            'Expected output to contain no null values.')        
        self.assertTrue(all(col == np.float64 for col in list(qso_all_seq1.dtypes)), 
            "Descriptor dataframe dtypes expected to be np.float64, got {}.".format(list(qso_all_seq1.dtypes)))
#2.)
        lag = 15
        qso_all_seq2 = protpy.quasi_sequence_order_all(self.protein_seq2, lag=lag)
        #testing first 10 columns 
        qso_all_seq2_expected_values = pd.DataFrame(columns=["QSO_SW1", "QSO_SW2", "QSO_SW3", "QSO_SW4", "QSO_SW5", "QSO_SW6", "QSO_SW7", "QSO_SW8", "QSO_SW9", "QSO_SW10"])
        qso_all_seq2_expected_values.loc[0] = [0.03651, 0.0, 0.023622, 0.015035, 0.013961, 0.048323, 0.00537, 0.011813, 0.03114, 0.027918]

        self.assertIsInstance(qso_all_seq2, pd.DataFrame, 
            'Expected output to be a DataFrame, got {}.'.format(type(qso_all_seq2)))
        self.assertEqual(qso_all_seq2.shape, (1, (lag+20)*2),
            'Expected output to be of shape 1 x 70, got {}.'.format((1, (lag+20)*2), qso_all_seq2.shape)) 
        self.assertTrue(qso_all_seq2.iloc[:, : 10].equals(qso_all_seq2_expected_values), 
            "Output did not match expected values, got\n{}.".format(qso_all_seq2.iloc[:, : 10]))
        for col in list(qso_all_seq2.columns):
            #check all columns follow pattern of SOCNX or SOCNXY where x & y integers between 0 and 9
            self.assertTrue((bool(re.match(r'QSO_SW[0-9]', col)) or bool(re.match(r'QSO_SW[0-9][0-9]', col)) or
                bool(re.match(r'QSO_Grant[0-9]', col)) or bool(re.match(r'QSO_Grant[0-9][0-9]', col))), 
                "Column name doesn't match expected regex pattern: {}.".format(col))     
        self.assertTrue(qso_all_seq2.any().isnull().sum()==0,
            'Expected output to contain no null values.')        
        self.assertTrue(all(col == np.float64 for col in list(qso_all_seq2.dtypes)), 
            "Descriptor dataframe dtypes expected to be np.float64, got {}.".format(list(qso_all_seq2.dtypes)))
#3.)
        lag = 10
        qso_all_seq3 = protpy.quasi_sequence_order_all(self.protein_seq3, lag=lag)
        #testing first 10 columns 
        qso_all_seq3_expected_values = pd.DataFrame(columns=["QSO_SW1", "QSO_SW2", "QSO_SW3", "QSO_SW4", "QSO_SW5", "QSO_SW6", "QSO_SW7", "QSO_SW8", "QSO_SW9", "QSO_SW10"])
        qso_all_seq3_expected_values.loc[0] = [0.310006, 0.23249, 0.077516, 0.23249, 0.310006, 0.155032, 0.0, 0.23249, 0.155032, 1.08505]

        self.assertIsInstance(qso_all_seq3, pd.DataFrame, 
            'Expected output to be a DataFrame, got {}.'.format(type(qso_all_seq3)))
        self.assertEqual(qso_all_seq3.shape, (1, (lag+20)*2),
            'Expected output to be of shape 1 x 60, got {}.'.format((1, (lag+20)*2), qso_all_seq3.shape)) 
        self.assertTrue(qso_all_seq3.iloc[:, : 10].equals(qso_all_seq3_expected_values), 
            "Output did not match expected values, got\n{}.".format(qso_all_seq3.iloc[:, : 10]))
        for col in list(qso_all_seq3.columns):
            #check all columns follow pattern of SOCNX or SOCNXY where x & y integers between 0 and 9
            self.assertTrue((bool(re.match(r'QSO_SW[0-9]', col)) or bool(re.match(r'QSO_SW[0-9][0-9]', col)) or
                bool(re.match(r'QSO_Grant[0-9]', col)) or bool(re.match(r'QSO_Grant[0-9][0-9]', col))), 
                "Column name doesn't match expected regex pattern: {}.".format(col))     
        self.assertTrue(qso_all_seq3.any().isnull().sum()==0,
            'Expected output to contain no null values.')        
        self.assertTrue(all(col == np.float64 for col in list(qso_all_seq3.dtypes)), 
            "Descriptor dataframe dtypes expected to be np.float64, got {}.".format(list(qso_all_seq3.dtypes)))
#4.)
        lag = 5
        qso_all_seq4 = protpy.quasi_sequence_order_all(self.protein_seq4, lag=lag)
        #testing first 10 columns 
        qso_all_seq4_expected_values = pd.DataFrame(columns=["QSO_SW1", "QSO_SW2", "QSO_SW3", "QSO_SW4", "QSO_SW5", "QSO_SW6", "QSO_SW7", "QSO_SW8", "QSO_SW9", "QSO_SW10"])
        qso_all_seq4_expected_values.loc[0] = [0.123287, 0.0, 0.079967, 0.039983, 0.04332, 0.143279, 0.013332, 0.046643, 0.103295, 0.089963]

        self.assertIsInstance(qso_all_seq4, pd.DataFrame, 
            'Expected output to be a DataFrame, got {}.'.format(type(qso_all_seq4)))
        self.assertEqual(qso_all_seq4.shape, (1, (lag+20)*2),
            'Expected output to be of shape 1 x 50, got {}.'.format((1, (lag+20)*2), qso_all_seq4.shape)) 
        self.assertTrue(qso_all_seq4.iloc[:, : 10].equals(qso_all_seq4_expected_values), 
            "Output did not match expected values, got\n{}.".format(qso_all_seq4.iloc[:, : 10]))
        for col in list(qso_all_seq4.columns):
            #check all columns follow pattern of SOCNX or SOCNXY where x & y integers between 0 and 9
            self.assertTrue((bool(re.match(r'QSO_SW[0-9]', col)) or bool(re.match(r'QSO_SW[0-9][0-9]', col)) or
                bool(re.match(r'QSO_Grant[0-9]', col)) or bool(re.match(r'QSO_Grant[0-9][0-9]', col))), 
                "Column name doesn't match expected regex pattern: {}.".format(col))     
        self.assertTrue(qso_all_seq4.any().isnull().sum()==0,
            'Expected output to contain no null values.')        
        self.assertTrue(all(col == np.float64 for col in list(qso_all_seq4.dtypes)), 
            "Descriptor dataframe dtypes expected to be np.float64, got {}.".format(list(qso_all_seq4.dtypes)))
#5.)
        invalid_seq5 = "ABCDEF"
        invalid_seq6 = "OOOOO"
        invalid_lag = "BLAHBLAH"
        with (self.assertRaises(ValueError)):
            quasi_sequence_order_seq5 = protpy.quasi_sequence_order_all(invalid_seq5, lag=30) 
            quasi_sequence_order_seq6 = protpy.quasi_sequence_order_all(invalid_seq6, lag=30) 
            quasi_sequence_order_seq7 = protpy.quasi_sequence_order_all(self.protein_seq3, lag=invalid_lag) 
#6.)
        invalid_seq7 = 1234
        invalid_seq8 = False
        with (self.assertRaises(TypeError)):
            quasi_sequence_order_seq8 = protpy.quasi_sequence_order_all(invalid_seq7, lag=30) 
            quasi_sequence_order_seq9 = protpy.quasi_sequence_order_all(invalid_seq8, lag=30) 