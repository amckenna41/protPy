###########################################################################
#############      protPy - Autocorrelation Module Tests      #############
###########################################################################

import pandas as pd
import numpy as np
import os
import unittest
import re
from Bio import SeqIO
unittest.TestLoader.sortTestMethodsUsing = None

import protpy as protpy

class ProtpyAutocorrelationTests(unittest.TestCase):
    """
    Test suite for testing autcorrelation module and functionality 
    in protpy package, including the MoreauBroto, Moran and Geary
    Autocorrelation descriptors.

    Test Cases
    ----------
    test_moreaubroto_autocorrelation:
        testing correct protpy MoreauBroto Autocorrelation functionality.
    test_moran_autocorrelation:
        testing correct protpy Moran Autocorrelation functionality.
    test_geary_autocorrelation:
        testing correct protpy Geary Autocorrelation functionality.
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
        self.amino_acids = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", 
            "Q", "R", "S", "T", "V", "W", "Y"]

    def test_moreaubroto_autocorrelation(self):
        """ Testing moreaubroto autocorrelation descriptor attributes and functionality. """
        moreaubroto_lag = list(range(5, 35, 5))
        properties = ["CIDH920105", "BHAR880101", "CHAM820101", "CHAM820102", 
            "CHOC760101", "BIGC670101", "CHAM810101", "DAYM780201"]
#1.)
        #testing using range of lag values
        for lag in moreaubroto_lag:
            moreaubroto_seq1 = protpy.moreaubroto_autocorrelation(self.protein_seq1, lag=lag, 
                properties=properties, normalize=True)
            self.assertIsInstance(moreaubroto_seq1, pd.DataFrame, 
                'Expected output to be of type DataFrame, got {}.'.format(type(moreaubroto_seq1)))
            self.assertEqual(moreaubroto_seq1.shape, (1, lag*len(properties)), 
                'Expected output to be of shape {}, got {}.'.format((1, lag*len(properties)), moreaubroto_seq1.shape)) 
            self.assertTrue(moreaubroto_seq1.any().isnull().sum()==0, 
                'Expected output to contain no null values.')        
            self.assertTrue(all(col == np.float64 for col in list(moreaubroto_seq1.dtypes)), 
                "Expected output values to be of datatype np.float64, got {}.".format(list(moreaubroto_seq1.dtypes)))

            #check all columns follow pattern of MBAutoX_Y where x is the asscession number of
            #the AAindex record and y is the count of the descriptor
            for col in list(moreaubroto_seq1.columns):
                self.assertTrue(bool(re.match(r"MBAuto_[A-Z0-9]{10}_[0-9]", col)), 
                    "Column name doesn't match expected regex pattern: {}.".format(col))     
#2.)
        #testing using range of lag values
        for lag in moreaubroto_lag:
            moreaubroto_seq2 = protpy.moreaubroto_autocorrelation(self.protein_seq2, lag=lag, 
                properties=properties, normalize=True)
            self.assertIsInstance(moreaubroto_seq2, pd.DataFrame,
                'Expected output to be of type DataFrame, got {}.'.format(type(moreaubroto_seq2)))
            self.assertEqual(moreaubroto_seq2.shape, (1, lag*len(properties)),
                'Expected output to be of shape {}, got {}.'.format((1, lag*len(properties)), moreaubroto_seq2.shape)) 
            self.assertTrue(moreaubroto_seq2.any().isnull().sum()==0, 
                'Expected output to contain no null values.')        
            self.assertTrue(all(col == np.float64 for col in list(moreaubroto_seq2.dtypes)), 
                "Expected output values to be of datatype np.float64, got {}.".format(list(moreaubroto_seq2.dtypes)))

            #check all columns follow pattern of MBAutoX_Y where x is the asscession number of
            #the AAindex record and y is the count of the descriptor
            for col in list(moreaubroto_seq2.columns):
                self.assertTrue(bool(re.match(r"MBAuto_[A-Z0-9]{10}_[0-9]", col)), 
                    "Column name doesn't match expected regex pattern: {}.".format(col))
#3.)
        #testing using range of lag values
        for lag in moreaubroto_lag:
            moreaubroto_seq3 = protpy.moreaubroto_autocorrelation(self.protein_seq3, lag=lag, 
                properties=properties, normalize=True)
            self.assertIsInstance(moreaubroto_seq3, pd.DataFrame, 
                'Expected output to be of type DataFrame, got {}.'.format(type(moreaubroto_seq3)))
            self.assertEqual(moreaubroto_seq3.shape, (1, lag*len(properties)), 
                'Expected output to be of shape {}, got {}.'.format((1, lag*len(properties)), moreaubroto_seq3.shape)) 
            self.assertTrue(moreaubroto_seq3.any().isnull().sum()==0, 
                'Expected output to contain no null values.')        
            self.assertTrue(all(col == np.float64 for col in list(moreaubroto_seq3.dtypes)), 
                "Expected output values to be of datatype np.float64, got {}.".format(list(moreaubroto_seq3.dtypes)))

            #check all columns follow pattern of MBAutoX_Y where x is the asscession number of
            # the AAindex record and y is the count of the descriptor
            for col in list(moreaubroto_seq3.columns):
                self.assertTrue(bool(re.match(r"MBAuto_[A-Z0-9]{10}_[0-9]", col)), 
                    "Column name doesn't match expected regex pattern: {}.".format(col))
#4.)
        #testing using range of lag values
        for lag in moreaubroto_lag:
            moreaubroto_seq4 = protpy.moreaubroto_autocorrelation(self.protein_seq4, lag=lag, 
                properties=properties, normalize=True)
            self.assertIsInstance(moreaubroto_seq4, pd.DataFrame, 
                'Expected output to be of type DataFrame, got {}.'.format(type(moreaubroto_seq4)))
            self.assertEqual(moreaubroto_seq4.shape, (1, lag*len(properties)), 
                'Expected output to be of shape {}, got {}.'.format((1, lag*len(properties)), moreaubroto_seq4.shape)) 
            self.assertTrue(moreaubroto_seq4.any().isnull().sum()==0,
                'Expected output to contain no null values.')        
            self.assertTrue(all(col == np.float64 for col in list(moreaubroto_seq4.dtypes)), 
                "Expected output values to be of datatype np.float64, got {}.".format(list(moreaubroto_seq4.dtypes)))

            #check all columns follow pattern of MBAutoX_Y where x is the asscession number of
            #the AAindex record and y is the count of the descriptor
            for col in list(moreaubroto_seq4.columns):
                self.assertTrue(bool(re.match(r"MBAuto_[A-Z0-9]{10}_[0-9]", col)), 
                    "Column name doesn't match expected regex pattern: {}.".format(col))
#5.)
        invalid_seq5 = "ABCDEF"
        invalid_seq6 = "OOOOO"
        with (self.assertRaises(ValueError)):
            moreaubroto_seq5 = protpy.moreaubroto_autocorrelation(invalid_seq5)
            moreaubroto_seq6 = protpy.moreaubroto_autocorrelation(invalid_seq6)
#6.)
        invalid_seq7 = 12345
        invalid_seq8 = False
        with (self.assertRaises(TypeError)):
            moreaubroto_seq7 = protpy.moreaubroto_autocorrelation(invalid_seq7)
            moreaubroto_seq8 = protpy.moreaubroto_autocorrelation(invalid_seq8)

    def test_moran_autocorrelation(self):
        """ Testing moran autocorrelation descriptor attributes and functionality. """
        moran_auto_lag = list(range(5, 35, 5))
        properties = ["CIDH920105", "BHAR880101", "CHAM820101", "CHAM820102", 
            "CHOC760101", "BIGC670101", "CHAM810101", "DAYM780201"]
#1.)
        #testing using range of lag values
        for lag in moran_auto_lag:
            moran_auto_seq1 = protpy.moran_autocorrelation(self.protein_seq1, lag=lag, 
                properties=properties, normalize=True)
            self.assertIsInstance(moran_auto_seq1, pd.DataFrame, 
                'Expected output to be of type DataFrame, got {}.'.format(type(moran_auto_seq1)))
            self.assertEqual(moran_auto_seq1.shape, (1, lag*len(properties)),
                'Expected output to be of shape {}, got {}.'.format((1, lag*len(properties)), moran_auto_seq1.shape)) 
            self.assertTrue(moran_auto_seq1.any().isnull().sum()==0,
                'Expected output to contain no null values.')        
            self.assertTrue(all(col == np.float64 for col in list(moran_auto_seq1.dtypes)), 
                "Expected output values to be of datatype np.float64, got {}.".format(list(moran_auto_seq1.dtypes)))

            #check all columns follow pattern of MAuto_X_Y where x is the asscession number of
            #the AAindex record and y is the count of the descriptor
            for col in list(moran_auto_seq1.columns):
                self.assertTrue(bool(re.match(r"MAuto_[A-Z0-9]{10}_[0-9]", col)), 
                    "Column name doesn't match expected regex pattern: {}.".format(col))
#2.)
        #testing using range of lag values
        for lag in moran_auto_lag:
            moran_auto_seq2 = protpy.moran_autocorrelation(self.protein_seq2, lag=lag, 
                properties=properties, normalize=True)
            self.assertIsInstance(moran_auto_seq2, pd.DataFrame, 
                'Expected output to be of type DataFrame, got {}.'.format(type(moran_auto_seq2)))
            self.assertEqual(moran_auto_seq2.shape, (1, lag*len(properties)), 
                'Expected output to be of shape {}, got {}.'.format((1, lag*len(properties)), moran_auto_seq2.shape)) 
            self.assertTrue(moran_auto_seq2.any().isnull().sum()==0,
                'Expected output to contain no null values.')        
            self.assertTrue(all(col == np.float64 for col in list(moran_auto_seq2.dtypes)), 
                "Expected output values to be of datatype np.float64, got {}.".format(list(moran_auto_seq2.dtypes)))

            #check all columns follow pattern of MAuto_X_Y where x is the asscession number of
            #   the AAindex record and y is the count of the descriptor
            for col in list(moran_auto_seq2.columns):
                self.assertTrue(bool(re.match(r"MAuto_[A-Z0-9]{10}_[0-9]", col)), 
                    "Column name doesn't match expected regex pattern: {}.".format(col))
#3.)
        #testing using range of lag values
        for lag in moran_auto_lag:
            moran_auto_seq3 = protpy.moran_autocorrelation(self.protein_seq3, lag=lag, 
                properties=properties, normalize=True)
            self.assertIsInstance(moran_auto_seq3, pd.DataFrame, 
                'Expected output to be of type DataFrame, got {}.'.format(type(moran_auto_seq3)))
            self.assertEqual(moran_auto_seq3.shape, (1, lag*len(properties)), 
                'Expected output to be of shape {}, got {}.'.format((1, lag*len(properties)), moran_auto_seq3.shape)) 
            self.assertTrue(moran_auto_seq3.any().isnull().sum()==0,
                'Expected output to contain no null values.')        
            self.assertTrue(all(col == np.float64 for col in list(moran_auto_seq3.dtypes)), 
                "Expected output values to be of datatype np.float64, got {}.".format(list(moran_auto_seq3.dtypes)))

            #check all columns follow pattern of MAuto_X_Y where x is the asscession number of
            #the AAindex record and y is the count of the descriptor
            for col in list(moran_auto_seq3.columns):
                self.assertTrue(bool(re.match(r"MAuto_[A-Z0-9]{10}_[0-9]", col)), 
                    "Column name doesn't match expected regex pattern: {}".format(col))
#4.)
        #testing using range of lag values
        for lag in moran_auto_lag:
            moran_auto_seq4 = protpy.moran_autocorrelation(self.protein_seq4, lag=lag, 
                properties=properties, normalize=True)
            self.assertIsInstance(moran_auto_seq4, pd.DataFrame, 
                'Expected output to be of type DataFrame, got {}.'.format(type(moran_auto_seq4)))
            self.assertEqual(moran_auto_seq4.shape, (1, lag*len(properties)), 
                'Expected output to be of shape {}, got {}.'.format((1, lag*len(properties)), moran_auto_seq4.shape)) 
            self.assertTrue(moran_auto_seq4.any().isnull().sum()==0,
                'Expected output to contain no null values.')        
            self.assertTrue(all(col == np.float64 for col in list(moran_auto_seq4.dtypes)), 
                "Expected output values to be of datatype np.float64, got {}.".format(list(moran_auto_seq4.dtypes)))

            #check all columns follow pattern of MAuto_X_Y where x is the asscession number of
            #the AAindex record and y is the count of the descriptor
            for col in list(moran_auto_seq4.columns):
                self.assertTrue(bool(re.match(r"MAuto_[A-Z0-9]{10}_[0-9]", col)), 
                    "Column name doesn't match expected regex pattern: {}".format(col))
#5.)
        invalid_seq5 = "ABCDEF"
        invalid_seq6 = "OOOOO"
        with (self.assertRaises(ValueError)):
            moran_auto_seq5 = protpy.moran_autocorrelation(invalid_seq5)
            moran_auto_seq6 = protpy.moran_autocorrelation(invalid_seq6)
#6.)
        invalid_seq7 = 12345
        invalid_seq8 = False
        with (self.assertRaises(TypeError)):
            moran_auto_seq7 = protpy.moran_autocorrelation(invalid_seq7)
            moran_auto_seq8 = protpy.moran_autocorrelation(invalid_seq8)

    def test_geary_autocorrelation(self):
        """ Testing geary autocorrelation descriptor attributes and functionality. """
        geary_auto_lag = list(range(5, 35, 5))
        properties = ["CIDH920105", "BHAR880101", "CHAM820101", "CHAM820102", 
            "CHOC760101", "BIGC670101", "CHAM810101", "DAYM780201"]
#1.)
        #testing using range of lag values
        for lag in geary_auto_lag:
            geary_auto_seq1 = protpy.geary_autocorrelation(self.protein_seq1, lag=lag, 
                properties=properties, normalize=True)
            self.assertIsInstance(geary_auto_seq1, pd.DataFrame, 
                'Expected output to be of type DataFrame, got {}.'.format(type(geary_auto_seq1)))
            self.assertEqual(geary_auto_seq1.shape, (1, lag*len(properties)), 
                'Expected output to be of shape {}, got {}.'.format((1, lag*len(properties)), geary_auto_seq1.shape)) 
            self.assertTrue(geary_auto_seq1.any().isnull().sum()==0,
                'Expected output to contain no null values.')        
            self.assertTrue(all(col == np.float64 for col in list(geary_auto_seq1.dtypes)), 
                "Expected output values to be of datatype np.float64, got {}.".format(list(geary_auto_seq1.dtypes)))

            #check all columns follow pattern of GAuto_X_Y where x is the asscession number of
            #the AAindex record and y is the count of the descriptor
            for col in list(geary_auto_seq1.columns):
                self.assertTrue(bool(re.match(r"GAuto_[A-Z0-9]{10}_[0-9]", col)), 
                    "Column name doesn't match expected regex pattern: {}.".format(col))
#2.)
        #testing using range of lag values
        for lag in geary_auto_lag:
            geary_auto_seq2 = protpy.geary_autocorrelation(self.protein_seq2, lag=lag, 
                properties=properties, normalize=True)
            self.assertIsInstance(geary_auto_seq2, pd.DataFrame, 
                'Expected output to be of type DataFrame, got {}.'.format(type(geary_auto_seq2)))
            self.assertEqual(geary_auto_seq2.shape, (1, lag*len(properties)), 
                'Expected output to be of shape {}, got {}.'.format((1, lag*len(properties)), geary_auto_seq2.shape)) 
            self.assertTrue(geary_auto_seq2.any().isnull().sum()==0,
                'Expected output to contain no null values.')        
            self.assertTrue(all(col == np.float64 for col in list(geary_auto_seq2.dtypes)), 
                "Expected output values to be of datatype np.float64, got {}.".format(list(geary_auto_seq2.dtypes)))

            #check all columns follow pattern of GAuto_X_Y where x is the asscession number of
            #the AAindex record and y is the count of the descriptor
            for col in list(geary_auto_seq2.columns):
                self.assertTrue(bool(re.match(r"GAuto_[A-Z0-9]{10}_[0-9]", col)), 
                    "Column name doesn't match expected regex pattern: {}.".format(col))
#3.)
        #testing using range of lag values
        for lag in geary_auto_lag:
            geary_auto_seq3 = protpy.geary_autocorrelation(self.protein_seq3, lag=lag, 
                properties=properties, normalize=True)
            self.assertIsInstance(geary_auto_seq3, pd.DataFrame, 
                'Expected output to be of type DataFrame, got {}.'.format(type(geary_auto_seq3)))
            self.assertEqual(geary_auto_seq3.shape, (1, lag*len(properties)), 
                'Expected output to be of shape {}, got {}.'.format((1, lag*len(properties)), geary_auto_seq3.shape)) 
            self.assertTrue(geary_auto_seq3.any().isnull().sum()==0,
                'Expected output to contain no null values.')        
            self.assertTrue(all(col == np.float64 for col in list(geary_auto_seq3.dtypes)), 
                "Expected output values to be of datatype np.float64, got {}.".format(list(geary_auto_seq3.dtypes)))

            #check all columns follow pattern of GAuto_X_Y where x is the asscession number of
            #the AAindex record and y is the count of the descriptor
            for col in list(geary_auto_seq3.columns):
                self.assertTrue(bool(re.match(r"GAuto_[A-Z0-9]{10}_[0-9]", col)), 
                    "Column name doesn't match expected regex pattern: {}".format(col))
#4.)
        #testing using range of lag values
        for lag in geary_auto_lag:
            geary_auto_seq4 = protpy.geary_autocorrelation(self.protein_seq4, lag=lag, 
                properties=properties, normalize=True)
            self.assertIsInstance(geary_auto_seq4, pd.DataFrame, 
                'Expected output to be of type DataFrame, got {}.'.format(type(geary_auto_seq4)))
            self.assertEqual(geary_auto_seq4.shape, (1, lag*len(properties)), 
                'Expected output to be of shape {}, got {}.'.format((1, lag*len(properties)), geary_auto_seq4.shape)) 
            self.assertTrue(geary_auto_seq4.any().isnull().sum()==0,
                'Expected output to contain no null values.')        
            self.assertTrue(all(col == np.float64 for col in list(geary_auto_seq4.dtypes)), 
                "Expected output values to be of datatype np.float64, got {}.".format(list(geary_auto_seq4.dtypes)))

            #check all columns follow pattern of GAuto_X_Y where x is the asscession number of
            #the AAindex record and y is the count of the descriptor
            for col in list(geary_auto_seq4.columns):
                self.assertTrue(bool(re.match(r"GAuto_[A-Z0-9]{10}_[0-9]", col)), 
                    "Column name doesn't match expected regex pattern: {}".format(col))
#5.)
        invalid_seq5 = "ABCDEF"
        invalid_seq6 = "OOOOO"
        with (self.assertRaises(ValueError)):
            geary_auto_seq5 = protpy.geary_autocorrelation(invalid_seq5)
            geary_auto_seq6 = protpy.geary_autocorrelation(invalid_seq6)
#6.)
        invalid_seq7 = 12345
        invalid_seq8 = True
        with (self.assertRaises(TypeError)):
            geary_auto_seq7 = protpy.geary_autocorrelation(invalid_seq7)
            geary_auto_seq8 = protpy.geary_autocorrelation(invalid_seq8)