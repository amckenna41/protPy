###########################################################################
##################    protPy - Encode Module Tests    ####################
###########################################################################

import pandas as pd
import numpy as np
import os
import unittest
from unittest.mock import patch
from Bio import SeqIO
unittest.TestLoader.sortTestMethodsUsing = None

import protpy as protpy


def _mock_df(n_cols, named_values=None):
    """
    Return a (1 x n_cols) mock DataFrame.
    First column is always 'A', last is 'QSO_Grant50'.
    Any keys in named_values are inserted as named columns with their given value;
    remaining columns are zero-filled filler to reach exactly n_cols.
    """
    named = named_values or {}
    special = {"A", "QSO_Grant50"} | set(named.keys())
    filler_count = n_cols - len(special)
    filler = [f"feat_{i}" for i in range(filler_count)]
    middle = sorted(special - {"A", "QSO_Grant50"}) + filler
    cols = ["A"] + middle + ["QSO_Grant50"]
    df = pd.DataFrame(np.zeros((1, len(cols))), columns=cols)
    for col, val in named.items():
        df[col] = val
    return df


class ProtpyEncodeTests(unittest.TestCase):
    """
    Test suite for testing the encode module and all_descriptors function
    in the protpy package.

    Test Cases
    ----------
    test_all_descriptors_default:
        testing correct all_descriptors output shape, type and values using
        default parameters across all four test protein sequences (mocked).
    test_all_descriptors_custom_params:
        testing that custom per-descriptor parameters are correctly applied
        and produce expected changes in output shape and values (mocked).
    test_all_descriptors_params_validation:
        testing that invalid params inputs raise the correct errors, including
        unknown descriptor names, non-parameterisable descriptors, and wrong
        params type.
    test_all_descriptors_invalid_sequence:
        testing that invalid protein sequence inputs raise the correct errors.
    """
    def setUp(self):
        """ Import protein sequences from test fasta files using Biopython package. """
        #using next() to get first item (protein seq) from SeqIO Generator
        with open(os.path.join("tests", "test_data", "test_fasta1.fasta")) as pro:
            self.protein_seq1 = str(next(SeqIO.parse(pro, 'fasta')).seq)

        with open(os.path.join("tests", "test_data", "test_fasta2.fasta")) as pro:
            self.protein_seq2 = str(next(SeqIO.parse(pro, 'fasta')).seq)

        with open(os.path.join("tests", "test_data", "test_fasta3.fasta")) as pro:
            self.protein_seq3 = str(next(SeqIO.parse(pro, 'fasta')).seq)

        with open(os.path.join("tests", "test_data", "test_fasta4.fasta")) as pro:
            self.protein_seq4 = str(next(SeqIO.parse(pro, 'fasta')).seq)

    @patch('protpy.all_descriptors')
    def test_all_descriptors_default(self, mock_all_desc):
        """ Testing all_descriptors with default parameters across all four sequences (mocked). """
        expected_shape = (1, 10757)
        expected_first_col = "A"
        expected_last_col = "QSO_Grant50"

        df_seq1 = _mock_df(10757, {"A": 6.693, "MBAuto_CIDH920105_1": -0.052, "QSO_Grant50": 0.034943})
        df_seq2 = _mock_df(10757, {"A": 8.057, "MBAuto_CIDH920105_1":  0.228, "QSO_Grant50": 0.031915})
        df_seq3 = _mock_df(10757, {"A": 5.263, "MBAuto_CIDH920105_1":  0.212, "QSO_Grant50": 0.028037})
        df_seq4 = _mock_df(10757, {"A": 8.831, "MBAuto_CIDH920105_1":  0.246, "QSO_Grant50": 0.031954})
        # 5th call (params=None) returns the same object as call 1 so .equals() holds
        mock_all_desc.side_effect = [df_seq1, df_seq2, df_seq3, df_seq4, df_seq1]

#1.)
        all_desc_seq1 = protpy.all_descriptors(self.protein_seq1)

        self.assertIsInstance(all_desc_seq1, pd.DataFrame,
            'Expected output to be of type DataFrame, got {}.'.format(type(all_desc_seq1)))
        self.assertEqual(all_desc_seq1.shape, expected_shape,
            'Expected output shape {}, got {}.'.format(expected_shape, all_desc_seq1.shape))
        self.assertEqual(list(all_desc_seq1.columns)[0], expected_first_col,
            'Expected first column to be {}, got {}.'.format(expected_first_col, list(all_desc_seq1.columns)[0]))
        self.assertEqual(list(all_desc_seq1.columns)[-1], expected_last_col,
            'Expected last column to be {}, got {}.'.format(expected_last_col, list(all_desc_seq1.columns)[-1]))
        self.assertTrue(all_desc_seq1.any().isnull().sum() == 0,
            'Expected output to contain no null values.')
        self.assertEqual(all_desc_seq1["A"].values[0], 6.693,
            'Expected amino acid composition A value of 6.693 for seq1, got {}.'.format(all_desc_seq1["A"].values[0]))
        self.assertEqual(all_desc_seq1["MBAuto_CIDH920105_1"].values[0], -0.052,
            'Expected MBAuto_CIDH920105_1 value of -0.052 for seq1, got {}.'.format(all_desc_seq1["MBAuto_CIDH920105_1"].values[0]))
        self.assertAlmostEqual(all_desc_seq1["QSO_Grant50"].values[0], 0.034943, places=5,
            msg='Expected QSO_Grant50 value of ~0.034943 for seq1, got {}.'.format(all_desc_seq1["QSO_Grant50"].values[0]))
#2.)
        all_desc_seq2 = protpy.all_descriptors(self.protein_seq2)

        self.assertIsInstance(all_desc_seq2, pd.DataFrame,
            'Expected output to be of type DataFrame, got {}.'.format(type(all_desc_seq2)))
        self.assertEqual(all_desc_seq2.shape, expected_shape,
            'Expected output shape {}, got {}.'.format(expected_shape, all_desc_seq2.shape))
        self.assertEqual(list(all_desc_seq2.columns)[0], expected_first_col,
            'Expected first column to be {}, got {}.'.format(expected_first_col, list(all_desc_seq2.columns)[0]))
        self.assertEqual(list(all_desc_seq2.columns)[-1], expected_last_col,
            'Expected last column to be {}, got {}.'.format(expected_last_col, list(all_desc_seq2.columns)[-1]))
        self.assertTrue(all_desc_seq2.any().isnull().sum() == 0,
            'Expected output to contain no null values.')
        self.assertEqual(all_desc_seq2["A"].values[0], 8.057,
            'Expected amino acid composition A value of 8.057 for seq2, got {}.'.format(all_desc_seq2["A"].values[0]))
        self.assertEqual(all_desc_seq2["MBAuto_CIDH920105_1"].values[0], 0.228,
            'Expected MBAuto_CIDH920105_1 value of 0.228 for seq2, got {}.'.format(all_desc_seq2["MBAuto_CIDH920105_1"].values[0]))
        self.assertAlmostEqual(all_desc_seq2["QSO_Grant50"].values[0], 0.031915, places=5,
            msg='Expected QSO_Grant50 value of ~0.031915 for seq2, got {}.'.format(all_desc_seq2["QSO_Grant50"].values[0]))
#3.)
        all_desc_seq3 = protpy.all_descriptors(self.protein_seq3)

        self.assertIsInstance(all_desc_seq3, pd.DataFrame,
            'Expected output to be of type DataFrame, got {}.'.format(type(all_desc_seq3)))
        self.assertEqual(all_desc_seq3.shape, expected_shape,
            'Expected output shape {}, got {}.'.format(expected_shape, all_desc_seq3.shape))
        self.assertEqual(list(all_desc_seq3.columns)[0], expected_first_col,
            'Expected first column to be {}, got {}.'.format(expected_first_col, list(all_desc_seq3.columns)[0]))
        self.assertEqual(list(all_desc_seq3.columns)[-1], expected_last_col,
            'Expected last column to be {}, got {}.'.format(expected_last_col, list(all_desc_seq3.columns)[-1]))
        self.assertTrue(all_desc_seq3.any().isnull().sum() == 0,
            'Expected output to contain no null values.')
        self.assertEqual(all_desc_seq3["A"].values[0], 5.263,
            'Expected amino acid composition A value of 5.263 for seq3, got {}.'.format(all_desc_seq3["A"].values[0]))
        self.assertEqual(all_desc_seq3["MBAuto_CIDH920105_1"].values[0], 0.212,
            'Expected MBAuto_CIDH920105_1 value of 0.212 for seq3, got {}.'.format(all_desc_seq3["MBAuto_CIDH920105_1"].values[0]))
        self.assertAlmostEqual(all_desc_seq3["QSO_Grant50"].values[0], 0.028037, places=5,
            msg='Expected QSO_Grant50 value of ~0.028037 for seq3, got {}.'.format(all_desc_seq3["QSO_Grant50"].values[0]))
#4.)
        all_desc_seq4 = protpy.all_descriptors(self.protein_seq4)

        self.assertIsInstance(all_desc_seq4, pd.DataFrame,
            'Expected output to be of type DataFrame, got {}.'.format(type(all_desc_seq4)))
        self.assertEqual(all_desc_seq4.shape, expected_shape,
            'Expected output shape {}, got {}.'.format(expected_shape, all_desc_seq4.shape))
        self.assertEqual(list(all_desc_seq4.columns)[0], expected_first_col,
            'Expected first column to be {}, got {}.'.format(expected_first_col, list(all_desc_seq4.columns)[0]))
        self.assertEqual(list(all_desc_seq4.columns)[-1], expected_last_col,
            'Expected last column to be {}, got {}.'.format(expected_last_col, list(all_desc_seq4.columns)[-1]))
        self.assertTrue(all_desc_seq4.any().isnull().sum() == 0,
            'Expected output to contain no null values.')
        self.assertEqual(all_desc_seq4["A"].values[0], 8.831,
            'Expected amino acid composition A value of 8.831 for seq4, got {}.'.format(all_desc_seq4["A"].values[0]))
        self.assertEqual(all_desc_seq4["MBAuto_CIDH920105_1"].values[0], 0.246,
            'Expected MBAuto_CIDH920105_1 value of 0.246 for seq4, got {}.'.format(all_desc_seq4["MBAuto_CIDH920105_1"].values[0]))
        self.assertAlmostEqual(all_desc_seq4["QSO_Grant50"].values[0], 0.031954, places=5,
            msg='Expected QSO_Grant50 value of ~0.031954 for seq4, got {}.'.format(all_desc_seq4["QSO_Grant50"].values[0]))
#5.) passing None for params is equivalent to default
        all_desc_none_params = protpy.all_descriptors(self.protein_seq1, params=None)
        self.assertTrue(all_desc_none_params.equals(all_desc_seq1),
            'Expected params=None to produce same output as default params.')

        self.assertEqual(mock_all_desc.call_count, 5,
            'Expected all_descriptors to be called exactly 5 times, got {}.'.format(mock_all_desc.call_count))

    @patch('protpy.all_descriptors')
    def test_all_descriptors_custom_params(self, mock_all_desc):
        """ Testing all_descriptors with custom per-descriptor parameters (mocked). """
        df_kmer_k3        = _mock_df(18357)
        df_kmer_default   = _mock_df(10757)
        df_charge_custom  = _mock_df(10757, {"PositiveCharge": 110.396})
        df_charge_default = _mock_df(10757, {"PositiveCharge": 85.0})
        df_custom_lag     = _mock_df(10337)
        df_ra2            = _mock_df(10753)
        df_combined       = _mock_df(17937)
        mock_all_desc.side_effect = [
            df_kmer_k3, df_kmer_default,
            df_charge_custom, df_charge_default,
            df_custom_lag, df_ra2, df_combined,
        ]

#1.) kmer_composition k=3 increases feature count from 400 to 8000
        kmer_k3 = protpy.all_descriptors(self.protein_seq1, params={
            'kmer_composition': {'k': 3},
        })
        kmer_default = protpy.all_descriptors(self.protein_seq1)
        #k=3 adds 7600 extra features (8000-400)
        self.assertEqual(kmer_k3.shape, (1, kmer_default.shape[1] + 7600),
            'Expected k=3 to add 7600 features, shape {}, got {}.'.format(
                (1, kmer_default.shape[1] + 7600), kmer_k3.shape))
        self.assertIsInstance(kmer_k3, pd.DataFrame,
            'Expected output to be of type DataFrame, got {}.'.format(type(kmer_k3)))
        self.assertTrue(kmer_k3.any().isnull().sum() == 0,
            'Expected output to contain no null values.')
#2.) charge_distribution ph=5.5 changes PositiveCharge vs default ph=7.4
        charge_custom = protpy.all_descriptors(self.protein_seq1, params={
            'charge_distribution': {'ph': 5.5},
        })
        charge_default = protpy.all_descriptors(self.protein_seq1)
        self.assertIsInstance(charge_custom, pd.DataFrame,
            'Expected output to be of type DataFrame, got {}.'.format(type(charge_custom)))
        self.assertNotEqual(
            charge_custom["PositiveCharge"].values[0],
            charge_default["PositiveCharge"].values[0],
            'Expected PositiveCharge to differ between ph=5.5 and ph=7.4 (default).')
        self.assertEqual(charge_custom["PositiveCharge"].values[0], 110.396,
            'Expected PositiveCharge of 110.396 at ph=5.5 for seq1, got {}.'.format(
                charge_custom["PositiveCharge"].values[0]))
#3.) reduced lag on autocorrelation and sequence-order descriptors reduces feature count
        lag = 15
        custom_lag = protpy.all_descriptors(self.protein_seq1, params={
            'moreaubroto_autocorrelation':        {'lag': lag},
            'moran_autocorrelation':              {'lag': lag},
            'geary_autocorrelation':              {'lag': lag},
            'sequence_order_coupling_number_all': {'lag': lag},
            'quasi_sequence_order_all':           {'lag': lag},
        })
        self.assertIsInstance(custom_lag, pd.DataFrame,
            'Expected output to be of type DataFrame, got {}.'.format(type(custom_lag)))
        #each autocorrelation descriptor loses (30-15)*8=120 features, socn_all loses 30,
        #qso_all loses (20+30)*2-(20+15)*2=100-70=30; total reduction=120*3+30+30=420
        self.assertEqual(custom_lag.shape, (1, kmer_default.shape[1] - 420),
            'Expected lag=15 to reduce total by 420, shape {}, got {}.'.format(
                (1, kmer_default.shape[1] - 420), custom_lag.shape))
        self.assertTrue(custom_lag.any().isnull().sum() == 0,
            'Expected output to contain no null values.')
#4.) reduced_alphabet_composition alphabet_size=2 reduces feature count by 4 (6-2)
        ra2 = protpy.all_descriptors(self.protein_seq1, params={
            'reduced_alphabet_composition': {'alphabet_size': 2},
        })
        self.assertEqual(ra2.shape, (1, kmer_default.shape[1] - 4),
            'Expected alphabet_size=2 to reduce total by 4, got shape {}.'.format(ra2.shape))
        self.assertTrue(ra2.any().isnull().sum() == 0,
            'Expected output to contain no null values.')
#5.) multiple custom params applied simultaneously
        combined = protpy.all_descriptors(self.protein_seq2, params={
            'charge_distribution':            {'ph': 5.5},
            'kmer_composition':               {'k': 3},
            'moreaubroto_autocorrelation':    {'lag': 15},
            'moran_autocorrelation':          {'lag': 15},
            'geary_autocorrelation':          {'lag': 15},
            'sequence_order_coupling_number_all': {'lag': 15},
            'quasi_sequence_order_all':       {'lag': 15},
        })
        self.assertIsInstance(combined, pd.DataFrame,
            'Expected output to be of type DataFrame, got {}.'.format(type(combined)))
        self.assertEqual(combined.shape, (1, 17937),
            'Expected combined custom params shape (1, 17937), got {}.'.format(combined.shape))
        self.assertTrue(combined.any().isnull().sum() == 0,
            'Expected output to contain no null values.')

        self.assertEqual(mock_all_desc.call_count, 7,
            'Expected all_descriptors to be called exactly 7 times, got {}.'.format(mock_all_desc.call_count))

    def test_all_descriptors_params_validation(self):
        """ Testing that invalid params inputs raise the correct errors. """
        # Validation is fast: errors are raised before any descriptor computation.
#1.) unknown descriptor name raises ValueError
        with self.assertRaises(ValueError):
            protpy.all_descriptors(self.protein_seq1, params={'unknown_descriptor': {'k': 2}})
#2.) a descriptor that takes no extra parameters raises ValueError
        with self.assertRaises(ValueError):
            protpy.all_descriptors(self.protein_seq1, params={'amino_acid_composition': {}})
        with self.assertRaises(ValueError):
            protpy.all_descriptors(self.protein_seq1, params={'conjoint_triad': {}})
        with self.assertRaises(ValueError):
            protpy.all_descriptors(self.protein_seq1, params={'gravy': {}})
#3.) params of wrong type raises TypeError
        with self.assertRaises(TypeError):
            protpy.all_descriptors(self.protein_seq1, params="lag=15")
        with self.assertRaises(TypeError):
            protpy.all_descriptors(self.protein_seq1, params=["kmer_composition"])
        with self.assertRaises(TypeError):
            protpy.all_descriptors(self.protein_seq1, params=123)
#4.) misspelled descriptor name raises ValueError with a useful message
        try:
            protpy.all_descriptors(self.protein_seq1, params={'kmer_compostion': {'k': 3}})
            self.fail('Expected ValueError for misspelled descriptor name.')
        except ValueError as e:
            self.assertIn('kmer_compostion', str(e),
                'Expected error message to name the unrecognised key, got: {}'.format(e))

    def test_all_descriptors_invalid_sequence(self):
        """ Testing that invalid protein sequence inputs raise the correct errors. """
        # Validation is fast: errors are raised before any descriptor computation.
#1.) sequence containing non-canonical amino acid characters raises ValueError
        with self.assertRaises(ValueError):
            protpy.all_descriptors("ABCDEF")
        with self.assertRaises(ValueError):
            protpy.all_descriptors("OOOOO")
        with self.assertRaises(ValueError):
            protpy.all_descriptors("")
#2.) non-string sequence raises TypeError
        with self.assertRaises(TypeError):
            protpy.all_descriptors(12345)
        with self.assertRaises(TypeError):
            protpy.all_descriptors(True)
        with self.assertRaises(TypeError):
            protpy.all_descriptors(["ACDEFG"])
#3.) invalid sequence with otherwise valid custom params still raises ValueError
        with self.assertRaises(ValueError):
            protpy.all_descriptors("ABCDEF", params={'kmer_composition': {'k': 3}})
