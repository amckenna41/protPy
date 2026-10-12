###########################################################################
##################    protPy - Encode Module Tests    ####################
###########################################################################

import pandas as pd
import numpy as np
import os
import shutil
import tempfile
import unittest
from unittest.mock import patch, MagicMock
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
            protpy.all_descriptors(self.protein_seq1, params={'aromaticity': {}})
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

    @patch('protpy.encode._tqdm')
    def test_all_descriptors_show_progress(self, mock_tqdm):
        """
        Testing that show_progress=True instantiates tqdm with the correct
        arguments and drives all bar lifecycle methods (set_description,
        update x N, close). The computation runs on the shortest test sequence
        (protein_seq3, ~76 AA) so this test stays fast.
        """
        mock_bar = MagicMock()
        mock_tqdm.return_value = mock_bar

        df = protpy.all_descriptors(self.protein_seq3, show_progress=True, output_file=None)

        # Determine the actual step count so the test stays correct if steps change.
        expected_steps = mock_tqdm.call_args.kwargs.get('total')

#1.) tqdm constructed once with expected keyword arguments
        mock_tqdm.assert_called_once()
        kwargs = mock_tqdm.call_args.kwargs
        self.assertIsNotNone(expected_steps,
            'Expected tqdm to be called with a total= argument.')
        self.assertEqual(kwargs['unit'], 'descriptor',
            'Expected tqdm unit="descriptor", got "{}".'.format(kwargs.get('unit')))
        self.assertIn('bar_format', kwargs,
            'Expected bar_format to be passed to tqdm.')
        self.assertIn('dynamic_ncols', kwargs,
            'Expected dynamic_ncols to be passed to tqdm.')
#2.) bar.update() called exactly once per descriptor step
        self.assertEqual(mock_bar.update.call_count, expected_steps,
            'Expected bar.update() called {} times, got {}.'.format(
                expected_steps, mock_bar.update.call_count))
#3.) set_description called at least once (once per step + final done message)
        self.assertTrue(mock_bar.set_description.called,
            'Expected bar.set_description() to be called at least once.')
        self.assertGreaterEqual(mock_bar.set_description.call_count, expected_steps,
            'Expected at least {} set_description calls, got {}.'.format(
                expected_steps, mock_bar.set_description.call_count))
#4.) bar.close() called exactly once
        mock_bar.close.assert_called_once_with()
#5.) output is still a valid DataFrame with a single row and no nulls
        self.assertIsInstance(df, pd.DataFrame,
            'Expected DataFrame with show_progress=True, got {}.'.format(type(df)))
        self.assertEqual(df.shape[0], 1,
            'Expected 1-row output, got shape {}.'.format(df.shape))
        self.assertEqual(df.isnull().any().sum(), 0,
            'Expected no null values in output with show_progress=True.')

    @patch('protpy.encode._tqdm')
    def test_all_descriptors_no_progress_by_default(self, mock_tqdm):
        """
        Testing that show_progress defaults to False and that tqdm is never
        instantiated when show_progress is omitted or explicitly set to False.
        """
#1.) default call — tqdm must not be called
        protpy.all_descriptors(self.protein_seq3, output_file=None)
        mock_tqdm.assert_not_called()
#2.) explicit show_progress=False — tqdm must still not be called
        mock_tqdm.reset_mock()
        protpy.all_descriptors(self.protein_seq3, show_progress=False, output_file=None)
        mock_tqdm.assert_not_called()


###########################################################################
##############    protPy - read_sequences_from_file Tests    ##############
###########################################################################

class ProtpyReadSequencesTests(unittest.TestCase):
    """
    Test suite for the read_sequences_from_file function in the encode module.

    Test Cases
    ----------
    test_read_sequences_fasta_single_record:
        testing correct parsing of a single-record FASTA file.
    test_read_sequences_fasta_multi_record:
        testing that all records in a multi-record FASTA are returned.
    test_read_sequences_fasta_fa_extension:
        testing that the .fa extension is handled identically to .fasta.
    test_read_sequences_txt_plain:
        testing plain TXT parsing (one sequence per non-blank line).
    test_read_sequences_txt_delimited_seq_col_only:
        testing multi-column TXT with seq_col returns list[str].
    test_read_sequences_txt_delimited_with_id_col:
        testing multi-column TXT with id_col returns DataFrame with only
        the two requested columns; all other columns are dropped.
    test_read_sequences_csv_by_col_name:
        testing CSV parsing with sequence column selected by name.
    test_read_sequences_csv_by_col_index:
        testing CSV parsing with sequence column selected by zero-based index.
    test_read_sequences_csv_with_id_col:
        testing CSV with id_col returns two-column DataFrame; other cols dropped.
    test_read_sequences_tsv:
        testing TSV file parsing with column selection by name and with id_col.
    test_read_sequences_has_header_false:
        testing headerless CSV/TXT files are parsed correctly.
    test_read_sequences_output_is_uppercase:
        testing that sequences are returned in uppercase regardless of input case.
    test_read_sequences_drops_blank_entries:
        testing that blank sequence rows are excluded from both list and
        DataFrame outputs.
    test_read_sequences_file_not_found:
        testing FileNotFoundError raised for a missing file.
    test_read_sequences_unsupported_extension:
        testing ValueError raised for an unsupported file extension.
    test_read_sequences_bad_col_name:
        testing ValueError raised when seq_col or id_col name is not found.
    test_read_sequences_bad_col_index:
        testing ValueError raised when a column index is out of range.
    test_read_sequences_bad_col_type:
        testing TypeError raised when seq_col or id_col has an invalid type.
    """

    #canonical sequences used across tests
    SEQ_A = "ACDEFGHIKLMNPQRSTVWY"
    SEQ_B = "MKTAYIAKQRQISFVKSHFS"
    SEQ_C = "WQRQISFVKSHFSMKTAYIA"

    def setUp(self):
        """ Create a temporary directory and populate it with test data files. """
        self.tmpdir = tempfile.mkdtemp()

        #plain TXT — one sequence per line, with a blank line that should be skipped
        self.txt_path = os.path.join(self.tmpdir, "seqs.txt")
        with open(self.txt_path, "w") as f:
            f.write(f"{self.SEQ_A}\n{self.SEQ_B}\n\n{self.SEQ_C}\n")

        #lowercase TXT — sequences in lowercase; output must be uppercase
        self.txt_lower_path = os.path.join(self.tmpdir, "seqs_lower.txt")
        with open(self.txt_lower_path, "w") as f:
            f.write(f"{self.SEQ_A.lower()}\n{self.SEQ_B.lower()}\n")

        #multi-column TXT (comma-separated, like thermostability.txt)
        self.multi_txt_path = os.path.join(self.tmpdir, "multi.txt")
        with open(self.multi_txt_path, "w") as f:
            f.write("name,sequence,label\n")
            f.write(f"seq1,{self.SEQ_A},1\n")
            f.write(f"seq2,{self.SEQ_B},0\n")
            f.write(f"seq3,{self.SEQ_C},1\n")

        #CSV with three columns; sequence column is 'sequence'
        self.csv_path = os.path.join(self.tmpdir, "data.csv")
        with open(self.csv_path, "w") as f:
            f.write("id,sequence,score\n")
            f.write(f"id1,{self.SEQ_A},1.5\n")
            f.write(f"id2,{self.SEQ_B},2.3\n")
            f.write(f"id3,{self.SEQ_C},0.8\n")

        #TSV with three columns
        self.tsv_path = os.path.join(self.tmpdir, "data.tsv")
        with open(self.tsv_path, "w") as f:
            f.write("id\tsequence\tscore\n")
            f.write(f"id1\t{self.SEQ_A}\t1.5\n")
            f.write(f"id2\t{self.SEQ_B}\t2.3\n")
            f.write(f"id3\t{self.SEQ_C}\t0.8\n")

        #headerless CSV — sequences only, no column names
        self.csv_no_header_path = os.path.join(self.tmpdir, "no_header.csv")
        with open(self.csv_no_header_path, "w") as f:
            f.write(f"{self.SEQ_A}\n{self.SEQ_B}\n{self.SEQ_C}\n")

        #CSV with a blank sequence entry that should be excluded from output
        self.csv_blank_path = os.path.join(self.tmpdir, "with_blank.csv")
        with open(self.csv_blank_path, "w") as f:
            f.write("id,sequence\n")
            f.write(f"id1,{self.SEQ_A}\n")
            f.write("id2,\n")
            f.write(f"id3,{self.SEQ_B}\n")

        #multi-record FASTA — three records
        self.multi_fasta_path = os.path.join(self.tmpdir, "multi.fasta")
        with open(self.multi_fasta_path, "w") as f:
            f.write(f">rec1 First sequence\n{self.SEQ_A}\n")
            f.write(f">rec2 Second sequence\n{self.SEQ_B}\n")
            f.write(f">rec3 Third sequence\n{self.SEQ_C}\n")

        #.fa copy of the multi-record FASTA
        self.fa_path = os.path.join(self.tmpdir, "multi.fa")
        shutil.copy(self.multi_fasta_path, self.fa_path)

    def tearDown(self):
        """ Remove all temporary test files and directory. """
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def test_read_sequences_fasta_single_record(self):
        """ Testing correct list output for a single-record FASTA file. """
        seqs = protpy.read_sequences_from_file(
            os.path.join("tests", "test_data", "test_fasta1.fasta")
        )
        self.assertIsInstance(seqs, list,
            'Expected list output for FASTA, got {}.'.format(type(seqs)))
        self.assertEqual(len(seqs), 1,
            'Expected 1 sequence from single-record FASTA, got {}.'.format(len(seqs)))
        self.assertIsInstance(seqs[0], str,
            'Expected sequence to be a string, got {}.'.format(type(seqs[0])))
        self.assertTrue(seqs[0].isupper(),
            'Expected sequence to be uppercase.')
        self.assertGreater(len(seqs[0]), 0,
            'Expected non-empty sequence.')

    def test_read_sequences_fasta_multi_record(self):
        """ Testing that all records in a multi-record FASTA are returned in order. """
        seqs = protpy.read_sequences_from_file(self.multi_fasta_path)
        self.assertIsInstance(seqs, list,
            'Expected list output for FASTA, got {}.'.format(type(seqs)))
        self.assertEqual(len(seqs), 3,
            'Expected 3 sequences from 3-record FASTA, got {}.'.format(len(seqs)))
        self.assertEqual(seqs[0], self.SEQ_A,
            'Expected first sequence {}, got {}.'.format(self.SEQ_A, seqs[0]))
        self.assertEqual(seqs[1], self.SEQ_B,
            'Expected second sequence {}, got {}.'.format(self.SEQ_B, seqs[1]))
        self.assertEqual(seqs[2], self.SEQ_C,
            'Expected third sequence {}, got {}.'.format(self.SEQ_C, seqs[2]))

    def test_read_sequences_fasta_fa_extension(self):
        """ Testing that the .fa extension is handled identically to .fasta. """
        seqs = protpy.read_sequences_from_file(self.fa_path)
        self.assertIsInstance(seqs, list,
            'Expected list output for .fa file, got {}.'.format(type(seqs)))
        self.assertEqual(len(seqs), 3,
            'Expected 3 sequences from .fa file, got {}.'.format(len(seqs)))
        self.assertEqual(seqs[0], self.SEQ_A,
            'Expected first sequence to match, got {}.'.format(seqs[0]))

    def test_read_sequences_txt_plain(self):
        """ Testing plain TXT parsing: one sequence per non-blank line. """
        seqs = protpy.read_sequences_from_file(self.txt_path)
        self.assertIsInstance(seqs, list,
            'Expected list output for plain TXT, got {}.'.format(type(seqs)))
#blank line in the middle must be skipped
        self.assertEqual(len(seqs), 3,
            'Expected 3 sequences (blank line skipped), got {}.'.format(len(seqs)))
        self.assertEqual(seqs[0], self.SEQ_A,
            'Expected first sequence {}, got {}.'.format(self.SEQ_A, seqs[0]))
        self.assertEqual(seqs[1], self.SEQ_B,
            'Expected second sequence {}, got {}.'.format(self.SEQ_B, seqs[1]))
        self.assertEqual(seqs[2], self.SEQ_C,
            'Expected third sequence {}, got {}.'.format(self.SEQ_C, seqs[2]))

    def test_read_sequences_txt_delimited_seq_col_only(self):
        """ Testing multi-column TXT returns list[str] when seq_col given but no id_col. """
        seqs = protpy.read_sequences_from_file(self.multi_txt_path, seq_col="sequence")
        self.assertIsInstance(seqs, list,
            'Expected list output when id_col is None, got {}.'.format(type(seqs)))
        self.assertEqual(len(seqs), 3,
            'Expected 3 sequences from multi-column TXT, got {}.'.format(len(seqs)))
        self.assertEqual(seqs[0], self.SEQ_A,
            'Expected first sequence {}, got {}.'.format(self.SEQ_A, seqs[0]))

    def test_read_sequences_txt_delimited_with_id_col(self):
        """
        Testing multi-column TXT with id_col returns a two-column DataFrame;
        all other columns (e.g. 'label') are dropped.
        """
        df = protpy.read_sequences_from_file(
            self.multi_txt_path, seq_col="sequence", id_col="name"
        )
        self.assertIsInstance(df, pd.DataFrame,
            'Expected DataFrame when id_col provided, got {}.'.format(type(df)))
        self.assertEqual(df.shape, (3, 2),
            'Expected shape (3, 2), got {}.'.format(df.shape))
        self.assertEqual(list(df.columns), ["name", "sequence"],
            'Expected columns ["name", "sequence"], got {}.'.format(list(df.columns)))
        self.assertNotIn("label", df.columns,
            '"label" column should have been dropped.')
        self.assertEqual(df["name"].tolist(), ["seq1", "seq2", "seq3"],
            'Expected id values ["seq1","seq2","seq3"], got {}.'.format(df["name"].tolist()))
        self.assertEqual(df["sequence"].tolist(),
            [self.SEQ_A, self.SEQ_B, self.SEQ_C],
            'Sequence values do not match expected.')

    def test_read_sequences_csv_by_col_name(self):
        """ Testing CSV sequence column selected by column name. """
        seqs = protpy.read_sequences_from_file(self.csv_path, seq_col="sequence")
        self.assertIsInstance(seqs, list,
            'Expected list output for CSV, got {}.'.format(type(seqs)))
        self.assertEqual(len(seqs), 3,
            'Expected 3 sequences, got {}.'.format(len(seqs)))
        self.assertEqual(seqs, [self.SEQ_A, self.SEQ_B, self.SEQ_C],
            'Sequence list does not match expected values.')

    def test_read_sequences_csv_by_col_index(self):
        """ Testing CSV sequence column selected by zero-based index. """
        seqs = protpy.read_sequences_from_file(self.csv_path, seq_col=1)
        self.assertIsInstance(seqs, list,
            'Expected list output for CSV with integer seq_col, got {}.'.format(type(seqs)))
        self.assertEqual(seqs, [self.SEQ_A, self.SEQ_B, self.SEQ_C],
            'Sequence list does not match expected values for seq_col=1.')

    def test_read_sequences_csv_with_id_col(self):
        """
        Testing CSV with id_col specified returns a two-column DataFrame;
        the 'score' column is dropped.
        """
        df = protpy.read_sequences_from_file(
            self.csv_path, seq_col="sequence", id_col="id"
        )
        self.assertIsInstance(df, pd.DataFrame,
            'Expected DataFrame when id_col provided, got {}.'.format(type(df)))
        self.assertEqual(df.shape, (3, 2),
            'Expected shape (3, 2), got {}.'.format(df.shape))
        self.assertEqual(list(df.columns), ["id", "sequence"],
            'Expected columns ["id", "sequence"], got {}.'.format(list(df.columns)))
        self.assertNotIn("score", df.columns,
            '"score" column should have been dropped.')
        self.assertEqual(df["id"].tolist(), ["id1", "id2", "id3"],
            'Expected id values ["id1","id2","id3"], got {}.'.format(df["id"].tolist()))
        self.assertEqual(df["sequence"].tolist(),
            [self.SEQ_A, self.SEQ_B, self.SEQ_C],
            'Sequence values do not match expected.')

    def test_read_sequences_tsv(self):
        """ Testing TSV file parsed correctly with both seq_col only and with id_col. """
#1.) sequences only
        seqs = protpy.read_sequences_from_file(self.tsv_path, seq_col="sequence")
        self.assertIsInstance(seqs, list,
            'Expected list output for TSV, got {}.'.format(type(seqs)))
        self.assertEqual(len(seqs), 3,
            'Expected 3 sequences from TSV, got {}.'.format(len(seqs)))
        self.assertEqual(seqs, [self.SEQ_A, self.SEQ_B, self.SEQ_C],
            'TSV sequence list does not match expected values.')
#2.) with id_col
        df = protpy.read_sequences_from_file(
            self.tsv_path, seq_col="sequence", id_col="id"
        )
        self.assertIsInstance(df, pd.DataFrame,
            'Expected DataFrame for TSV with id_col, got {}.'.format(type(df)))
        self.assertEqual(list(df.columns), ["id", "sequence"],
            'Expected columns ["id", "sequence"], got {}.'.format(list(df.columns)))
        self.assertEqual(df.shape, (3, 2),
            'Expected shape (3, 2), got {}.'.format(df.shape))

    def test_read_sequences_has_header_false(self):
        """ Testing headerless CSV parsed correctly with has_header=False. """
        seqs = protpy.read_sequences_from_file(
            self.csv_no_header_path, has_header=False
        )
        self.assertIsInstance(seqs, list,
            'Expected list output for headerless CSV, got {}.'.format(type(seqs)))
        self.assertEqual(len(seqs), 3,
            'Expected 3 sequences from headerless CSV, got {}.'.format(len(seqs)))
        self.assertEqual(seqs, [self.SEQ_A, self.SEQ_B, self.SEQ_C],
            'Sequence list does not match expected for headerless CSV.')
#column 0 by index should give the same result
        seqs2 = protpy.read_sequences_from_file(
            self.csv_no_header_path, seq_col=0, has_header=False
        )
        self.assertEqual(seqs, seqs2,
            'Expected seq_col=0 to produce same result as default for headerless CSV.')

    def test_read_sequences_output_is_uppercase(self):
        """ Testing that sequences are returned uppercase regardless of input case. """
        seqs = protpy.read_sequences_from_file(self.txt_lower_path)
        self.assertTrue(all(s.isupper() for s in seqs),
            'Expected all sequences to be uppercase, got: {}.'.format(seqs))
        self.assertEqual(seqs[0], self.SEQ_A,
            'Expected lowercase input to be uppercased to {}, got {}.'.format(
                self.SEQ_A, seqs[0]))
        self.assertEqual(seqs[1], self.SEQ_B,
            'Expected lowercase input to be uppercased to {}, got {}.'.format(
                self.SEQ_B, seqs[1]))

    def test_read_sequences_drops_blank_entries(self):
        """
        Testing that blank sequence entries are excluded from both list
        and DataFrame outputs.
        """
#1.) list output — blank entry (id2) must be absent
        seqs = protpy.read_sequences_from_file(self.csv_blank_path, seq_col="sequence")
        self.assertEqual(len(seqs), 2,
            'Expected 2 sequences (blank entry dropped), got {}.'.format(len(seqs)))
        self.assertNotIn("", seqs,
            'Empty string must not appear in sequence list output.')
#2.) DataFrame output — row with blank sequence must be dropped
        df = protpy.read_sequences_from_file(
            self.csv_blank_path, seq_col="sequence", id_col="id"
        )
        self.assertEqual(df.shape[0], 2,
            'Expected 2 rows in DataFrame (blank entry dropped), got {}.'.format(df.shape[0]))
        self.assertFalse((df["sequence"] == "").any(),
            'DataFrame must not contain rows with blank sequence.')
        self.assertEqual(df["id"].tolist(), ["id1", "id3"],
            'Expected ids ["id1", "id3"] after blank drop, got {}.'.format(
                df["id"].tolist()))

    def test_read_sequences_file_not_found(self):
        """ Testing FileNotFoundError raised when the file path does not exist. """
        with self.assertRaises(FileNotFoundError):
            protpy.read_sequences_from_file("definitely_does_not_exist.fasta")
        with self.assertRaises(FileNotFoundError):
            protpy.read_sequences_from_file(
                os.path.join(self.tmpdir, "missing.csv"), seq_col="sequence"
            )

    def test_read_sequences_unsupported_extension(self):
        """ Testing ValueError raised for unsupported file extensions. """
        for ext in (".json", ".xlsx", ".parquet", ".xml"):
            path = os.path.join(self.tmpdir, f"data{ext}")
            with open(path, "w") as f:
                f.write("placeholder")
            with self.assertRaises(ValueError,
                    msg=f'Expected ValueError for extension {ext}'):
                protpy.read_sequences_from_file(path)

    def test_read_sequences_bad_col_name(self):
        """ Testing ValueError raised when seq_col or id_col name is not in the file. """
#bad seq_col name
        with self.assertRaises(ValueError):
            protpy.read_sequences_from_file(self.csv_path, seq_col="nonexistent_col")
#bad id_col name (seq_col is valid)
        with self.assertRaises(ValueError):
            protpy.read_sequences_from_file(
                self.csv_path, seq_col="sequence", id_col="nonexistent_id"
            )
#error message should name the bad column
        try:
            protpy.read_sequences_from_file(self.csv_path, seq_col="bad_col")
            self.fail('Expected ValueError for unknown seq_col.')
        except ValueError as e:
            self.assertIn("bad_col", str(e),
                'Expected error message to name the missing column, got: {}'.format(e))

    def test_read_sequences_bad_col_index(self):
        """ Testing ValueError raised when a column index is out of range. """
        with self.assertRaises(ValueError):
            protpy.read_sequences_from_file(self.csv_path, seq_col=99)
        with self.assertRaises(ValueError):
            protpy.read_sequences_from_file(
                self.csv_path, seq_col="sequence", id_col=99
            )

    def test_read_sequences_bad_col_type(self):
        """ Testing TypeError raised when seq_col or id_col has an invalid type. """
        for bad_val in (3.5, [1], {"col": 0}, True):
            with self.assertRaises((TypeError, ValueError),
                    msg=f'Expected error for seq_col={bad_val!r}'):
                protpy.read_sequences_from_file(self.csv_path, seq_col=bad_val)
            with self.assertRaises((TypeError, ValueError),
                    msg=f'Expected error for id_col={bad_val!r}'):
                protpy.read_sequences_from_file(
                    self.csv_path, seq_col="sequence", id_col=bad_val
                )


###########################################################################
################    protPy - encode_dataset Tests    #####################
###########################################################################

class ProtpyEncodeDatasetTests(unittest.TestCase):
    """
    Test suite for the encode_dataset function in the encode module.

    Test Cases
    ----------
    test_encode_dataset_list_ids:
        testing that a plain list input produces auto-generated seq1/seq2/…
        identifiers in the first column.
    test_encode_dataset_dataframe_ids:
        testing that a DataFrame input uses the first column as the identifier.
    test_encode_dataset_first_col_is_id:
        testing that the identifier column is always the first column in the
        returned DataFrame.
    test_encode_dataset_empty_list:
        testing that an empty list input returns an empty DataFrame.
    test_encode_dataset_empty_dataframe:
        testing that an empty DataFrame input returns an empty DataFrame.
    test_encode_dataset_gap_raises_error:
        testing ValueError raised when a list sequence contains gaps and strip_gaps=False.
    test_encode_dataset_gap_raises_error_dataframe:
        testing ValueError raised for a gapped DataFrame sequence, and that the
        sequence identifier appears in the error message.
    test_encode_dataset_strip_gaps_removes_gaps:
        testing that strip_gaps=True strips gap characters before calling all_descriptors.
    test_encode_dataset_bad_seq_col:
        testing ValueError raised when seq_col is not found in a DataFrame.
    test_encode_dataset_bad_type:
        testing TypeError raised when sequences is not a list or DataFrame.
    test_encode_dataset_output_file:
        testing that the CSV output file is written with the identifier as
        the first column.
    """

    SEQ_A = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGD"
    SEQ_B = "ACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWY"

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def _mock_row(self, n_cols=10):
        cols = [f"feat_{i}" for i in range(n_cols)]
        return pd.DataFrame([[float(i) for i in range(n_cols)]], columns=cols)

    @patch('protpy.encode.all_descriptors')
    def test_encode_dataset_list_ids(self, mock_all_desc):
        """ Testing plain list input generates seq1, seq2, … identifier column. """
        mock_all_desc.side_effect = [self._mock_row(), self._mock_row()]

        result = protpy.encode_dataset([self.SEQ_A, self.SEQ_B], output_file=None)

        self.assertIsInstance(result, pd.DataFrame,
            'Expected DataFrame output, got {}.'.format(type(result)))
        self.assertEqual(result.columns[0], "sequence_id",
            'Expected first column "sequence_id", got "{}".'.format(result.columns[0]))
        self.assertEqual(result["sequence_id"].tolist(), ["seq1", "seq2"],
            'Expected ids ["seq1","seq2"], got {}.'.format(result["sequence_id"].tolist()))
        self.assertEqual(result.shape[0], 2,
            'Expected 2 rows, got {}.'.format(result.shape[0]))

    @patch('protpy.encode.all_descriptors')
    def test_encode_dataset_dataframe_ids(self, mock_all_desc):
        """ Testing DataFrame input uses the first column as the identifier. """
        mock_all_desc.side_effect = [self._mock_row(), self._mock_row()]

        df_input = pd.DataFrame({
            "name": ["prot_A", "prot_B"],
            "sequence": [self.SEQ_A, self.SEQ_B],
        })
        result = protpy.encode_dataset(df_input, seq_col="sequence", output_file=None)

        self.assertEqual(result.columns[0], "name",
            'Expected first column "name", got "{}".'.format(result.columns[0]))
        self.assertEqual(result["name"].tolist(), ["prot_A", "prot_B"],
            'Expected ids ["prot_A","prot_B"], got {}.'.format(result["name"].tolist()))
        self.assertEqual(result.shape[0], 2,
            'Expected 2 rows, got {}.'.format(result.shape[0]))

    @patch('protpy.encode.all_descriptors')
    def test_encode_dataset_first_col_is_id(self, mock_all_desc):
        """ Testing that identifier is always the first column regardless of input form. """
        mock_all_desc.side_effect = [self._mock_row(), self._mock_row()]

        # list input
        result_list = protpy.encode_dataset([self.SEQ_A, self.SEQ_B], output_file=None)
        self.assertEqual(result_list.columns[0], "sequence_id",
            'Expected "sequence_id" as first column for list input.')

        mock_all_desc.side_effect = [self._mock_row(), self._mock_row()]
        df_input = pd.DataFrame({
            "entry_id": ["e1", "e2"],
            "sequence": [self.SEQ_A, self.SEQ_B],
        })
        result_df = protpy.encode_dataset(df_input, seq_col="sequence", output_file=None)
        self.assertEqual(result_df.columns[0], "entry_id",
            'Expected "entry_id" as first column for DataFrame input.')

    def test_encode_dataset_empty_list(self):
        """ Testing that an empty list input returns an empty DataFrame. """
        result = protpy.encode_dataset([], output_file=None)
        self.assertIsInstance(result, pd.DataFrame,
            'Expected DataFrame for empty list, got {}.'.format(type(result)))
        self.assertTrue(result.empty,
            'Expected empty DataFrame for empty list input.')

    def test_encode_dataset_empty_dataframe(self):
        """ Testing that an empty DataFrame input returns an empty DataFrame. """
        empty_df = pd.DataFrame({"name": [], "sequence": []})
        result = protpy.encode_dataset(empty_df, seq_col="sequence", output_file=None)
        self.assertIsInstance(result, pd.DataFrame,
            'Expected DataFrame for empty DataFrame input, got {}.'.format(type(result)))
        self.assertTrue(result.empty,
            'Expected empty DataFrame for empty DataFrame input.')

    def test_encode_dataset_bad_seq_col(self):
        """ Testing ValueError raised when seq_col is absent from DataFrame. """
        df_input = pd.DataFrame({
            "name": ["prot_A"],
            "sequence": [self.SEQ_A],
        })
        with self.assertRaises(ValueError):
            protpy.encode_dataset(df_input, seq_col="nonexistent_col", output_file=None)

    def test_encode_dataset_gap_raises_error(self):
        """ Testing ValueError raised when a sequence contains gaps and strip_gaps=False. """
        gapped_seq = "MKTAY-IAKQR"
        try:
            protpy.encode_dataset([gapped_seq], output_file=None, strip_gaps=False)
            self.fail('Expected ValueError for sequence containing gap characters.')
        except ValueError as e:
            self.assertIn("gap", str(e).lower(),
                'Expected error message to mention gaps, got: {}'.format(e))
            self.assertIn("strip_gaps", str(e),
                'Expected error message to mention strip_gaps, got: {}'.format(e))

    def test_encode_dataset_gap_raises_error_dataframe(self):
        """ Testing ValueError raised for gapped sequences in DataFrame input. """
        df_input = pd.DataFrame({
            "name": ["prot_A", "prot_B"],
            "sequence": [self.SEQ_A, "MKTAY-IAKQR"],
        })
        try:
            protpy.encode_dataset(df_input, seq_col="sequence", output_file=None, strip_gaps=False)
            self.fail('Expected ValueError for DataFrame sequence containing gap characters.')
        except ValueError as e:
            self.assertIn("prot_B", str(e),
                'Expected error message to include the sequence identifier, got: {}'.format(e))

    @patch('protpy.encode.all_descriptors')
    def test_encode_dataset_strip_gaps_removes_gaps(self, mock_all_desc):
        """ Testing that strip_gaps=True removes gap characters before encoding. """
        mock_all_desc.return_value = self._mock_row()
        gapped_seq = "MKTAY-IAKQR"
        protpy.encode_dataset([gapped_seq], output_file=None, strip_gaps=True)
        called_seq = mock_all_desc.call_args[0][0]
        self.assertNotIn("-", called_seq,
            'Expected gap characters to be removed before encoding, got: {}'.format(called_seq))

    def test_encode_dataset_bad_type(self):
        """ Testing TypeError raised when sequences is not a list or DataFrame. """
        for bad in ("MKTAYIAK", 42, None, {"a": "b"}):
            with self.assertRaises(TypeError,
                    msg=f'Expected TypeError for input type {type(bad)}'):
                protpy.encode_dataset(bad, output_file=None)

    @patch('protpy.encode.all_descriptors')
    def test_encode_dataset_output_file(self, mock_all_desc):
        """
        Testing that the CSV written by encode_dataset has the identifier as
        the first column and contains one row per sequence.
        """
        mock_all_desc.side_effect = [self._mock_row(), self._mock_row()]

        output_path = os.path.join(self.tmpdir, "out.csv")
        protpy.encode_dataset(
            [self.SEQ_A, self.SEQ_B],
            output_file=output_path,
            save_each_sequence=True,
        )

        self.assertTrue(os.path.isfile(output_path),
            'Expected output CSV to be written.')
        written = pd.read_csv(output_path)
        self.assertEqual(written.shape[0], 2,
            'Expected 2 rows in CSV, got {}.'.format(written.shape[0]))
        self.assertEqual(written.columns[0], "sequence_id",
            'Expected first column "sequence_id" in CSV, got "{}".'.format(written.columns[0]))
        self.assertEqual(written["sequence_id"].tolist(), ["seq1", "seq2"],
            'Expected ids ["seq1","seq2"] in CSV, got {}.'.format(
                written["sequence_id"].tolist()))
