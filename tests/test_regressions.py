"""Regression tests for bugs fixed in the post-1.4.2 audit."""
import os
import tempfile
import unittest

import pandas as pd
from aaindex import aaindex1

import protpy
from protpy.ctd import ctd_properties

SEQ = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGD"


class ProtpyRegressionTests(unittest.TestCase):

    def test_overlapping_counts(self):
        self.assertEqual(protpy.dipeptide_composition("AAA")["AA"].iloc[0], 100.0)
        self.assertEqual(protpy.tripeptide_composition("AAAA")["AAA"].iloc[0], 100.0)
        self.assertEqual(protpy.conjoint_triad("AAAAA")["conj_triad_111"].iloc[0], 3)

    def test_ctd_groups_partition_all_amino_acids(self):
        for name, prop in ctd_properties.items():
            groups = prop["1"] + prop["2"] + prop["3"]
            self.assertEqual(sorted(groups), sorted(protpy.amino_acids), name)

    def test_moran_matches_reference_formula(self):
        vals = {aa: v for aa, v in aaindex1["CIDH920105"].values.items() if aa in protpy.amino_acids}
        mean = sum(vals.values()) / 20
        std = (sum((v - mean) ** 2 for v in vals.values()) / 20) ** 0.5
        p = [(vals[aa] - mean) / std for aa in SEQ]
        n, pbar = len(p), sum(p) / len(p)
        var = sum((x - pbar) ** 2 for x in p) / n
        ref = sum((p[i] - pbar) * (p[i + 1] - pbar) for i in range(n - 1)) / (n - 1) / var
        got = protpy.moran_autocorrelation(SEQ, lag=1, properties=["CIDH920105"]).iloc[0, 0]
        self.assertAlmostEqual(got, ref, places=4)

    def test_short_and_repetitive_sequences_do_not_crash(self):
        funcs = [protpy.dipeptide_composition, protpy.ctd_transition, protpy.moreaubroto_autocorrelation,
                 protpy.moran_autocorrelation, protpy.geary_autocorrelation,
                 protpy.pseudo_amino_acid_composition, protpy.amphiphilic_pseudo_amino_acid_composition,
                 protpy.quasi_sequence_order]
        for seq in ["A", "AAAAAAAAAA", "ACDEFGHIKL"]:
            for func in funcs:
                df = func(seq)
                self.assertFalse(df.isnull().any().any(), f"{func.__name__}({seq})")

    def test_invalid_params_raise(self):
        with self.assertRaises(ValueError):
            protpy.moran_autocorrelation(SEQ, lag=-1)
        with self.assertRaises(ValueError):
            protpy.quasi_sequence_order(SEQ, weight="x")
        with self.assertRaises(ValueError):
            protpy.kmer_composition(SEQ, k=0)
        with self.assertRaises(ValueError):
            protpy.aggregation_propensity(SEQ, window=0)

    def test_misc_descriptor_fixes(self):
        self.assertEqual(protpy.motif_composition("DAAALLK")["dileucine_sorting"].iloc[0], 1)
        self.assertGreater(protpy.isoelectric_point("RRRRRRRRRR").iloc[0, 0], 12.0)
        self.assertEqual(protpy.molecular_weight("G").iloc[0, 0], 75.067)
        self.assertEqual(protpy.hydrophobic_moment(SEQ).shape, (1, 2))

    def test_ctd_all_ctd_param(self):
        df = protpy.all_descriptors(SEQ, params={"ctd_": {"all_ctd": False}}, _show_descriptor_list=False)
        self.assertEqual(df.filter(like="CTD_").shape[1], 21)

    def test_encode_dataset_file_safety(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "out.csv")
            pd.DataFrame({"x": [1]}).to_csv(path, index=False)
            with self.assertRaises(ValueError):
                protpy.encode_dataset([SEQ], output_file=path, verbose=False)
            self.assertEqual(list(pd.read_csv(path).columns), ["x"])  # untouched

            path = os.path.join(tmp, "enc.csv")
            protpy.encode_dataset([SEQ, SEQ], output_file=path, verbose=False)
            protpy.encode_dataset([SEQ], output_file=path, verbose=False)  # re-run seq1
            ids = pd.read_csv(path)["sequence_id"].tolist()
            self.assertEqual(sorted(ids), ["seq1", "seq2"])


if __name__ == "__main__":
    unittest.main()
