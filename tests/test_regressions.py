"""Regression tests for bugs fixed in the post-1.4.2 audit."""
import os
import random
import tempfile
import unittest

import pandas as pd
from aaindex import aaindex1

import protpy
from protpy.ctd import ctd_properties
from protpy.complexity import _lz76
from protpy._helpers import _scale_values

SEQ = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGD"


class ProtpyRegressionTests(unittest.TestCase):

    def test_overlapping_counts(self):
        self.assertAlmostEqual(protpy.dipeptide_composition("AAA")["AA"].iloc[0], 100.0, 3)
        self.assertAlmostEqual(protpy.tripeptide_composition("AAAA")["AAA"].iloc[0], 100.0, 3)
        self.assertEqual(protpy.conjoint_triad("AAAAA", normalize=False)["conj_triad_111"].iloc[0], 3)

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
        lagged = {protpy.moreaubroto_autocorrelation: "lag", protpy.moran_autocorrelation: "lag",
                  protpy.geary_autocorrelation: "lag", protpy.pseudo_amino_acid_composition: "lamda",
                  protpy.amphiphilic_pseudo_amino_acid_composition: "lamda", protpy.quasi_sequence_order: "lag"}
        for seq in ["A", "AAAAAAAAAA", "ACDEFGHIKL"]:
            for func in (protpy.dipeptide_composition, protpy.ctd_transition):
                self.assertFalse(func(seq).isnull().any().any(), f"{func.__name__}({seq})")
            for func, name in lagged.items():
                #v1.6.0 lag policy: a lag that is not less than the length raises instead of silently becoming 30
                with self.assertRaises(ValueError):
                    func(seq, **{name: len(seq)})
                if len(seq) > 1:
                    df = func(seq, **{name: len(seq) - 1})
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
        self.assertAlmostEqual(protpy.molecular_weight("G").iloc[0, 0], 75.067, 3)
        self.assertEqual(protpy.hydrophobic_moment(SEQ).shape, (1, 2))

    def test_ctd_all_ctd_param(self):
        df = protpy.all_descriptors(SEQ, params={"ctd": {"all_ctd": False}}, _show_descriptor_list=False)
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


class ProtpyFinalReviewRegressionTests(unittest.TestCase):
    """Bugs fixed in the v1.6.0 final review."""

    def test_lz76_matches_kaspar_schuster(self):
        def kaspar_schuster(s):
            #reference implementation (Kaspar & Schuster 1987)
            n = len(s)
            if n == 1:
                return 1
            c, i, k, l, k_max = 1, 0, 1, 1, 1
            while True:
                if s[i + k - 1] == s[l + k - 1]:
                    k += 1
                    if l + k > n:
                        return c + 1
                else:
                    k_max = max(k, k_max)
                    i += 1
                    if i == l:
                        c += 1
                        l += k_max
                        if l + 1 > n:
                            return c
                        i, k, k_max = 0, 1, 1
                    else:
                        k = 1
        rng = random.Random(0)
        seqs = ["A", "AB", "AAAA", "ABAB" * 10, "0001101001000101"]
        seqs += ["".join(rng.choice(a) for _ in range(rng.randint(1, 300))) for a in ("AB", "ACDEFGHIKLMNPQRSTVWY") for _ in range(200)]
        for s in seqs:
            self.assertEqual(_lz76(s), kaspar_schuster(s), s)
        #the published example: 0001101001000101 parses as 0.001.10.100.1000.101 (6 components)
        self.assertEqual(_lz76("0001101001000101"), 6)

    def test_numeric_parameter_validation(self):
        with self.assertRaisesRegex(ValueError, "must be an integer"):
            protpy.kmer_composition(SEQ, k=2.5)
        with self.assertRaisesRegex(ValueError, "must be an integer"):
            protpy.autocorrelation(SEQ, lag=2.7)
        self.assertEqual(protpy.kmer_composition(SEQ, k=2.0).shape, (1, 400))
        for weight in (float("nan"), float("inf"), -1):
            with self.assertRaisesRegex(ValueError, "finite non-negative"):
                protpy.pseudo_amino_acid_composition(SEQ, lamda=5, weight=weight)
        for ph in (float("nan"), -5, 15, "7"):
            with self.assertRaisesRegex(ValueError, "between 0 and 14"):
                protpy.charge_distribution(SEQ, ph=ph)

    def test_invalid_motif_regex(self):
        with self.assertRaisesRegex(ValueError, "Invalid regex pattern"):
            protpy.motif_composition(SEQ, motifs=["(("])

    def test_scale_values_returns_a_copy(self):
        values = _scale_values("kyte_doolittle")
        before = protpy.gravy(SEQ).iloc[0, 0]
        values["A"] = 999.0
        self.assertEqual(protpy.gravy(SEQ).iloc[0, 0], before)

    def test_encode_dataset_rejects_non_string_sequence(self):
        with self.assertRaisesRegex(TypeError, "Sequence 'seq2' must be a string"):
            protpy.encode_dataset([SEQ, None], output_file=None, verbose=False)

    def test_fasta_ids_stop_codon_and_csv_escaping(self):
        with tempfile.TemporaryDirectory() as tmp:
            fasta = os.path.join(tmp, "t.fasta")
            with open(fasta, "w") as fh:
                fh.write(">sp|P1|A one\r\n" + SEQ + "*\r\n>empty\r\n>=HYPERLINK(1)\r\n" + SEQ.lower() + "\r\n")
            self.assertEqual(protpy.read_sequences_from_file(fasta), [SEQ, SEQ])
            df = protpy.read_sequences_from_file(fasta, id_col="name")
            self.assertEqual(list(df.columns), ["name", "sequence"])
            self.assertEqual(df["name"].tolist(), ["sp|P1|A", "=HYPERLINK(1)"])
            out = os.path.join(tmp, "out.csv")
            for _ in range(2):  #a re-run replaces the rows instead of duplicating them
                result = protpy.encode_dataset(df, output_file=out, verbose=False)
            self.assertEqual(result["name"].tolist(), ["sp|P1|A", "=HYPERLINK(1)"])
            with open(out) as fh:
                first_cells = [line.split(",")[0] for line in fh.read().splitlines()]
            self.assertEqual(first_cells, ["name", "sp|P1|A", "'=HYPERLINK(1)"])


if __name__ == "__main__":
    unittest.main()
