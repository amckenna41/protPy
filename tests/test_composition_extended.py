###########################################################################
#########   protPy - v1.5.0 Composition Descriptor Tests (Phase 1)  #######
###########################################################################

import itertools
import math
import unittest

import numpy as np

import protpy
from protpy.ctd import ctd_property_sets

SEQ = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQ"
AA = protpy.amino_acids


class ProtpyCompositionExtendedTests(unittest.TestCase):

    def _check_invalid(self, func, **bad):
        with self.assertRaises(TypeError):
            func(123)
        with self.assertRaises(ValueError):
            func("")
        with self.assertRaises(ValueError):
            func("ACDXZ")
        for k, v in bad.items():
            with self.assertRaises(ValueError, msg=f"{func.__name__}({k}={v})"):
                func(SEQ, **{k: v})

    def test_cksaap(self):
        df = protpy.cksaap(SEQ)
        self.assertEqual(df.shape, (1, 2400))
        self.assertEqual(df.columns[0], "CKSAAP_AA_gap0")
        self.assertEqual(df.columns[-1], "CKSAAP_YY_gap5")
        self.assertTrue(all(dt == np.float64 for dt in df.dtypes))
        #gap=0 equals the dipeptide composition (which is rounded to 2 dp)
        dpc = protpy.dipeptide_composition(SEQ)
        gap0 = protpy.cksaap(SEQ, gap=0)
        np.testing.assert_allclose(gap0.values[0], dpc.values[0], atol=0.006)
        self.assertEqual(protpy.cksaap("AAAA", gap=1)["CKSAAP_AA_gap1"].iloc[0], 100.0)
        #spacing with no pairs -> zeros
        self.assertEqual(protpy.cksaap("AC", gap=2).filter(like="_gap2").values.sum(), 0)
        self._check_invalid(protpy.cksaap, gap=-1)

    def test_cksaagp(self):
        df = protpy.cksaagp(SEQ)
        self.assertEqual(df.shape, (1, 150))
        self.assertEqual(df.columns[0], "CKSAAGP_aliphatic.aliphatic_gap0")
        for g in range(6):
            self.assertAlmostEqual(df.filter(like=f"_gap{g}").values.sum(), 100, delta=0.05)
        #G and A are both aliphatic
        self.assertEqual(protpy.cksaagp("GAG", gap=0)["CKSAAGP_aliphatic.aliphatic_gap0"].iloc[0], 100.0)
        self._check_invalid(protpy.cksaagp, gap="x")

    def test_dde(self):
        df = protpy.dde("ACDEFGHIKL")
        self.assertEqual(df.shape, (1, 400))
        self.assertEqual(df.columns[1], "DDE_AC")
        #hand computation of the AC cell: one AC dipeptide, A has 4 codons, C has 2
        dc, tm = 1 / 9, (4 / 61) * (2 / 61)
        expected = (dc - tm) / math.sqrt(tm * (1 - tm) / 9)
        self.assertAlmostEqual(df["DDE_AC"].iloc[0], round(expected, 3))
        #an absent dipeptide is negative
        self.assertLess(df["DDE_AA"].iloc[0], 0)
        with self.assertRaises(ValueError):
            protpy.dde("A")
        self._check_invalid(protpy.dde)

    def test_asdc(self):
        df = protpy.asdc("ACA")
        self.assertEqual(df.shape, (1, 400))
        for col in ("ASDC_AC", "ASDC_AA", "ASDC_CA"):
            self.assertEqual(df[col].iloc[0], 33.333)
        self.assertAlmostEqual(protpy.asdc(SEQ).values.sum(), 100, delta=0.05)
        #brute force on SEQ
        counts = {}
        for i, j in itertools.combinations(range(len(SEQ)), 2):
            counts[SEQ[i] + SEQ[j]] = counts.get(SEQ[i] + SEQ[j], 0) + 1
        total = len(SEQ) * (len(SEQ) - 1) / 2
        self.assertEqual(protpy.asdc(SEQ)["ASDC_KV"].iloc[0], round(counts["KV"] / total * 100, 3))
        self.assertEqual(protpy.asdc("A").values.sum(), 0)
        self._check_invalid(protpy.asdc)

    def test_distance_pair(self):
        df = protpy.distance_pair(SEQ)
        self.assertEqual(df.shape, (1, 14 + 14 * 14 * 3))
        self.assertEqual(df.columns[0], "DP_IMV")
        self.assertIn("DP_IMV.L_d1", df.columns)
        self.assertEqual(protpy.distance_pair(SEQ, rng=3, alphabet="cp13").shape, (1, 13 + 169 * 3))
        self.assertEqual(protpy.distance_pair(SEQ, rng=1, alphabet="cp20").shape, (1, 420))
        #composition block sums to 100; each distance block sums to 100
        self.assertAlmostEqual(df.iloc[0, :14].sum(), 100, delta=0.05)
        self.assertAlmostEqual(df.filter(like="_d2").values.sum(), 100, delta=0.05)
        #I, M and V share a cp14 cluster
        self.assertEqual(protpy.distance_pair("IMV", rng=2)["DP_IMV.IMV_d2"].iloc[0], 100.0)
        self._check_invalid(protpy.distance_pair, rng=-1, alphabet="cp15")

    def test_reduced_alphabets_partition(self):
        from protpy.composition import _named_alphabet
        for name in ("ifeature5", "cp13", "cp14", "cp19", "cp20"):
            residues = "".join(g for _, g in _named_alphabet(name))
            self.assertEqual(sorted(residues), AA, name)
        self.assertEqual(len(_named_alphabet("cp14")), 14)

    def test_natural_vector(self):
        df = protpy.natural_vector("ACA")
        self.assertEqual(df.shape, (1, 60))
        #A at positions 1 and 3: n=2, mu=2, D2 = ((1-2)^2 + (3-2)^2) / (2*3)
        self.assertEqual(df["NV_n_A"].iloc[0], 2)
        self.assertEqual(df["NV_mu_A"].iloc[0], 2.0)
        self.assertEqual(df["NV_D2_A"].iloc[0], round(2 / 6, 3))
        self.assertEqual(df["NV_mu_C"].iloc[0], 2.0)
        self.assertEqual(df["NV_D2_C"].iloc[0], 0.0)
        self.assertEqual(df["NV_mu_W"].iloc[0], 0.0)
        self._check_invalid(protpy.natural_vector)

    def test_atomic_composition(self):
        df = protpy.atomic_composition("GG", subtract_water=True)
        self.assertEqual(list(df.columns), ["ATC_C", "ATC_H", "ATC_N", "ATC_O", "ATC_S"])
        #C4H8N2O3 -> 17 atoms
        np.testing.assert_allclose(df.values[0], np.round(np.array([4, 8, 2, 3, 0]) / 17 * 100, 3))
        #default: free amino acids, no water removed (Pfeature convention): 2 x C2H5NO2
        np.testing.assert_allclose(protpy.atomic_composition("GG").values[0], [20, 50, 10, 20, 0])
        self.assertGreater(protpy.atomic_composition("MC")["ATC_S"].iloc[0], 0)
        self.assertAlmostEqual(protpy.atomic_composition(SEQ).values.sum(), 100, delta=0.01)
        self._check_invalid(protpy.atomic_composition)

    #extensions to existing functions

    def test_reduced_alphabet_extension(self):
        #defaults unchanged
        self.assertEqual(list(protpy.reduced_alphabet_composition(SEQ).columns), [f"ReducedAlphabet_{i}" for i in range(1, 7)])
        gaac = protpy.reduced_alphabet_composition("GAFKD", alphabet="ifeature5")
        self.assertEqual(list(gaac.columns), ["GAAC_aliphatic", "GAAC_aromatic", "GAAC_positive", "GAAC_negative", "GAAC_uncharged"])
        self.assertEqual(gaac.values[0].tolist(), [40.0, 20.0, 20.0, 20.0, 0.0])
        self.assertEqual(protpy.reduced_alphabet_composition(SEQ, alphabet="ifeature5", k=2).shape, (1, 25))
        gtpc = protpy.reduced_alphabet_composition(SEQ, alphabet="ifeature5", k=3)
        self.assertEqual(gtpc.shape, (1, 125))
        self.assertEqual(gtpc.columns[0], "GTPC_aliphatic.aliphatic.aliphatic")
        self.assertAlmostEqual(gtpc.values.sum(), 100, delta=0.1)
        with self.assertRaises(ValueError):
            protpy.reduced_alphabet_composition(SEQ, k=4)
        with self.assertRaises(ValueError):
            protpy.reduced_alphabet_composition(SEQ, alphabet="nope")

    def test_ksctriad(self):
        self.assertEqual(protpy.conjoint_triad(SEQ).shape, (1, 343))
        df = protpy.conjoint_triad("AAAAAAA", gap=1)
        self.assertEqual(df.shape, (1, 686))
        #gap 0: 5 contiguous triads; gap 1: (i, i+2, i+4) -> 3 triads
        self.assertEqual(df["conj_triad_111_gap0"].iloc[0], 5)
        self.assertEqual(df["conj_triad_111_gap1"].iloc[0], 3)
        with self.assertRaises(ValueError):
            protpy.conjoint_triad(SEQ, gap=-1)

    def test_ctd_ifeature13(self):
        for name, prop in ctd_property_sets["ifeature13"].items():
            groups = prop["1"] + prop["2"] + prop["3"]
            self.assertEqual(sorted(groups), AA, name)
        self.assertEqual(len(ctd_property_sets["ifeature13"]), 13)
        self.assertEqual(protpy.ctd_(SEQ).shape, (1, 147))
        self.assertEqual(protpy.ctd_(SEQ, property_set="ifeature13").shape, (1, 273))
        self.assertEqual(protpy.ctd_composition(SEQ, property="hydrophobicity_ENGD860101").shape, (1, 3))
        with self.assertRaises(ValueError):
            protpy.ctd_(SEQ, property_set="nope")

    def test_region(self):
        seq = "A" * 25 + "C" * 10 + "D" * 25
        df = protpy.amino_acid_composition(seq, region="split")
        self.assertEqual(df.shape, (1, 60))
        self.assertEqual((df["NT_A"].iloc[0], df["MID_C"].iloc[0], df["CT_D"].iloc[0]), (100.0, 100.0, 100.0))
        self.assertEqual(protpy.amino_acid_composition(seq, region="nterm", n=5)["NT_A"].iloc[0], 100.0)
        self.assertEqual(protpy.dipeptide_composition(seq, region="cterm").shape, (1, 400))
        self.assertEqual(protpy.dipeptide_composition(seq, region="split").columns[400], "MID_AA")
        with self.assertRaises(ValueError):
            protpy.amino_acid_composition("A" * 50, region="split")
        with self.assertRaises(ValueError):
            protpy.amino_acid_composition(seq, region="middle")


if __name__ == "__main__":
    unittest.main()
