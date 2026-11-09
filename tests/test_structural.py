###########################################################################
###############  protPy - Structural Module Tests (Phase 7)     ###########
###########################################################################

import math
import os
import unittest

import numpy as np
from Bio import SeqIO

import protpy

SEQ = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGD"
with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "test_data", "structural_reference.fasta")) as _f:
    REF = {r.id: str(r.seq) for r in SeqIO.parse(_f, "fasta")}


class ProtpyStructuralTests(unittest.TestCase):

    def test_invalid(self):
        for name in ("aggrescan", "charge_hydropathy", "uversky_charge_hydropathy", "fold_index", "disorder_propensity",
                     "transmembrane_features", "signal_peptide_features", "functional_site_counts", "plaac",
                     "antigenicity_surface", "cysteine_features"):
            func = getattr(protpy, name)
            with self.subTest(func=name):
                with self.assertRaises(TypeError):
                    func(["A"])
                with self.assertRaises(ValueError):
                    func("")
                with self.assertRaises(ValueError):
                    func("ABX")

    def test_aggrescan(self):
        self.assertEqual(protpy.aggrescan(REF["ABETA42"]).shape, (1, 5))
        df = protpy.aggrescan(REF["ABETA42"], full=True)
        self.assertEqual(df.shape, (1, 12))
        self.assertTrue(all(c.startswith("AGG_") for c in df.columns))
        #Abeta42 has two hot spots: the central hydrophobic cluster and the C-terminus [Table 1 of the paper]
        self.assertEqual(df["AGG_nHS"].iloc[0], 2)
        self.assertAlmostEqual(df["AGG_NnHS"].iloc[0], 2 / 42 * 100)
        #poly-I is all hot spot; poly-D none; proline breaks hot spots
        self.assertEqual(protpy.aggrescan("I" * 12, full=True)["AGG_nHS"].iloc[0], 1)
        self.assertEqual(protpy.aggrescan("D" * 12, full=True)["AGG_nHS"].iloc[0], 0)
        self.assertEqual(protpy.aggrescan("IIIIIPIIIII", full=True)["AGG_nHS"].iloc[0], 2)
        self.assertEqual(protpy.aggrescan("IIIIPIIII", full=True)["AGG_nHS"].iloc[0], 0)
        #a3vSA is the mean propensity: I = 1.822
        self.assertAlmostEqual(protpy.aggrescan("I" * 6)["AGG_a3vSA"].iloc[0], 1.822, 3)
        self.assertFalse(protpy.aggrescan("ACD").isnull().any().any())
        with self.assertRaises(ValueError):
            protpy.aggrescan(SEQ, window=6)
        #psico/server window rule: 5 (<= 75), 7 (<= 175), 9 (<= 300), 11
        for length, window in ((42, 5), (100, 7), (250, 9), (400, 11)):
            seq = (SEQ * 10)[:length]
            self.assertEqual(protpy.aggrescan(seq, window="psico").values.tolist(),
                             protpy.aggrescan(seq, window=window).values.tolist())

    def test_uversky_fold_index(self):
        df = protpy.uversky_charge_hydropathy("KKKKKEEEEE")
        h = sum((-3.9 + 4.5) / 9 for _ in range(5)) / 10 + sum((-3.5 + 4.5) / 9 for _ in range(5)) / 10
        self.assertAlmostEqual(df["UV_MeanH"].iloc[0], h)
        self.assertAlmostEqual(df["UV_MeanR"].iloc[0], 0.0, 3)
        self.assertEqual(df["UV_Unfolded"].iloc[0], 1)
        self.assertEqual(protpy.uversky_charge_hydropathy("LIVFLIVF")["UV_Unfolded"].iloc[0], 0)
        #window >= L: FoldIndex equals the Uversky score
        fi = protpy.fold_index(SEQ, window=100)
        self.assertAlmostEqual(fi["FI_mean"].iloc[0], protpy.uversky_charge_hydropathy(SEQ)["UV_HQ"].iloc[0])
        self.assertAlmostEqual(protpy.charge_hydropathy(SEQ)["CH_HQ"].iloc[0], fi["FI_mean"].iloc[0])
        seq = "LIVF" * 20 + "EKSPEKSQEDKPS" * 6
        df = protpy.fold_index(seq, window=21)
        self.assertTrue(0 < df["FI_frac_disordered"].iloc[0] < 1)
        self.assertGreater(df["FI_longest_disordered"].iloc[0], 21)
        #each window's score belongs to its centre residue (v1.6.0): the disordered fraction is close
        #to the length of the disordered half, not inflated by a whole window
        self.assertAlmostEqual(df["FI_frac_disordered"].iloc[0], 78 / len(seq), delta=0.06)
        ch = protpy.charge_hydropathy(seq, window=21)
        self.assertEqual(list(ch.columns), [f"CH_w21_{s}" for s in ("mean", "max", "min", "std")])
        self.assertAlmostEqual(ch["CH_w21_mean"].iloc[0], df["FI_mean"].iloc[0])

    def test_disorder_propensity(self):
        df = protpy.disorder_propensity("WP")
        self.assertTrue(np.allclose(df.values[0], [0.5, 0.5, (-0.884 + 0.987) / 2]))

    def test_transmembrane(self):
        df = protpy.transmembrane_features(REF["GYPA_HUMAN"])
        self.assertGreaterEqual(df["TM_n_segments"].iloc[0], 1)
        self.assertGreater(df["TM_max_KD"].iloc[0], 1.6)
        self.assertEqual(protpy.transmembrane_features(SEQ)["TM_n_segments"].iloc[0], 0)
        #K/R before the segment, none after -> positive bias (and the reverse)
        tm = "S" * 10 + "L" * 22 + "S" * 10
        self.assertGreater(protpy.transmembrane_features("KRKRK" + tm + "SSSSS")["TM_positive_inside_bias"].iloc[0], 0)
        self.assertLess(protpy.transmembrane_features("SSSSS" + tm + "KRKRK")["TM_positive_inside_bias"].iloc[0], 0)

    def test_signal_peptide(self):
        df = protpy.signal_peptide_features("MKKRLLLLLLLLAAAASA" + SEQ)
        self.assertEqual(df["SP_nregion_charge"].iloc[0], 3)
        self.assertGreater(df["SP_hregion_maxKD"].iloc[0], 2.5)

    def test_functional_sites(self):
        df = protpy.functional_site_counts("NKSANPTANGTAKDEL")
        #N-glyc: NKSA yes, NPTA no (P at X), NGTA yes
        self.assertEqual(df["FS_N_glycosylation"].iloc[0], 2)
        self.assertEqual(df["FS_ER_retention"].iloc[0], 1)
        self.assertEqual(protpy.functional_site_counts("AAASKL")["FS_PTS1"].iloc[0], 1)
        self.assertEqual(protpy.functional_site_counts("AKPAKRA")["FS_tryptic_site"].iloc[0], 2)
        #a C-terminal K/R is not a cleavage site (v1.6.0)
        self.assertEqual(protpy.functional_site_counts("AKAK")["FS_tryptic_site"].iloc[0], 1)
        #the deprecated function is motif_composition with the PROSITE library
        mc = protpy.motif_composition("NKSANPTANGTAKDEL", motifs=["N_glycosylation", "ER_retention", "PTS1",
                                                              "tryptic_site"], per_100=True)
        self.assertEqual(mc.values.tolist(), df.values.tolist())
        self.assertAlmostEqual(df["FS_ER_retention_per100"].iloc[0], 100 / 16)

    def test_plaac(self):
        #reference output of PLAAC 1.1.0 for MOT3 (cli/example/MOT3-candidates-plaac-1.1.0.tsv)
        df = protpy.plaac(REF["MOT3"])
        self.assertAlmostEqual(df["PLAAC_LLR"].iloc[0], 40.751, 3)
        self.assertAlmostEqual(df["PLAAC_COREscore"].iloc[0], 40.751, 3)
        self.assertEqual(df["PLAAC_PrD_length"].iloc[0], 296)
        self.assertEqual(df["PLAAC_has_PrD"].iloc[0], 1)
        df = protpy.plaac(SEQ)
        self.assertTrue(math.isnan(df["PLAAC_LLR"].iloc[0]))
        self.assertEqual(protpy.plaac(SEQ, fill_value=0.0)["PLAAC_COREscore"].iloc[0], 0.0)
        self.assertEqual(df["PLAAC_has_PrD"].iloc[0], 0)
        self.assertEqual(protpy.plaac(REF["GYPA_HUMAN"])["PLAAC_has_PrD"].iloc[0], 0)

    def test_antigenicity(self):
        df = protpy.antigenicity_surface(SEQ)
        #exactly two hydrophobicity_profile calls
        hp = protpy.hydrophobicity_profile(SEQ, scale="parker", window=7)
        self.assertEqual(df.iloc[0, :4].tolist(), hp.iloc[0].tolist())
        self.assertEqual(list(df.columns), [f"{s}_{t}" for s in ("PARKER", "JANIN") for t in ("mean", "max", "min", "std")])
        self.assertAlmostEqual(protpy.antigenicity_surface("D" * 10)["PARKER_mean"].iloc[0], 10.0, 3)

    def test_cysteine(self):
        df = protpy.cysteine_features("CAACAAACAC")
        self.assertTrue(np.allclose(df.values[0], [4, 0.4, 1, 2, 1, 3.0]))
        self.assertAlmostEqual(protpy.cysteine_features("AAA")["Cys_mean_spacing"].iloc[0], 0.0, 3)


if __name__ == "__main__":
    unittest.main()
