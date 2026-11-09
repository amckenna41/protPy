###########################################################################
###########   protPy - Physicochemical Module Tests (Phases 3-4)   ########
###########################################################################

import math
import unittest

import numpy as np

import protpy

SEQ = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGD"

#(Glu-Lys)25 sequence variants and kappa values from Das & Pappu 2013, Fig. 1
#(the regularly patterned variants, each 25 E + 25 K)
DAS_PAPPU_SV = {
    "sv1": ("EK" * 25, 0.0009),
    "sv2": ("EEEKKK" * 8 + "EK", 0.0025),
    "sv7": ("EEEEKKKK" * 6 + "EK", 0.0450),
    "sv8": ("KKKKEEEE" * 6 + "KE", 0.0450),
    "sv19": ("EEEEEKKKKK" * 5, 0.1941),
    "sv30": ("E" * 25 + "K" * 25, 1.0000),
}


class ProtpyChargeTests(unittest.TestCase):

    def _check_invalid(self, func):
        with self.assertRaises(TypeError):
            func(5)
        with self.assertRaises(ValueError):
            func("")
        with self.assertRaises(ValueError):
            func("AXB")

    def test_charge_fractions(self):
        df = protpy.charge_fractions("KKDD")
        self.assertEqual(list(df.columns), ["Charge_FCR", "Charge_FPlus", "Charge_FMinus", "Charge_NCPR", "Charge_Sigma"])
        self.assertEqual(df.values[0].tolist(), [1.0, 0.5, 0.5, 0.0, 0.0])
        self.assertAlmostEqual(protpy.charge_fractions("KKKA")["Charge_NCPR"].iloc[0], 0.75)
        #H is neutral; sigma = NCPR^2 / FCR
        df = protpy.charge_fractions("KKKDAAHH")
        self.assertAlmostEqual(df["Charge_FCR"].iloc[0], 0.5, 3)
        self.assertAlmostEqual(df["Charge_Sigma"].iloc[0], 0.25 ** 2 / 0.5)
        self.assertAlmostEqual(protpy.charge_fractions("AAAA")["Charge_Sigma"].iloc[0], 0.0, 3)
        self._check_invalid(protpy.charge_fractions)

    def test_ncpr(self):
        with self.assertWarns(DeprecationWarning):
            protpy.ncpr("KKKA")
        self.assertAlmostEqual(protpy.ncpr("KKKA")["NCPR"].iloc[0], 0.75, 3)
        self.assertAlmostEqual(protpy.ncpr("KRDE")["NCPR"].iloc[0], 0.0, 3)
        self._check_invalid(protpy.ncpr)

    def test_kappa_published_values(self):
        for name, (seq, ref) in DAS_PAPPU_SV.items():
            with self.subTest(sv=name):
                self.assertEqual((len(seq), seq.count("E"), seq.count("K")), (50, 25, 25))
                self.assertAlmostEqual(protpy.kappa(seq)["Kappa"].iloc[0], ref, delta=0.002)

    def test_kappa_matches_localcider_single_sign(self):
        #reference values from localCIDER 0.1.21 SequenceParameters(seq).get_kappa()
        for seq, ref in (("KKKKAAAAKKKKGGGGKKKKAAAA", 0.1724), ("MKKKAAAGSSKRRAAAPGGK", 0.2776),
                         ("DDDDAAASSSDDGGGAAAD", 0.3006)):
            with self.subTest(seq=seq):
                self.assertAlmostEqual(protpy.kappa(seq)["Kappa"].iloc[0], ref, delta=0.003)

    def test_kappa_delta_max_matches_brute_force(self):
        #fast delta_max == evaluating every segregated arrangement over the whole sequence
        from protpy.physicochemical import _delta, _delta_max
        def brute(n_pos, n_neg, L, g):
            sigma = (n_pos - n_neg) ** 2 / (L * (n_pos + n_neg))
            n0, best = L - n_pos - n_neg, 0.0
            for a in range(n0 + 1):
                for b in range(n0 - a + 1):
                    pos = np.zeros(L); neg = np.zeros(L)
                    pos[a:a + n_pos] = 1; neg[a + n_pos + b:a + n_pos + b + n_neg] = 1
                    best = max(best, float(_delta(pos, neg, g, sigma)))
            return best
        for case in ((3, 2, 12, 5), (4, 0, 10, 3), (1, 1, 9, 6), (5, 5, 14, 5), (2, 3, 7, 1)):
            with self.subTest(case=case):
                self.assertAlmostEqual(_delta_max(*case), brute(*case), 12)

    def test_kappa(self):
        df = protpy.kappa(SEQ)
        self.assertEqual(list(df.columns), ["Kappa"])
        self.assertTrue(0 <= df["Kappa"].iloc[0] <= 1)
        #documented edge cases: no charge -> NaN; one charge type is scored like localCIDER
        self.assertTrue(np.isnan(protpy.kappa("AAAAAAA")["Kappa"].iloc[0]))
        self.assertGreater(protpy.kappa("KKAAKK")["Kappa"].iloc[0], 0.0)
        #short sequence: blob reduced to the length, no crash
        self.assertFalse(protpy.kappa("KD").isnull().any().any())
        self.assertEqual(protpy.kappa(SEQ, blob_sizes=(5,)).shape, (1, 1))
        with self.assertRaises(ValueError):
            protpy.kappa(SEQ, blob_sizes=())
        self._check_invalid(protpy.kappa)

    def test_scd(self):
        self.assertAlmostEqual(protpy.sequence_charge_decoration("AAAA")["Charge_SCD"].iloc[0], 0.0, 3)
        #brute force
        q = [1 if aa in "KR" else -1 if aa in "DE" else 0 for aa in SEQ]
        ref = sum(q[m] * q[n] * math.sqrt(m - n) for m in range(len(q)) for n in range(m)) / len(q)
        self.assertAlmostEqual(protpy.sequence_charge_decoration(SEQ)["Charge_SCD"].iloc[0], ref)
        #segregated charges are more negative than mixed
        self.assertLess(protpy.sequence_charge_decoration("E" * 25 + "K" * 25)["Charge_SCD"].iloc[0],
                        protpy.sequence_charge_decoration("EK" * 25)["Charge_SCD"].iloc[0])
        self._check_invalid(protpy.sequence_charge_decoration)

    def test_charge_density(self):
        with self.assertWarns(DeprecationWarning):
            df = protpy.charge_density(SEQ, ph=7.0)
        net = protpy.charge_distribution(SEQ, ph=7.0)["NetCharge"].iloc[0]
        mw = protpy.molecular_weight(SEQ)["MolecularWeight"].iloc[0]
        self.assertAlmostEqual(df["ChargeDensity"].iloc[0], net / mw, places=6)
        self._check_invalid(protpy.charge_density)

    def test_pka_scale(self):
        #Bjellqvist is the default since v1.6.0
        self.assertEqual(protpy.isoelectric_point(SEQ).iloc[0, 0], protpy.isoelectric_point(SEQ, pka_scale="bjellqvist").iloc[0, 0])
        values = {s: protpy.isoelectric_point(SEQ, pka_scale=s).iloc[0, 0] for s in ("emboss", "protpy", "bjellqvist", "ipc_protein")}
        self.assertEqual(len(set(values.values())), 4)
        for s in ("emboss", "bjellqvist", "ipc_protein"):
            self.assertEqual(protpy.charge_distribution(SEQ, pka_scale=s).shape, (1, 3))
        with self.assertRaises(ValueError):
            protpy.isoelectric_point(SEQ, pka_scale="nope")
        with self.assertRaises(ValueError):
            protpy.charge_distribution(SEQ, pka_scale="nope")

    def test_bjellqvist_matches_biopython(self):
        from Bio.SeqUtils.IsoelectricPoint import IsoelectricPoint
        for seq in (SEQ, "ACDEFGHIKLMNPQRSTVWY", "PEPTIDEKRD", "MAAAKE"):
            with self.subTest(seq=seq):
                self.assertAlmostEqual(protpy.isoelectric_point(seq, pka_scale="bjellqvist").iloc[0, 0],
                                       IsoelectricPoint(seq).pi(), delta=0.001)


class ProtpyHydrophobicityStabilityTests(unittest.TestCase):

    def test_hydrophobicity_profile(self):
        df = protpy.hydrophobicity_profile(SEQ)
        self.assertEqual(list(df.columns), [f"HP_kyte_doolittle_w9_{s}" for s in ("mean", "max", "min", "std")])
        #segment = residues 6-45 (windows with >= 10 L); K+R in the 15 before = 3, in the 15 after = 1
        seg = protpy.hydrophobicity_profile("KKK" + "G" * 12 + "L" * 22 + "G" * 12 + "R", window=19, threshold=1.6, segments=True)
        self.assertEqual(seg.shape, (1, 8))
        self.assertEqual(seg["HP_kyte_doolittle_w19_n_segments"].iloc[0], 1)
        self.assertEqual(seg["HP_kyte_doolittle_w19_positive_inside_bias"].iloc[0], 2)
        self.assertAlmostEqual(seg["HP_kyte_doolittle_w19_frac_in_segments"].iloc[0], 40 / 50)
        with self.assertRaises(ValueError):
            protpy.hydrophobicity_profile(SEQ, segments=True)
        self.assertEqual(protpy.hydrophobicity_profile(SEQ, threshold=1.6).shape, (1, 5))
        self.assertEqual(protpy.hydrophobicity_profile(SEQ, scale="HOPT810101").shape, (1, 4))
        #window longer than sequence -> one window = GRAVY
        self.assertAlmostEqual(protpy.hydrophobicity_profile("ACD", window=9)["HP_kyte_doolittle_w9_mean"].iloc[0],
                               protpy.gravy("ACD")["GRAVY"].iloc[0])
        with self.assertRaises(ValueError):
            protpy.hydrophobicity_profile(SEQ, edge=2)
        with self.assertRaises(ValueError):
            protpy.hydrophobicity_profile(SEQ, scale="vhse")
        with self.assertRaises(ValueError):
            protpy.hydrophobicity_profile(SEQ, window=0)

    def test_hydrophobicity_profile_matches_biopython(self):
        from Bio.SeqUtils.ProtParam import ProteinAnalysis
        from Bio.SeqUtils.ProtParamData import kd
        for edge in (1.0, 0.4):
            ref = ProteinAnalysis(SEQ).protein_scale(kd, 9, edge)
            df = protpy.hydrophobicity_profile(SEQ, window=9, edge=edge)
            self.assertAlmostEqual(df["HP_kyte_doolittle_w9_mean"].iloc[0], np.mean(ref), 2)
            self.assertAlmostEqual(df["HP_kyte_doolittle_w9_max"].iloc[0], max(ref), 9)

    def test_hydrophobic_ratio(self):
        #APD example: TP10 is 61% hydrophobic (13 of 21 residues)
        self.assertAlmostEqual(protpy.hydrophobic_ratio("AGYLLGKINLKALAALAKKIL")["HydrophobicRatio"].iloc[0], 13 / 21 * 100)
        self.assertAlmostEqual(protpy.hydrophobic_ratio("GGKK")["HydrophobicRatio"].iloc[0], 0.0, 3)

    def test_half_life(self):
        df = protpy.half_life("MKT")
        self.assertEqual(df.values[0].tolist(), [30.0, 20.0, 10.0])
        self.assertTrue(np.allclose(protpy.half_life("RAA").values[0], [1.0, 2 / 60, 2 / 60]))
        self.assertTrue(np.isnan(protpy.half_life("PAA")["HalfLife_Ecoli_h"].iloc[0]))

    def test_ivywrel_cvp(self):
        df = protpy.thermostability_indices("IVYWRELA")
        self.assertEqual(list(df.columns), ["IVYWREL", "CvP_Bias"])
        self.assertAlmostEqual(df["IVYWREL"].iloc[0], 7 / 8)
        self.assertEqual(protpy.ivywrel_index("IVYWRELA")["IVYWREL"].iloc[0], 7 / 8)
        self.assertAlmostEqual(protpy.cvp_bias("DEKRNQST")["CvP_Bias"].iloc[0], 0.0, 3)
        self.assertAlmostEqual(protpy.cvp_bias("DDAA")["CvP_Bias"].iloc[0], 50.0, 3)

    def test_flexibility(self):
        df = protpy.flexibility(SEQ)
        self.assertEqual(list(df.columns), ["Flex_mean", "Flex_max", "Flex_min", "Flex_std"])
        #Biopython uses the same table and weights (but an off-by-one centre and drops the last window)
        from Bio.SeqUtils.ProtParam import ProteinAnalysis
        ref = ProteinAnalysis(SEQ).flexibility()
        self.assertAlmostEqual(df["Flex_mean"].iloc[0], np.mean(ref), delta=0.01)
        #the deprecated alias is exactly the profile engine with the Vihinen weighting
        hp = protpy.hydrophobicity_profile(SEQ, scale="vihinen_flexibility", window=9, edge=0.25)
        self.assertEqual(df.values.tolist(), hp.values.tolist())

    def test_extensions(self):
        self.assertEqual(list(protpy.gravy(SEQ).columns), ["GRAVY"])
        hw = protpy.gravy("AD", scale="hopp_woods")
        self.assertEqual(list(hw.columns), ["GRAVY_hopp_woods"])
        self.assertEqual(hw.iloc[0, 0], (-0.5 + 3.0) / 2)
        self.assertEqual(protpy.hydrophobic_moment(SEQ, angle="helix").values.tolist(),
                         protpy.hydrophobic_moment(SEQ).values.tolist())
        self.assertEqual(protpy.hydrophobic_moment(SEQ, angle="sheet").values.tolist(),
                         protpy.hydrophobic_moment(SEQ, angle=160).values.tolist())
        self.assertNotEqual(protpy.hydrophobic_moment(SEQ, scale="eisenberg_arb1984").values.tolist(),
                            protpy.hydrophobic_moment(SEQ).values.tolist())
        #the default table is available by name
        self.assertTrue(np.allclose(protpy.hydrophobic_moment(SEQ, scale="eisenberg_jmb1984").values,
                                    protpy.hydrophobic_moment(SEQ).values))
        with self.assertWarns(DeprecationWarning):
            protpy.hydrophobic_moment(SEQ, scale="eisenberg")
        with self.assertRaises(ValueError):
            protpy.hydrophobic_moment(SEQ, angle="coil")
        ext = protpy.extinction_coefficient(SEQ + "W", absorbance=True)
        self.assertEqual(ext.shape, (1, 4))
        mw = protpy.molecular_weight(SEQ + "W").iloc[0, 0]
        self.assertAlmostEqual(ext["Abs01_Reduced"].iloc[0], ext["ExtCoeff_Reduced"].iloc[0] / mw)
        self.assertEqual(protpy.extinction_coefficient(SEQ).shape, (1, 2))


if __name__ == "__main__":
    unittest.main()
