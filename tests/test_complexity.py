###########################################################################
###############    protPy - Complexity Module Tests (Phase 5)  ############
###########################################################################

import math
import random
import unittest

import protpy
from protpy.complexity import _lz76

SEQ = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGD"


class ProtpyComplexityTests(unittest.TestCase):

    def _check_invalid(self, func):
        with self.assertRaises(TypeError):
            func(None)
        with self.assertRaises(ValueError):
            func("")
        with self.assertRaises(ValueError):
            func("AJX")

    def test_wootton_federhen(self):
        df = protpy.wootton_federhen_complexity(SEQ)
        self.assertEqual(list(df.columns), ["WF_K1_mean", "WF_K2_mean", "WF_K2_min", "WF_LowComplexity_frac"])
        #poly-A: K1 = K2 = 0, everything low complexity
        self.assertEqual(protpy.wootton_federhen_complexity("A" * 30).values[0].tolist(), [0.0, 0.0, 0.0, 1.0])
        #12 distinct residues in one window: K2 = log2(12), K1 = log20(12!)/12
        df = protpy.wootton_federhen_complexity("ACDEFGHIKLMN")
        self.assertEqual(df["WF_K2_mean"].iloc[0], round(math.log2(12), 3))
        self.assertEqual(df["WF_K1_mean"].iloc[0], round(math.log(math.factorial(12), 20) / 12, 3))
        self.assertEqual(df["WF_LowComplexity_frac"].iloc[0], 0.0)
        #a poly-Q tract inside a complex sequence is partly flagged
        frac = protpy.wootton_federhen_complexity(SEQ + "Q" * 20 + SEQ)["WF_LowComplexity_frac"].iloc[0]
        self.assertTrue(0 < frac < 1)
        self._check_invalid(protpy.wootton_federhen_complexity)

    def test_linguistic_complexity(self):
        random.seed(0)
        rand = "".join(random.choice(protpy.amino_acids) for _ in range(50))
        self.assertGreater(protpy.linguistic_complexity(rand).iloc[0, 0], 0.95)
        self.assertLess(protpy.linguistic_complexity("A" * 50).iloc[0, 0], 0.1)
        #"ACA": k=1 {A,C} of 2 possible (min(20,3)=3)... observed 2+2+1=5, possible 3+2+1=6
        self.assertEqual(protpy.linguistic_complexity("ACA").iloc[0, 0], round(5 / 6, 3))
        self.assertEqual(protpy.linguistic_complexity("ACA", max_k=1).iloc[0, 0], round(2 / 3, 3))
        self._check_invalid(protpy.linguistic_complexity)

    def test_lempel_ziv(self):
        #Lempel & Ziv (1976) example 0001101001000101 parses as 0.001.10.100.1000.101 -> 6
        self.assertEqual(_lz76("0001101001000101".replace("0", "A").replace("1", "C")), 6)
        self.assertEqual(_lz76("A" * 20), 2)
        df = protpy.lempel_ziv_complexity(SEQ)
        self.assertEqual(list(df.columns), ["LZ_Complexity", "LZ_Normalized"])
        c = df["LZ_Complexity"].iloc[0]
        self.assertEqual(df["LZ_Normalized"].iloc[0], round(c * math.log(len(SEQ), 20) / len(SEQ), 3))
        self.assertEqual(protpy.lempel_ziv_complexity("A").values[0].tolist(), [1, 0.0])
        self._check_invalid(protpy.lempel_ziv_complexity)

    def test_shannon_reexport(self):
        from protpy import complexity
        self.assertIs(complexity.shannon_entropy, protpy.shannon_entropy)


if __name__ == "__main__":
    unittest.main()
