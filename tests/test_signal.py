###########################################################################
###############      protPy - Signal Module Tests (Phase 6)     ###########
###########################################################################

import math
import unittest

import numpy as np

import protpy

SEQ = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGD"
#a two-valued scale turns A/C patterns into square waves
WAVE = {aa: (1.0 if aa == "A" else -1.0) for aa in protpy.amino_acids}


class ProtpySignalTests(unittest.TestCase):

    def test_fft(self):
        df = protpy.fft_descriptors(SEQ)
        self.assertEqual(list(df.columns), [f"FFT_eiip_bin{i}" for i in range(1, 11)])
        self.assertTrue((df.values >= 0).all())
        #period-2 signal: all energy at f = 0.5 (last bin)
        df = protpy.fft_descriptors("AC" * 32, scale=WAVE)
        self.assertEqual(list(df.columns)[0], "FFT_custom_bin1")
        self.assertGreater(df["FFT_custom_bin10"].iloc[0], 0)
        self.assertEqual(df.iloc[0, :9].sum(), 0)
        #period-8 square wave: fundamental at f = 0.125 -> bin 3 holds the most power
        df = protpy.fft_descriptors("AAAACCCC" * 8, scale=WAVE)
        self.assertEqual(int(np.argmax(df.values[0])), 2)
        #constant signal has no power
        self.assertEqual(protpy.fft_descriptors("A" * 20).values.sum(), 0)
        self.assertEqual(protpy.fft_descriptors(SEQ, n_bins=4).shape, (1, 4))
        with self.assertRaises(ValueError):
            protpy.fft_descriptors(SEQ, n_bins=0)
        with self.assertRaises(ValueError):
            protpy.fft_descriptors(SEQ, scale="vhse")

    def test_dct(self):
        df = protpy.dct_descriptors(SEQ)
        self.assertEqual(list(df.columns), [f"DCT_kyte_doolittle_{i}" for i in range(10)])
        #explicit O(L^2) DCT-II
        kd = protpy.composition.kyte_doolittle_
        x = np.array([kd[aa] for aa in SEQ]); x = x - x.mean()
        L = len(x)
        for k in range(10):
            ref = sum(x[n] * math.cos(math.pi * k * (2 * n + 1) / (2 * L)) for n in range(L))
            ref *= math.sqrt((1 if k == 0 else 2) / L)
            self.assertAlmostEqual(df.iloc[0, k], round(ref, 3), places=3)
        #orthonormal: all L coefficients preserve the signal energy
        full = protpy.dct_descriptors(SEQ, n_coeffs=L)
        self.assertAlmostEqual((full.values ** 2).sum(), (x ** 2).sum(), delta=0.05)
        with self.assertRaises(ValueError):
            protpy.dct_descriptors("ACD", n_coeffs=4)

    def test_invalid(self):
        for func in (protpy.fft_descriptors, protpy.dct_descriptors):
            with self.assertRaises(TypeError):
                func(3.0)
            with self.assertRaises(ValueError):
                func("")
            with self.assertRaises(ValueError):
                func("AZ")


if __name__ == "__main__":
    unittest.main()
