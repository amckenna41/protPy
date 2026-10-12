###########################################################################
###############      protPy - Scales Module Tests (Phase 2)   #############
###########################################################################

import unittest

import numpy as np

import protpy
from protpy._helpers import _read_json
from protpy.scales import factor_scales_

SEQ = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGD"
DIMS = {"z_scales_5": 5, "vhse": 8, "atchley_factors": 5, "t_scales": 5, "st_scales": 8,
        "protfp": 8, "fasgai": 6, "blosum_indices": 10, "ms_whim": 3}


class ProtpyScalesTests(unittest.TestCase):

    def test_tables(self):
        tables = _read_json("aa_scales.json")
        for name, d in DIMS.items():
            values = tables[name]["values"]
            self.assertEqual(sorted(values), protpy.amino_acids, name)
            self.assertTrue(all(len(v) == d for v in values.values()), name)

    def test_wrappers(self):
        tables = _read_json("aa_scales.json")
        for name, d in DIMS.items():
            func = getattr(protpy, name)
            df = func(SEQ)
            abbr = factor_scales_[name]
            self.assertEqual(df.shape, (1, d), name)
            self.assertEqual(list(df.columns), [f"{abbr}_{j}" for j in range(1, d + 1)])
            self.assertTrue(all(dt == np.float64 for dt in df.dtypes))
            #poly-A gives A's row
            np.testing.assert_allclose(func("AAAAA").values[0], np.round(tables[name]["values"]["A"], 3))
            #domains
            df3 = func(SEQ[:30], domains=3)
            self.assertEqual(df3.shape, (1, 3 * d))
            self.assertEqual(df3.columns[d], f"{abbr}_1_d2")
            #mean of the 3 equal domains is the whole-sequence mean
            np.testing.assert_allclose(df3.values[0].reshape(3, d).mean(axis=0), func(SEQ[:30]).values[0], atol=2e-3)

    def test_invalid(self):
        for func in (protpy.vhse, protpy.auto_cross_covariance, protpy.physicochemical_distance_transformation):
            with self.assertRaises(TypeError):
                func(1)
            with self.assertRaises(ValueError):
                func("")
            with self.assertRaises(ValueError):
                func("AXC")
        with self.assertRaises(ValueError):
            protpy.vhse("ACD", domains=4)
        with self.assertRaises(ValueError):
            protpy.vhse("ACD", domains=0)
        with self.assertRaises(ValueError):
            protpy.aa_scale_descriptors(SEQ, "kidera")

    def test_auto_cross_covariance(self):
        df = protpy.auto_cross_covariance(SEQ)
        self.assertEqual(df.shape, (1, 25 * 5))
        self.assertEqual(protpy.auto_cross_covariance(SEQ, mode="ac").shape, (1, 25))
        self.assertEqual(protpy.auto_cross_covariance(SEQ, mode="cc").shape, (1, 100))
        #brute force
        z = _read_json("aa_scales.json")["z_scales_5"]["values"]
        L, lag, j, k = len(SEQ), 3, 1, 4
        ref = sum(z[SEQ[i]][j] * z[SEQ[i + lag]][k] for i in range(L - lag)) / (L - lag)
        self.assertAlmostEqual(df[f"ACC_{j + 1}_{k + 1}_lag{lag}"].iloc[0], round(ref, 3))
        #AAindex codes are accepted and standardised
        self.assertEqual(protpy.auto_cross_covariance(SEQ, scale=["KYTJ820101", "HOPT810101"], lag=2).shape, (1, 8))
        with self.assertRaises(ValueError):
            protpy.auto_cross_covariance("ACD", lag=3)
        with self.assertRaises(ValueError):
            protpy.auto_cross_covariance(SEQ, mode="x")

    def test_pdt(self):
        df = protpy.physicochemical_distance_transformation(SEQ)
        self.assertEqual(df.shape, (1, 24))
        self.assertEqual(df.columns[0], "PDT_CIDH920105_1")
        self.assertTrue((protpy.physicochemical_distance_transformation("A" * 10).values == 0).all())
        self.assertTrue((df.values >= 0).all())
        with self.assertRaises(ValueError):
            protpy.physicochemical_distance_transformation("AC", lag=2)
        with self.assertRaises(ValueError):
            protpy.physicochemical_distance_transformation(SEQ, properties=[])

    def test_cross_check_peptides(self):
        """Optional cross-check against peptides.py (GPL; test-only, never imported by protpy)."""
        try:
            import peptides.tables as T
        except ImportError:
            self.skipTest("peptides not installed")
        tables = _read_json("aa_scales.json")
        for name, ref in (("ms_whim", T.MSWHIM), ("t_scales", T.T_SCALES), ("protfp", T.PROTFP)):
            keys = sorted(ref, key=lambda k: int("".join(c for c in k if c.isdigit())))
            for i, key in enumerate(keys):
                for aa, row in tables[name]["values"].items():
                    self.assertAlmostEqual(row[i], ref[key][aa], places=6, msg=f"{name} {key} {aa}")


if __name__ == "__main__":
    unittest.main()
