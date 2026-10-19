###########################################################################
###############      protPy - Internal Helpers Tests        ###############
###########################################################################

import os
import re
import unittest

import numpy as np

import protpy
from protpy import _helpers as h
from protpy.utils import _DESCRIPTOR_REGISTRY

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


class ProtpyHelpersTests(unittest.TestCase):

    def test_load_json_table(self):
        t = h._load_json_table("reduced_alphabets.json", "ifeature5")
        self.assertEqual(set(t), {"name", "values", "reference", "doi", "source_note"})
        t["values"]["aromatic"] = "X"  #a copy: the cache must be unaffected
        self.assertEqual(h._load_json_table("reduced_alphabets.json", "ifeature5")["values"]["aromatic"], "FYW")
        with self.assertRaisesRegex(KeyError, "Valid tables"):
            h._load_json_table("reduced_alphabets.json", "nope")

    def test_every_data_table_has_provenance(self):
        data = os.path.join(ROOT, "protpy", "data")
        for fname in ("aa_scales.json", "reduced_alphabets.json", "residue_tables.json", "structural_params.json"):
            if not os.path.exists(os.path.join(data, fname)):
                continue
            for key, table in h._read_json(fname).items():
                with self.subTest(table=f"{fname}:{key}"):
                    self.assertTrue(table["reference"])
                    self.assertTrue(table["source_note"])
                    self.assertIn("doi", table)

    def test_property_matrix(self):
        m = h._property_matrix("ACA", "KYTJ820101")
        np.testing.assert_allclose(m[:, 0], [1.8, 2.5, 1.8])
        m = h._property_matrix("AC", ["KYTJ820101", "HOPT810101"])
        self.assertEqual(m.shape, (2, 2))
        m = h._property_matrix("AAA", {aa: [1.0, float(i)] for i, aa in enumerate(protpy.amino_acids)}, standardize=True)
        np.testing.assert_allclose(m[:, 0], 0)  #constant column -> 0, no division by zero
        self.assertAlmostEqual(m[0, 1], (0 - 9.5) / np.std(np.arange(20)))
        with self.assertRaises(ValueError):
            h._property_matrix("A", "NOT_A_SCALE")
        with self.assertRaises(ValueError):
            h._property_matrix("A", {"A": 1.0})

    def test_sliding_windows(self):
        np.testing.assert_allclose(h._sliding_windows([1, 2, 3, 4, 5], 3), [2, 3, 4])
        #window longer than profile -> one window
        np.testing.assert_allclose(h._sliding_windows([1, 2, 3], 9), [2])
        #edge weighting: weights 0.5,1,0.5 for window 3
        np.testing.assert_allclose(h._sliding_windows([0, 4, 0], 3, edge=0.5), [2.0])
        #window 5, edge 0.2 -> weights 0.2,0.6,1,0.6,0.2
        w = np.array([0.2, 0.6, 1, 0.6, 0.2])
        x = np.array([1., 2, 3, 4, 7])
        np.testing.assert_allclose(h._sliding_windows(x, 5, edge=0.2), [(w * x).sum() / w.sum()])

    def test_profile_stats(self):
        out = h._profile_stats(np.array([1., 2, 3]), "P", threshold=1.5)
        self.assertEqual(out, {"P_mean": 2.0, "P_max": 3.0, "P_min": 1.0, "P_std": np.std([1, 2, 3]),
                               "P_frac_above": 2 / 3})

    def test_counts(self):
        self.assertEqual(h._gapped_pair_counts("AAAA", 1)[("A", "A")], 2)
        self.assertEqual(h._gapped_pair_counts("ACDE", 0, {"A": "x", "C": "x", "D": "y", "E": "y"}),
                         {("x", "x"): 1, ("x", "y"): 1, ("y", "y"): 1})
        self.assertEqual(h._kmer_counts("AAAA", 3)[("A", "A", "A")], 2)
        self.assertEqual(sum(h._kmer_counts("A", 2).values()), 0)

    def test_charge_and_runs(self):
        np.testing.assert_array_equal(h._charge_vector("KRDEHA"), [1, 1, -1, -1, 0, 0])
        np.testing.assert_array_equal(h._charge_vector("H", his=True), [1])
        self.assertEqual(h._runs([1, 1, 0, 1, 0, 1, 1, 1]), [(0, 2), (3, 4), (5, 8)])
        self.assertEqual(h._runs([1, 1, 0, 1, 0, 1, 1, 1], min_len=2), [(0, 2), (5, 8)])
        self.assertEqual(h._runs([0, 0]), [])

    def test_to_df(self):
        df = h._to_df([1.23456, np.float64(2.0004), np.int64(3), True], ["a", "b", "c", "d"])
        self.assertEqual(df.shape, (1, 4))
        self.assertEqual(df.iloc[0].tolist(), [1.235, 2.0, 3, 1])


class ProtpyRegistryTests(unittest.TestCase):

    def test_every_descriptor_registered_with_reference(self):
        for name in protpy.descriptor_names:
            with self.subTest(name=name):
                self.assertIn(name, _DESCRIPTOR_REGISTRY)
                meta = _DESCRIPTOR_REGISTRY[name]
                self.assertTrue(meta["reference"])
                self.assertIn("doi", meta)
                self.assertIn("in_default_set", meta)

    def test_registry_parameters_match_signatures(self):
        import inspect
        for name, meta in _DESCRIPTOR_REGISTRY.items():
            with self.subTest(name=name):
                self.assertEqual(list(inspect.signature(getattr(protpy, name)).parameters), list(meta["parameters"]))

    def test_registry_dois_in_bib(self):
        with open(os.path.join(ROOT, "docs", "references.bib")) as f:
            bib = {d.lower() for d in re.findall(r"doi\s*=\s*\{([^}]*)\}", f.read(), flags=re.I)}
        for name, meta in _DESCRIPTOR_REGISTRY.items():
            if meta["doi"] is None:
                continue  #textbooks / sources without a DOI
            with self.subTest(name=name):
                self.assertTrue(meta["doi"].lower() in bib, meta["doi"])


if __name__ == "__main__":
    unittest.main()
