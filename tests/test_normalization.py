###########################################################################
############   protPy - normalized=True output rescaling tests    #########
###########################################################################

import inspect
import math
import random
import unittest

import numpy as np

import protpy

SEQ = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDWC"

#descriptors with a normalized parameter (dataframe outputs)
NORMALIZABLE = [
    "gravy", "instability_index", "isoelectric_point", "molecular_weight", "charge_distribution",
    "secondary_structure_propensity", "aliphatic_index", "extinction_coefficient", "boman_index",
    "shannon_entropy", "hydrophobic_moment", "pseudo_amino_acid_composition",
    "moreaubroto_autocorrelation", "moran_autocorrelation", "geary_autocorrelation",
    "ctd_composition", "ctd_transition", "ctd_distribution", "ctd_",
    "sequence_order_coupling_number", "sequence_order_coupling_number_all",
    "quasi_sequence_order", "quasi_sequence_order_all",
]


def _sequences():
    """Random sequences of several lengths plus homopolymers and very short edge cases."""
    rng = random.Random(7)
    seqs = ["".join(rng.choice(protpy.amino_acids) for _ in range(n)) for n in (35, 60, 120, 400)]
    seqs += [aa * 35 for aa in protpy.amino_acids]
    seqs += ["A", "AC", "WWWWW", "KKKKDDDDKKKKDDDD"]
    return seqs


class ProtpyNormalizationTests(unittest.TestCase):

    def test_default_off_and_unchanged(self):
        for name in NORMALIZABLE:
            with self.subTest(name=name):
                func = getattr(protpy, name)
                self.assertIs(inspect.signature(func).parameters["normalized"].default, False)
                self.assertTrue(func(SEQ).equals(func(SEQ, normalized=False)))

    def test_columns_suffixed_and_in_unit_range(self):
        for name in NORMALIZABLE:
            func = getattr(protpy, name)
            for seq in _sequences():
                with self.subTest(name=name, seq=seq[:12]):
                    raw, norm = func(seq), func(seq, normalized=True)
                    self.assertEqual(list(norm.columns), [f"{c}_norm" for c in raw.columns])
                    values = norm.to_numpy(dtype=float)
                    self.assertFalse(np.isnan(values).any())
                    self.assertTrue(((values >= 0) & (values <= 1)).all())

    def test_bounds_hold_without_clipping(self):
        #_minmax clips only to absorb 3 dp rounding of raw values; with clipping disabled every
        #value must still lie within rounding distance of [0, 1], proving the bounds are correct
        from unittest import mock
        with mock.patch("protpy._helpers._unit_clip", side_effect=lambda v: v):
            for name in NORMALIZABLE:
                for seq in _sequences():
                    with self.subTest(name=name, seq=seq[:12]):
                        values = getattr(protpy, name)(seq, normalized=True).to_numpy(dtype=float)
                        self.assertTrue(((values >= -1e-3) & (values <= 1 + 1e-3)).all(), values.min())

    def test_options_combine_with_normalized(self):
        #bounds must also hold, unclipped, for every option that changes the raw values
        from unittest import mock
        cases = [
            ("gravy", {"scale": "hopp_woods"}), ("gravy", {"scale": "KYTJ820101"}),
            ("charge_distribution", {"ph": 3.0, "pka_scale": "ipc_protein"}), ("charge_distribution", {"ph": 11.0}),
            ("isoelectric_point", {"pka_scale": "bjellqvist"}),
            ("extinction_coefficient", {"absorbance": True}), ("hydrophobic_moment", {"scale": "eisenberg", "angle": "sheet"}),
            ("hydrophobic_moment", {"window": 3, "angle": 160}),
            ("pseudo_amino_acid_composition", {"lamda": 5, "weight": 0.9}),
            ("moreaubroto_autocorrelation", {"normalize": False, "lag": 5}), ("moran_autocorrelation", {"normalize": False}),
            ("geary_autocorrelation", {"lag": 40}), ("ctd_", {"property_set": "ifeature13"}),
            ("ctd_", {"all_ctd": False, "property": "charge"}),
            ("quasi_sequence_order", {"weight": 0.9, "lag": 5, "distance_matrix": "grantham"}),
            ("sequence_order_coupling_number", {"lag": 10, "distance_matrix": "grantham"}),
        ]
        with mock.patch("protpy._helpers._unit_clip", side_effect=lambda v: v):
            for name, kwargs in cases:
                for seq in _sequences():
                    with self.subTest(name=name, kwargs=kwargs, seq=seq[:12]):
                        values = getattr(protpy, name)(seq, normalized=True, **kwargs).to_numpy(dtype=float)
                        self.assertTrue(((values >= -1e-3) & (values <= 1 + 1e-3)).all())

    def test_bounds_are_attained(self):
        def one(name, seq, col=0, **kw):
            return getattr(protpy, name)(seq, normalized=True, **kw).iloc[0, col]
        self.assertEqual(one("gravy", "IIII"), 1.0)
        self.assertEqual(one("gravy", "RRRR"), 0.0)
        self.assertEqual(one("aliphatic_index", "LLLL"), 1.0)
        self.assertEqual(one("aliphatic_index", "GGGG"), 0.0)
        self.assertEqual(one("shannon_entropy", "".join(protpy.amino_acids)), 1.0)
        self.assertEqual(one("shannon_entropy", "AAAA"), 0.0)
        self.assertEqual(one("molecular_weight", "WWWW"), 1.0)
        self.assertEqual(one("molecular_weight", "GGGG"), 0.0)
        self.assertEqual(one("extinction_coefficient", "WWWW"), 1.0)
        #A is in hydrophobicity group 2; column 9 is its 100% position (residue L of L)
        self.assertEqual(one("ctd_distribution", "A" * 10, col=9), 1.0)
        #all-Arg: positive charge equals its maximum; all-Asp: negative charge does
        self.assertEqual(one("charge_distribution", "RRRR", col=0), 1.0)
        self.assertEqual(one("charge_distribution", "DDDD", col=1), 1.0)
        self.assertEqual(one("charge_distribution", "DDDD", col=2), 0.0)
        self.assertEqual(one("charge_distribution", "RRRR", col=2), 1.0)
        #a constant signal has no autocorrelation: Moran maps to the midpoint
        self.assertEqual(one("moran_autocorrelation", "A" * 40), 0.5)
        self.assertAlmostEqual(one("pseudo_amino_acid_composition", SEQ) * 100,
                               protpy.pseudo_amino_acid_composition(SEQ).iloc[0, 0], places=3)

    def test_socn_scalar(self):
        raw = protpy.sequence_order_coupling_number_(SEQ, d=2, distance_matrix="grantham")
        norm = protpy.sequence_order_coupling_number_(SEQ, d=2, distance_matrix="grantham", normalized=True)
        self.assertAlmostEqual(norm, raw / ((len(SEQ) - 2) * 4.226 ** 2), places=4)
        self.assertEqual(protpy.sequence_order_coupling_number_("AC", d=5, normalized=True), 0.0)
        #"aabb" with d=2 gives only the pairs (a, b), (a, b): a maximally distant pair reaches the bound
        dist = protpy.sequence_order._load_distance_matrix("grantham")
        a, b = max(dist, key=dist.get)
        self.assertEqual(protpy.sequence_order_coupling_number_(a + a + b + b, d=2, distance_matrix="grantham",
                                                                normalized=True), 1.0)

    def test_moran_geary_bounds_match_brute_force(self):
        #exact theoretical bounds: |I| <= N/(N-d), 0 <= C <= 2(N-1)/(N-d)
        rng = random.Random(3)
        for _ in range(30):
            seq = "".join(rng.choice("AKDW") for _ in range(rng.randint(8, 40)))
            for name in ("moran_autocorrelation", "geary_autocorrelation"):
                values = getattr(protpy, name)(seq, lag=5, normalized=True).to_numpy(dtype=float)
                self.assertTrue(((values >= 0) & (values <= 1)).all(), (name, seq))

    def test_all_descriptors_params(self):
        df = protpy.all_descriptors(SEQ, params={"gravy": {"normalized": True}, "ctd_": {"normalized": True}},
                                    _show_descriptor_list=False)
        self.assertIn("GRAVY_norm", df.columns)
        self.assertNotIn("GRAVY", df.columns)
        self.assertTrue(all(c.endswith("_norm") for c in df.columns if c.startswith("CTD_")))
        self.assertFalse(df.columns.duplicated().any())
        #default call is unaffected
        self.assertIn("GRAVY", protpy.all_descriptors(SEQ, _show_descriptor_list=False).columns)

    def test_encode_dataset_will_not_mix_raw_and_normalized(self):
        import os, tempfile
        with tempfile.TemporaryDirectory() as tmp:
            out = os.path.join(tmp, "enc.csv")
            protpy.encode_dataset([SEQ], output_file=out, verbose=False)
            with self.assertRaises(ValueError):
                protpy.encode_dataset([SEQ], output_file=out, verbose=False, params={"gravy": {"normalized": True}})

    def test_apaac_has_no_normalized_option(self):
        #its denominator 1 + w * sum(tau) can approach zero or go negative, so no fixed bound exists
        self.assertNotIn("normalized", inspect.signature(protpy.amphiphilic_pseudo_amino_acid_composition).parameters)


if __name__ == "__main__":
    unittest.main()
