###########################################################################
#########   protPy - descriptor_set / include / exclude (Phase 8)   #######
###########################################################################

import os
import tempfile
import unittest

import protpy
from protpy.encode import _set_members
from protpy.utils import _DESCRIPTOR_REGISTRY

SEQ = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDAKLVVSTEQW"


class ProtpyDescriptorSetTests(unittest.TestCase):

    def _all(self, seq=SEQ, **kw):
        return protpy.all_descriptors(seq, _show_descriptor_list=False, **kw)

    def test_classic_is_default(self):
        self.assertTrue(self._all().equals(self._all(descriptor_set="classic")))

    def test_sets_are_nested_and_unique(self):
        classic, extended, full = self._all(), self._all(descriptor_set="extended"), self._all(descriptor_set="all")
        self.assertLess(classic.shape[1], extended.shape[1])
        self.assertLess(extended.shape[1], full.shape[1])
        self.assertEqual(list(extended.columns[:classic.shape[1]]), list(classic.columns))
        self.assertTrue(set(extended.columns) <= set(full.columns))
        for df in (classic, extended, full):
            self.assertFalse(df.columns.duplicated().any(), df.columns[df.columns.duplicated()].tolist())

    def test_set_membership_from_registry(self):
        extended = _set_members("extended")
        self.assertIn("kappa", extended)
        self.assertIn("plaac", extended)
        for name in ("cksaap", "dde", "asdc", "distance_pair", "auto_cross_covariance"):
            self.assertNotIn(name, extended)
            self.assertIn(name, _set_members("all"))
        self.assertEqual(_set_members("classic"), [])
        #extended descriptors have at most 150 columns by default
        for name in extended:
            meta = _DESCRIPTOR_REGISTRY[name]
            func = getattr(protpy, name)
            self.assertLessEqual(func(SEQ).shape[1], 150, name)

    def test_include_exclude_params(self):
        df = self._all(include=["kappa"], exclude=["tripeptide_composition", "ctd_"],
                       params={"kappa": {"blob_sizes": (5,)}})
        self.assertIn("Kappa_5", df.columns)
        self.assertNotIn("Kappa_6", df.columns)
        self.assertNotIn("AAA", df.columns)
        self.assertFalse(any(c.startswith("CTD_") for c in df.columns))
        #including a member of the set again does not duplicate it
        self.assertEqual(self._all(descriptor_set="extended", include=["kappa"]).shape,
                         self._all(descriptor_set="extended").shape)
        with self.assertRaises(ValueError):
            self._all(descriptor_set="everything")
        with self.assertRaises(ValueError):
            self._all(include=["not_a_descriptor"])
        with self.assertRaises(ValueError):
            self._all(params={"cysteine_features": {"x": 1}})

    def test_short_sequence_extended(self):
        short = SEQ[:31]
        df = self._all(short, descriptor_set="extended")
        self.assertEqual(df.shape[0], 1)
        #only documented NaNs: PLAAC LLR/COREscore for sequences shorter than the core length
        self.assertEqual(sorted(df.columns[df.isnull().any()]), ["PLAAC_COREscore", "PLAAC_LLR"])

    def test_encode_dataset_passthrough(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = os.path.join(tmp, "enc.csv")
            df = protpy.encode_dataset([SEQ, SEQ[::-1]], output_file=out, verbose=False,
                                       descriptor_set="extended", exclude=["tripeptide_composition"])
            self.assertEqual(df.shape[0], 2)
            self.assertIn("Kappa", df.columns)
            self.assertNotIn("AAA", df.columns)


if __name__ == "__main__":
    unittest.main()
