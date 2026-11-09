###########################################################################
#########   protPy - descriptor_set / include / exclude (Phase 8)   #######
###########################################################################

import os
import tempfile
import unittest

import numpy as np

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
        classic, extended, full = (_set_members(s) for s in ("classic", "extended", "all"))
        #nested sets
        self.assertTrue(set(classic) < set(extended) < set(full))
        self.assertIn("geary_autocorrelation", classic)
        self.assertIn("kappa", extended)
        self.assertIn("plaac", extended)
        for name in ("cksaap", "asdc", "distance_pair", "auto_cross_covariance",
                     "moreaubroto_autocorrelation", "moran_autocorrelation", "sequence_order_coupling_number"):
            self.assertNotIn(name, extended)
            self.assertIn(name, full)
        #opt-in only: TPC (8000 columns), DDE, the kmer engine, aliases of members
        for name in ("tripeptide_composition", "dde", "kmer_composition", "ctd_", "aa_scale_descriptors"):
            self.assertNotIn(name, full)
        #extended descriptors have at most 150 columns by default
        for name in set(extended) - set(classic):
            self.assertLessEqual(getattr(protpy, name)(SEQ).shape[1], 150, name)

    def test_classic_bundle_contents(self):
        df = self._all()
        #AAC 20 + DPC 400 + 14 scalars + PAAC 50 + APAAC 80 + Geary 240 + CTriad 343 + CTD 147 + QSO 100 + motifs 8
        self.assertEqual(df.shape, (1, 1402))
        self.assertIn("ExtCoeff_Oxidized", df.columns)
        self.assertNotIn("ExtCoeff_Reduced", df.columns)
        self.assertEqual(df.filter(like="QSO_").shape[1], 100)
        self.assertFalse(any(c.startswith(("MBAuto_", "MAuto_", "SOCN_")) for c in df.columns))
        self.assertNotIn("AAA", df.columns)
        #rounding is opt-in
        rounded = self._all(decimals=2)
        self.assertTrue((rounded.round(2) == rounded).all().all())

    def test_include_exclude_params(self):
        df = self._all(include=["kappa", "tripeptide_composition"], exclude=["ctd"],
                       params={"kappa": {"blob_sizes": (5,)}, "extinction_coefficient": {}})
        self.assertIn("Kappa", df.columns)
        self.assertIn("AAA", df.columns)
        #user params for a bundled descriptor give its full output
        self.assertIn("ExtCoeff_Reduced", df.columns)
        self.assertFalse(any(c.startswith("CTD_") for c in df.columns))
        #including a member of the set again does not duplicate it
        self.assertEqual(self._all(descriptor_set="extended", include=["kappa"]).shape,
                         self._all(descriptor_set="extended").shape)
        with self.assertRaises(ValueError):
            self._all(descriptor_set="everything")
        with self.assertRaises(ValueError):
            self._all(include=["not_a_descriptor"])
        with self.assertRaises(ValueError):
            self._all(params={"aromaticity": {"x": 1}})
        #deprecated functions are not registered, so they cannot be included
        with self.assertRaises(ValueError):
            self._all(include=["ncpr"])

    def test_short_sequence_extended(self):
        short = SEQ[:31]
        df = self._all(short, descriptor_set="extended")
        self.assertEqual(df.shape[0], 1)
        #only documented NaNs: PLAAC LLR/COREscore for sequences shorter than the core length
        self.assertEqual(sorted(df.columns[df.isnull().any()]), ["PLAAC_COREscore", "PLAAC_LLR"])
        #a sequence too short for the default lags is computed with lag = length - 1, the
        #higher-lag columns are NaN and the columns match a long sequence (v1.6.0)
        full = self._all(SEQ, descriptor_set="all")
        for length in (1, 2, 5, 20, 30):
            with self.subTest(length=length):
                with self.assertWarnsRegex(UserWarning, "too short for the lags"):
                    df = self._all(SEQ[:length], descriptor_set="all")
                self.assertEqual(list(df.columns), list(full.columns))
        with self.assertWarns(UserWarning):
            df = self._all(SEQ[:20])
        geary = protpy.geary_autocorrelation(SEQ[:20], lag=19)
        self.assertTrue(np.allclose(df[geary.columns].values, geary.values))
        self.assertTrue(df[[f"GAuto_CIDH920105_{d}" for d in range(20, 31)]].isna().all(axis=None))
        #smaller lags avoid the NaN columns
        lags = {"pseudo_amino_acid_composition": {"lamda": 5}, "amphiphilic_pseudo_amino_acid_composition": {"lamda": 5, "weight": 0.1},
                "geary_autocorrelation": {"lag": 5}, "quasi_sequence_order": {"lag": 5}}
        self.assertFalse(self._all(SEQ[:20], params=lags).isna().any(axis=None))
        #invalid parameters still raise for a short sequence
        with self.assertRaises(ValueError):
            self._all(SEQ[:20], params={"geary_autocorrelation": {"lag": -1}})
        #mixed lengths in one dataset share the columns
        X = protpy.encode_dataset([SEQ[:8], SEQ], output_file=None, verbose=False)
        self.assertEqual(X.shape, (2, 1 + self._all(SEQ).shape[1]))

    def test_params_include_exclude_checks(self):
        #params for a descriptor that is not computed warn; a non-dict value raises
        with self.assertWarnsRegex(UserWarning, "not computed"):
            self._all(params={"kappa": {"blob_sizes": (5,)}})
        with self.assertRaises(TypeError):
            self._all(params={"pseudo_amino_acid_composition": 3})
        #a single name may be passed as a str
        self.assertIn("Kappa", self._all(include="kappa").columns)
        self.assertNotIn("GRAVY", self._all(exclude="gravy").columns)

    def test_registered_descriptors_on_very_short_sequences(self):
        #every registered descriptor either works or raises a ValueError on tiny inputs
        #(lag/lamda/gap too large for the sequence, or a length requirement)
        for name in _DESCRIPTOR_REGISTRY:
            func = getattr(protpy, name)
            kwargs = {"scale": "atchley_factors"} if name == "aa_scale_descriptors" else {}
            for seq in ("A", "AC", "MKTAYIAK"):
                with self.subTest(name=name, seq=seq):
                    try:
                        df = func(seq, **kwargs)
                    except ValueError as e:
                        self.assertRegex(str(e), "less than the sequence length|at least 2 residues|exceed")
                        continue
                    self.assertEqual(df.shape[0], 1)

    def test_encode_dataset_passthrough(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = os.path.join(tmp, "enc.csv")
            df = protpy.encode_dataset([SEQ, SEQ[::-1]], output_file=out, verbose=False,
                                       descriptor_set="extended", include=["tripeptide_composition"],
                                       decimals=3)
            self.assertEqual(df.shape[0], 2)
            self.assertIn("Kappa", df.columns)
            self.assertIn("AAA", df.columns)


if __name__ == "__main__":
    unittest.main()
