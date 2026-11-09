###########################################################################
###############   protPy - Backwards Compatibility Tests    ###############
###########################################################################

import gzip
import json
import os
import sys
import unittest

from Bio import SeqIO

import protpy

TEST_DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "test_data")
GOLDEN = os.path.join(TEST_DATA, "golden_v1_6_0.json.gz")

#every registered descriptor's default output is pinned (the v1.4.2 snapshot was re-baselined in
#v1.6.0, whose audit fixes deliberately change values; see CHANGELOG.md). Regenerate only for a
#deliberate, documented change.
DESCRIPTORS = list(protpy.utils._DESCRIPTOR_REGISTRY)
KWARGS = {"aa_scale_descriptors": {"scale": "atchley_factors"}}


def _sequences() -> dict:
    seqs = {}
    for name in sorted(os.listdir(TEST_DATA)):
        #only the original v1.4.2 test sequences are pinned
        if name.startswith("test_fasta") and name.endswith(".fasta"):
            with open(os.path.join(TEST_DATA, name)) as f:
                seqs[name] = str(next(SeqIO.parse(f, "fasta")).seq)
    return seqs


def _record(df) -> dict:
    #json round-trips python floats exactly, so exact comparison is safe
    if not hasattr(df, "columns"):  #scalar-returning descriptor
        return {"type": type(df).__name__, "value": df}
    values = [df.iloc[0, i] for i in range(df.shape[1])]
    return {"columns": [str(c) for c in df.columns], "dtypes": [str(d) for d in df.dtypes],
            "values": [v.item() if hasattr(v, "item") else v for v in values], "shape": list(df.shape)}


def compute_snapshot() -> dict:
    """Compute every default output that the golden snapshot pins."""
    out = {}
    for fname, seq in _sequences().items():
        for descriptor_set in ("classic", "extended", "all"):
            out[f"{fname}::all_descriptors::{descriptor_set}"] = _record(
                protpy.all_descriptors(seq, descriptor_set=descriptor_set, _show_descriptor_list=False))
        for desc in DESCRIPTORS:
            out[f"{fname}::{desc}"] = _record(getattr(protpy, desc)(seq, **KWARGS.get(desc, {})))
    return out


class ProtpyBackwardsCompatTests(unittest.TestCase):
    """Default descriptor outputs must be identical to the v1.6.0 baseline."""

    def test_golden_snapshot(self):
        with gzip.open(GOLDEN, "rt") as f:
            golden = json.load(f)
        current = json.loads(json.dumps(compute_snapshot()))
        self.assertEqual(sorted(golden), sorted(current))
        for key, expected in golden.items():
            with self.subTest(key=key):
                #exact columns/shape; values to 1e-9 (unrounded floats can differ in the last bits across platforms)
                self.assertEqual({k: v for k, v in current[key].items() if k != "values"},
                                 {k: v for k, v in expected.items() if k != "values"})
                for a, b in zip(current[key].get("values", []), expected.get("values", [])):
                    if isinstance(b, float) and b != b:  #NaN
                        self.assertNotEqual(a, a)
                    else:
                        self.assertAlmostEqual(a, b, delta=1e-9 * max(1.0, abs(b)))


if __name__ == "__main__":
    #regenerate only when deliberately re-baselining: python tests/test_backwards_compat.py --regen
    if "--regen" in sys.argv:
        with gzip.open(GOLDEN, "wt") as f:
            json.dump(compute_snapshot(), f)
        print(f"wrote {GOLDEN}")
    else:
        unittest.main()
