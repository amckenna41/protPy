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
GOLDEN = os.path.join(TEST_DATA, "golden_v1_4_2.json.gz")

#public descriptors that existed in v1.4.2; their default outputs must never change
V1_4_2_DESCRIPTORS = [
    "amino_acid_composition", "dipeptide_composition", "tripeptide_composition",
    "gravy", "aromaticity", "instability_index", "isoelectric_point", "molecular_weight",
    "charge_distribution", "hydrophobic_polar_charged_composition", "secondary_structure_propensity",
    "kmer_composition", "reduced_alphabet_composition", "motif_composition",
    "amino_acid_pair_composition", "aliphatic_index", "extinction_coefficient", "boman_index",
    "aggregation_propensity", "hydrophobic_moment",
    "moreaubroto_autocorrelation", "moran_autocorrelation", "geary_autocorrelation",
    "pseudo_amino_acid_composition", "amphiphilic_pseudo_amino_acid_composition",
    "conjoint_triad", "ctd_composition", "ctd_transition", "ctd_distribution", "ctd_",
    "sequence_order_coupling_number_", "sequence_order_coupling_number",
    "sequence_order_coupling_number_all", "quasi_sequence_order", "quasi_sequence_order_all",
    "shannon_entropy",
]


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
        out[f"{fname}::all_descriptors"] = _record(protpy.all_descriptors(seq, _show_descriptor_list=False))
        for desc in V1_4_2_DESCRIPTORS:
            out[f"{fname}::{desc}"] = _record(getattr(protpy, desc)(seq))
    return out


class ProtpyBackwardsCompatTests(unittest.TestCase):
    """Default descriptor outputs must be identical to the v1.4.2 (post-audit) baseline."""

    def test_golden_snapshot(self):
        with gzip.open(GOLDEN, "rt") as f:
            golden = json.load(f)
        current = json.loads(json.dumps(compute_snapshot()))
        self.assertEqual(sorted(golden), sorted(current))
        for key, expected in golden.items():
            with self.subTest(key=key):
                self.assertEqual(current[key], expected)


if __name__ == "__main__":
    #regenerate only when deliberately re-baselining: python tests/test_backwards_compat.py --regen
    if "--regen" in sys.argv:
        with gzip.open(GOLDEN, "wt") as f:
            json.dump(compute_snapshot(), f)
        print(f"wrote {GOLDEN}")
    else:
        unittest.main()
