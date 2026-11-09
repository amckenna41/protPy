###########################################################################
###########   protPy - Reference Values and Identities (v1.6.0)   ########
###########################################################################
"""
Checks against independent references rather than protpy's own earlier outputs
(the v1.6.0 descriptor audit found that frozen-output tests had pinned bugs, e.g.
the instability-index table and the QSO x100 error).
"""

import random
import unittest
import warnings

import numpy as np
from Bio.SeqUtils.ProtParam import ProteinAnalysis

import protpy

random.seed(0)
SEQS = ["".join(random.choice(protpy.amino_acids) for _ in range(random.randint(40, 300))) for _ in range(25)]


class ProtpyReferenceTests(unittest.TestCase):

    def test_protparam_matches_biopython(self):
        for seq in SEQS:
            pa = ProteinAnalysis(seq)
            with self.subTest(seq=seq[:10]):
                #Biopython's DIWV table has the typo FY = 33.601 (protpy: 33.60)
                self.assertAlmostEqual(protpy.instability_index(seq).iloc[0, 0], pa.instability_index(), delta=1e-3)
                self.assertAlmostEqual(protpy.isoelectric_point(seq).iloc[0, 0], pa.isoelectric_point(), delta=1e-3)
                self.assertAlmostEqual(protpy.gravy(seq).iloc[0, 0], pa.gravy(), 9)
                self.assertAlmostEqual(protpy.aromaticity(seq).iloc[0, 0], pa.aromaticity(), 9)

    def test_composition_sums(self):
        for seq in SEQS[:5]:
            self.assertAlmostEqual(protpy.quasi_sequence_order(seq).values.sum(), 1.0, 9)
            self.assertAlmostEqual(protpy.pseudo_amino_acid_composition(seq).values.sum(), 100.0, 9)
            self.assertAlmostEqual(protpy.amphiphilic_pseudo_amino_acid_composition(seq).values.sum(), 100.0, 9)
            self.assertAlmostEqual(protpy.kmer_composition(seq, k=1, alphabet=3).values.sum(), 100.0, 9)

    def test_audit_identities(self):
        seq = SEQS[0]
        #Moreau-Broto == auto covariance (mode="ac") on the same standardised properties
        props = ["CIDH920105", "BHAR880101"]
        mb = protpy.moreaubroto_autocorrelation(seq, lag=5, properties=props).values[0]
        acc = protpy.auto_cross_covariance(seq, scale=props, lag=5, mode="ac").values[0]
        #MB is property-major, ACC lag-major
        self.assertTrue(np.allclose(mb, acc.reshape(5, 2).T.ravel()))
        #PDT == 2 s^2 Geary (s^2 = sample variance of the standardised property over the sequence)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            pdt = protpy.physicochemical_distance_transformation(seq, lag=3, properties=["CIDH920105"]).values[0]
        geary = protpy.geary_autocorrelation(seq, lag=3, properties=["CIDH920105"]).values[0]
        from protpy._helpers import _property_matrix
        s2 = _property_matrix(seq, ["CIDH920105"], standardize=True)[:, 0].var(ddof=1)
        self.assertTrue(np.allclose(pdt, 2 * s2 * geary))
        #DPC == kmer(k=2) == CKSAAP gap 0
        dpc = protpy.dipeptide_composition(seq).values[0]
        self.assertTrue(np.allclose(dpc, protpy.kmer_composition(seq, k=2).values[0]))
        self.assertTrue(np.allclose(dpc, protpy.cksaap(seq, gap=0).values[0]))

    def test_ctd_distribution_quantiles_are_ordered(self):
        for seq in SEQS + ["AKAAAAAAAAAAKAAAAAAAAAAAK", "KAAA", "AAAK"]:
            d = protpy.ctd(seq, parts="D").values[0].reshape(-1, 5)
            for row in d:
                if row.any():
                    self.assertTrue(np.all(np.diff(row) >= 0), (seq[:10], row))

    def test_kappa_das_pappu(self):
        #published kappa values (Das & Pappu 2013, Fig. 1)
        for seq, ref in (("EK" * 25, 0.0009), ("EEEEEKKKKK" * 5, 0.1941), ("E" * 25 + "K" * 25, 1.0)):
            self.assertAlmostEqual(protpy.kappa(seq).iloc[0, 0], ref, delta=0.002)


if __name__ == "__main__":
    unittest.main()
