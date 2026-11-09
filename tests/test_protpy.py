###########################################################################
#################            protPy Module Tests          #################
###########################################################################

import unittest
unittest.TestLoader.sortTestMethodsUsing = None

import protpy as protpy

class ProtpyTests(unittest.TestCase):
    """
    Test suite for testing protpy package and its variables and metadata.

    Test Cases
    ----------
    test_protpy_metadata:
        testing pypi software package metadata.
    test_valid_descriptors:
        testing correct list of valid descriptors.
    """
    def test_protpy_metadata(self):
        """ Testing correct protpy version and metadata. """
        self.assertEqual(protpy.__version__, "1.6.0",
            "protpy version is not correct, expected version 1.6.0, got version {}".format(protpy.__version__))
        self.assertEqual(protpy.__name__, "protpy", 
            "protpy software name is not correct, expected protpy, got {}.".format(protpy.__name__))
        self.assertEqual(protpy.__author__, "AJ McKenna: https://github.com/amckenna41", 
            "protpy author is not correct, expected AJ McKenna, got {}.".format(protpy.__author__))
        self.assertEqual(protpy.__authorEmail__, "amckenna41@qub.ac.uk", 
            "protpy author email is not correct, expected amckenna41@qub.ac.uk, got {}.".format(protpy.__authorEmail__))
        self.assertEqual(protpy.__url__, "https://github.com/amckenna41/protPy", 
            "protpy repo URL is not correct, expected https://github.com/amckenna41/protPy, got {}.".format(protpy.__url__))
        self.assertEqual(protpy.__download_url__, "https://github.com/amckenna41/protPy/archive/refs/heads/main.zip", 
            "protpy repo download URL is not correct, expected https://github.com/amckenna41/protPy/archive/refs/heads/main.zip, got {}.".format(protpy.__download_url__))
        self.assertEqual(protpy.__status__, "Production", 
            "protpy status is not correct, expected Production got {}.".format(protpy.__status__))
        self.assertEqual(protpy.__license__, "MIT", 
            "protpy license type is not correct, expected MIT, got {}.".format(protpy.__license__))
        self.assertEqual(protpy.__maintainer__, "AJ McKenna", 
            "protpy maintainer is not correct, expected AJ McKenna, got {}.".format(protpy.__license__))
        self.assertEqual(protpy.__keywords__, ["bioinformatics", "protein engineering", 
            "python", "pypi", "machine learning", "aaindex", "protein descriptors", 
            "physicochemical descriptors", "biochemical descriptors", "structural descriptors", "pySAR"], 
            "protpy keywords is not correct, got\n{}.".format(protpy.__keywords__))

    def test_valid_descriptors(self):
        """ Testing the registered descriptor names and the deprecated v1.6.0 functions. """
        registered = ["amino_acid_composition", "dipeptide_composition", "tripeptide_composition",
            "gravy", "aromaticity", "instability_index", "isoelectric_point", "molecular_weight",
            "charge_distribution", "kmer_composition", "motif_composition", "aliphatic_index",
            "extinction_coefficient", "boman_index", "hydrophobic_moment", "shannon_entropy",
            "autocorrelation", "moreaubroto_autocorrelation", "moran_autocorrelation", "geary_autocorrelation",
            "pseudo_amino_acid_composition", "amphiphilic_pseudo_amino_acid_composition", "conjoint_triad",
            "ctd", "ctd_", "ctd_composition", "ctd_transition", "ctd_distribution",
            "sequence_order_coupling_number", "quasi_sequence_order"]
        for name in registered:
            self.assertIn(name, protpy.descriptor_names, f"{name} should be a registered descriptor.")
        #deprecated in v1.6.0: still callable (with a DeprecationWarning) but no longer registered
        deprecated = ["hydrophobic_polar_charged_composition", "secondary_structure_propensity",
            "reduced_alphabet_composition", "amino_acid_pair_composition", "aggregation_propensity",
            "sequence_order_coupling_number_all", "quasi_sequence_order_all", "ncpr", "charge_density",
            "hydrophobic_ratio", "half_life", "ivywrel_index", "cvp_bias", "flexibility", "linguistic_complexity",
            "fft_descriptors", "dct_descriptors", "uversky_charge_hydropathy", "fold_index", "disorder_propensity",
            "transmembrane_features", "signal_peptide_features", "functional_site_counts", "antigenicity_surface",
            "cysteine_features", "z_scales_5", "vhse", "atchley_factors", "t_scales", "st_scales", "protfp",
            "fasgai", "blosum_indices", "ms_whim", "physicochemical_distance_transformation"]
        seq = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDAKLVVSTEQW"
        for name in deprecated:
            with self.subTest(name=name):
                self.assertNotIn(name, protpy.descriptor_names)
                with self.assertWarns(DeprecationWarning):
                    getattr(protpy, name)(seq)
        #helpers no longer used by protpy are deprecated (and not exported)
        from protpy import composition, sequence_order
        for func, args in ((composition.sequence_order_correlation_factor, (seq, 1, [composition.hydrophobicity_])),
                           (composition.amphiphilic_sequence_order_correlation_factor, (seq, 1)),
                           (sequence_order.sequence_order_coupling_number_, (seq, 1))):
            with self.assertWarns(DeprecationWarning):
                func(*args)
        #internal helpers are no longer exported from the protpy namespace
        for name in ("sequence_order_coupling_number_", "str_to_num", "normalize_property",
                     "sequence_order_correlation_factor", "amphiphilic_sequence_order_correlation_factor"):
            self.assertFalse(hasattr(protpy, name), name)
        #permanent literature aliases give the same output as their family function
        self.assertTrue(protpy.tripeptide_composition(seq).equals(protpy.kmer_composition(seq, k=3)))
        self.assertTrue(protpy.dde(seq).equals(protpy.dipeptide_composition(seq, normalize="dde")))
        self.assertTrue(protpy.cksaagp(seq).equals(protpy.cksaap(seq, alphabet="ifeature5")))
        self.assertTrue(protpy.geary_autocorrelation(seq).equals(protpy.autocorrelation(seq, method="geary")))
        self.assertTrue(protpy.ctd_(seq).equals(protpy.ctd(seq)))
        self.assertTrue(protpy.ctd_composition(seq).equals(protpy.ctd(seq, all_ctd=False, parts="C")))
        self.assertAlmostEqual(protpy.boman_index(seq).iloc[0, 0], -protpy.gravy(seq, scale="RADA880101").iloc[0, 0])
