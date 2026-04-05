###########################################################################
###############       protPy - Utils Module Tests           ###############
###########################################################################

import io
import sys
import unittest
unittest.TestLoader.sortTestMethodsUsing = None

import protpy
from protpy.utils import (
    _DESCRIPTOR_REGISTRY,
    _SPEED_SYMBOL,
    get_descriptor_info,
    list_descriptors,
    print_descriptor_info,
)


class ProtpyUtilsTests(unittest.TestCase):
    """
    Test suite for testing the utils module and its three public functions:
    get_descriptor_info, list_descriptors, and print_descriptor_info.

    Test Cases
    ----------
    test_registry_completeness:
        All 36 descriptor names from all_descriptors appear in _DESCRIPTOR_REGISTRY.
    test_registry_schema:
        Every registry entry contains the required metadata keys with correct types.
    test_get_descriptor_info_valid:
        Known descriptor returns a dict with the correct "name" key and expected fields.
    test_get_descriptor_info_case_insensitive:
        get_descriptor_info is case-insensitive and strips surrounding whitespace.
    test_get_descriptor_info_returns_copy:
        The returned dict is a copy; mutating it does not affect the registry.
    test_get_descriptor_info_unknown_raises:
        An unknown name raises KeyError.
    test_get_descriptor_info_suggestion_hint:
        A near-miss name includes a close-match suggestion in the KeyError message.
    test_list_descriptors_all:
        list_descriptors() with no argument returns a sorted list of all 36 descriptors.
    test_list_descriptors_by_category:
        Filtering by each category returns a non-empty subset with correct category values.
    test_list_descriptors_case_insensitive:
        Category filter is case-insensitive and strips whitespace.
    test_list_descriptors_unknown_category:
        An unrecognised category returns an empty list (not an error).
    test_print_descriptor_info_output:
        print_descriptor_info writes expected fields to stdout.
    test_print_descriptor_info_unknown_raises:
        print_descriptor_info propagates the KeyError from get_descriptor_info.
    test_speed_symbol_mapping:
        _SPEED_SYMBOL covers Fast, Moderate, and Slow.
    """

    # ── Shared known-good descriptors for quick look-ups ─────────────────────
    KNOWN_NAMES = [
        "amino_acid_composition",
        "dipeptide_composition",
        "tripeptide_composition",
        "moreaubroto_autocorrelation",
        "moran_autocorrelation",
        "geary_autocorrelation",
        "conjoint_triad",
        "ctd_composition",
        "ctd_transition",
        "ctd_distribution",
        "ctd_",
        "sequence_order_coupling_number",
        "quasi_sequence_order",
        "pseudo_amino_acid_composition",
        "amphiphilic_pseudo_amino_acid_composition",
    ]

    # Required keys that every registry entry must contain
    REQUIRED_KEYS = {
        "description", "module", "category",
        "abbreviation", "parameters", "output_shape", "speed", "reference",
    }

    VALID_CATEGORIES = {"Composition", "Autocorrelation", "Conjoint Triad", "CTD", "Sequence Order"}
    VALID_SPEEDS = {"Fast", "Moderate", "Slow"}

# ── Registry integrity ────────────────────────────────────────────────────────

    def test_registry_completeness(self):
        """All 36 descriptor names exported by protpy appear in _DESCRIPTOR_REGISTRY."""
        missing = [n for n in protpy.all_descriptors if n not in _DESCRIPTOR_REGISTRY]
        self.assertEqual(
            missing, [],
            f"Descriptors in all_descriptors but missing from registry: {missing}",
        )
        self.assertEqual(len(_DESCRIPTOR_REGISTRY), len(protpy.all_descriptors),
            "Registry length doesn't match len(all_descriptors).")

    def test_registry_schema(self):
        """Every registry entry has all required keys with the correct types."""
        for name, entry in _DESCRIPTOR_REGISTRY.items():
            with self.subTest(descriptor=name):
                # All required keys present
                missing_keys = self.REQUIRED_KEYS - entry.keys()
                self.assertFalse(missing_keys,
                    f"Entry '{name}' is missing keys: {missing_keys}")

                # Value type checks
                self.assertIsInstance(entry["description"],   str,  f"'{name}' description must be str")
                self.assertIsInstance(entry["module"],        str,  f"'{name}' module must be str")
                self.assertIsInstance(entry["category"],      str,  f"'{name}' category must be str")
                self.assertIsInstance(entry["abbreviation"],  str,  f"'{name}' abbreviation must be str")
                self.assertIsInstance(entry["parameters"],    dict, f"'{name}' parameters must be dict")
                self.assertIsInstance(entry["output_shape"],  str,  f"'{name}' output_shape must be str")
                self.assertIsInstance(entry["speed"],         str,  f"'{name}' speed must be str")
                self.assertIsInstance(entry["reference"],     str,  f"'{name}' reference must be str")

                # Category must be an exact match from the controlled vocabulary
                self.assertIn(entry["category"], self.VALID_CATEGORIES,
                    f"'{name}' has unrecognised category '{entry['category']}'")
                # Speed must contain at least one base speed word (allows compound values
                # like "Fast (k≤3) / Slow (k=4)" used for descriptors with variable cost)
                self.assertTrue(
                    any(s in entry["speed"] for s in self.VALID_SPEEDS),
                    f"'{name}' has unrecognised speed '{entry['speed']}'"
                )

# ── get_descriptor_info ───────────────────────────────────────────────────────

    def test_get_descriptor_info_valid(self):
        """get_descriptor_info returns a complete dict for every known descriptor."""
        for name in self.KNOWN_NAMES:
            with self.subTest(descriptor=name):
                info = get_descriptor_info(name)
                self.assertIsInstance(info, dict)
                # "name" key is injected by get_descriptor_info
                self.assertEqual(info["name"], name)
                # All registry keys are present in the returned dict
                for key in self.REQUIRED_KEYS:
                    self.assertIn(key, info, f"Key '{key}' missing from get_descriptor_info('{name}')")

    def test_get_descriptor_info_specific_values(self):
        """Spot-check specific metadata values for amino_acid_composition."""
        info = get_descriptor_info("amino_acid_composition")
        self.assertEqual(info["category"],     "Composition")
        self.assertEqual(info["module"],       "composition")
        self.assertEqual(info["abbreviation"], "AAComp")
        self.assertEqual(info["output_shape"], "1 x 20")
        self.assertEqual(info["speed"],        "Fast")

    def test_get_descriptor_info_case_insensitive(self):
        """Descriptor name lookup ignores case and strips whitespace."""
        variants = [
            "Amino_Acid_Composition",
            "AMINO_ACID_COMPOSITION",
            "  amino_acid_composition  ",
            "Amino_acid_COMPOSITION",
        ]
        for variant in variants:
            with self.subTest(variant=variant):
                info = get_descriptor_info(variant)
                self.assertEqual(info["name"], "amino_acid_composition")

    def test_get_descriptor_info_returns_copy(self):
        """Mutating the returned dict must not affect the registry."""
        info = get_descriptor_info("gravy")
        original_desc = _DESCRIPTOR_REGISTRY["gravy"]["description"]
        info["description"] = "MUTATED"
        self.assertEqual(
            _DESCRIPTOR_REGISTRY["gravy"]["description"],
            original_desc,
            "Registry was mutated through the returned dict.",
        )

    def test_get_descriptor_info_unknown_raises(self):
        """An unrecognised name raises KeyError."""
        with self.assertRaises(KeyError):
            get_descriptor_info("not_a_real_descriptor")

    def test_get_descriptor_info_suggestion_hint(self):
        """A near-miss name includes a suggestion in the KeyError message."""
        with self.assertRaises(KeyError) as ctx:
            get_descriptor_info("amino_acid_compositon")  # deliberate typo
        self.assertIn("amino_acid_composition", str(ctx.exception),
            "Expected close-match suggestion to appear in error message.")

    def test_get_descriptor_info_empty_string_raises(self):
        """An empty string raises KeyError."""
        with self.assertRaises(KeyError):
            get_descriptor_info("")

    def test_get_descriptor_info_non_string_raises(self):
        """Passing a non-string raises AttributeError (no .strip() on int/None)."""
        for bad_input in [123, None, [], True]:
            with self.subTest(value=bad_input):
                with self.assertRaises((AttributeError, KeyError)):
                    get_descriptor_info(bad_input)

# ── list_descriptors ──────────────────────────────────────────────────────────

    def test_list_descriptors_all(self):
        """No-arg call returns a sorted list containing all 36 descriptors."""
        descriptors = list_descriptors()
        self.assertIsInstance(descriptors, list)
        self.assertEqual(len(descriptors), len(_DESCRIPTOR_REGISTRY))
        # Result is sorted
        self.assertEqual(descriptors, sorted(descriptors))
        # Every name in all_descriptors appears in the result
        for name in protpy.all_descriptors:
            self.assertIn(name, descriptors)

    def test_list_descriptors_by_category(self):
        """Each valid category returns a non-empty, sorted list of names."""
        for cat in self.VALID_CATEGORIES:
            with self.subTest(category=cat):
                result = list_descriptors(cat)
                self.assertIsInstance(result, list)
                self.assertGreater(len(result), 0, f"Expected at least one descriptor in '{cat}'")
                self.assertEqual(result, sorted(result), "Result should be sorted.")
                # Every returned name actually belongs to this category
                for name in result:
                    self.assertEqual(
                        _DESCRIPTOR_REGISTRY[name]["category"], cat,
                        f"'{name}' returned for category '{cat}' but has "
                        f"category '{_DESCRIPTOR_REGISTRY[name]['category']}'",
                    )

    def test_list_descriptors_case_insensitive(self):
        """Category filter ignores case and leading/trailing whitespace."""
        expected = list_descriptors("Composition")
        for variant in ["composition", "COMPOSITION", "  Composition  "]:
            with self.subTest(variant=variant):
                self.assertEqual(list_descriptors(variant), expected)

    def test_list_descriptors_no_overlap_between_categories(self):
        """The union of all per-category lists equals the full list, with no duplicates."""
        full = list_descriptors()
        combined = []
        for cat in self.VALID_CATEGORIES:
            combined.extend(list_descriptors(cat))
        self.assertEqual(sorted(combined), full)

    def test_list_descriptors_unknown_category(self):
        """An unrecognised category returns an empty list without raising."""
        result = list_descriptors("NotACategory")
        self.assertEqual(result, [])

    def test_list_descriptors_none_same_as_no_arg(self):
        """Explicitly passing None is identical to calling with no argument."""
        self.assertEqual(list_descriptors(None), list_descriptors())

# ── print_descriptor_info ─────────────────────────────────────────────────────

    def test_print_descriptor_info_output(self):
        """print_descriptor_info writes name, category, speed, and output shape to stdout."""
        captured = io.StringIO()
        sys.stdout = captured
        try:
            print_descriptor_info("gravy")
        finally:
            sys.stdout = sys.__stdout__

        output = captured.getvalue()
        self.assertIn("gravy",       output)
        self.assertIn("Composition", output)
        self.assertIn("Fast",        output)
        self.assertIn("1 x 1",       output)

    def test_print_descriptor_info_contains_parameters(self):
        """print_descriptor_info includes parameter descriptions in its output."""
        captured = io.StringIO()
        sys.stdout = captured
        try:
            print_descriptor_info("moreaubroto_autocorrelation")
        finally:
            sys.stdout = sys.__stdout__

        output = captured.getvalue()
        # "Parameters:" section header must appear
        self.assertIn("Parameters", output)
        # "sequence" parameter must be listed
        self.assertIn("sequence", output)

    def test_print_descriptor_info_contains_reference(self):
        """print_descriptor_info includes a reference string in its output."""
        captured = io.StringIO()
        sys.stdout = captured
        try:
            print_descriptor_info("conjoint_triad")
        finally:
            sys.stdout = sys.__stdout__

        self.assertIn("Reference", captured.getvalue())

    def test_print_descriptor_info_unknown_raises(self):
        """print_descriptor_info propagates KeyError for unknown names."""
        with self.assertRaises(KeyError):
            print_descriptor_info("not_a_real_descriptor")

    def test_print_descriptor_info_case_insensitive(self):
        """print_descriptor_info accepts mixed-case and whitespace-padded names."""
        captured = io.StringIO()
        sys.stdout = captured
        try:
            print_descriptor_info("  GRAVY  ")
        finally:
            sys.stdout = sys.__stdout__

        self.assertIn("gravy", captured.getvalue())

# ── _SPEED_SYMBOL internal constant ──────────────────────────────────────────

    def test_speed_symbol_mapping(self):
        """_SPEED_SYMBOL covers all valid speed values."""
        for speed in self.VALID_SPEEDS:
            self.assertIn(speed, _SPEED_SYMBOL,
                f"_SPEED_SYMBOL is missing an entry for speed '{speed}'")
            self.assertIsInstance(_SPEED_SYMBOL[speed], str)


if __name__ == "__main__":
    unittest.main()
