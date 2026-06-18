# Changelog

All notable changes to **protpy** are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [Unreleased]

---

## [1.4.1] - 2026-06-18

### Added
- `encode_dataset(sequences, seq_col, params, show_progress, output_file, save_each_sequence, strip_gaps)` — new public function in `protpy/encode.py` that computes all descriptors for every sequence in a dataset (plain `list[str]` or `pd.DataFrame`) and returns a single DataFrame whose first column is a unique sequence identifier. Auto-generates `seq1, seq2, …` identifiers for list input; uses the first column of the DataFrame as the identifier when a DataFrame is supplied. Each completed row is flushed to the output CSV immediately so that interrupted runs preserve progress. Raises `ValueError` for gap characters unless `strip_gaps=True`.
- `read_sequences_from_file(filepath, seq_col, id_col, has_header)` — new public function in `protpy/encode.py` for parsing protein sequences from FASTA (`.fasta` / `.fa`), plain TXT, CSV, and TSV files. Format is inferred from the file extension. Returns `list[str]` when `id_col` is omitted; returns a two-column `pd.DataFrame` (identifier + sequence) when `id_col` is provided. Handles multi-line FASTA records, blank-line skipping, uppercase normalisation, headerless files, column selection by name or zero-based integer index, and auto-detection of the delimiter (`,` vs `\t`) for TXT files with `seq_col`.
- `all_descriptors` new parameters:
  - `show_progress: bool` (default `False`) — displays a live `tqdm` progress bar naming the descriptor currently being computed and tracking progress across all steps.
  - `output_file: str | None` (default `"descriptors.csv"`) — path to a CSV file where results are written; pass `None` to disable all file output.
  - `save_each_sequence: bool` (default `True`) — flushes the output CSV after each sequence completes, preserving progress across interrupted multi-sequence loops.
  - `_show_descriptor_list: bool` (internal, default `True`) — suppresses the descriptor name print when called from `encode_dataset` for sequences after the first.
- `encode_dataset` and `read_sequences_from_file` exported from the top-level `protpy` namespace.
- `tqdm>=4.0.0` added as a runtime dependency in `pyproject.toml` and `requirements.txt`.
- `tests/test_encode.py` — extended with three new test classes:
  - `ProtpyEncodeTests` (existing): two new test cases added — `test_all_descriptors_show_progress` (verifies `tqdm` is instantiated with correct arguments and all bar lifecycle methods are called) and `test_all_descriptors_no_progress_by_default` (verifies `tqdm` is never called when `show_progress` defaults to `False` or is set explicitly).
  - `ProtpyReadSequencesTests` (new) — 18 test cases covering all supported file formats (FASTA single/multi-record, `.fa` extension, plain TXT, multi-column TXT, CSV by name and index, TSV, headerless files), uppercase normalisation, blank-entry exclusion, `FileNotFoundError`, unsupported extension, bad column name/index/type errors.
  - `ProtpyEncodeDatasetTests` (new) — 11 test cases covering list and DataFrame inputs, auto-generated vs named identifiers, first-column guarantee, empty inputs, gap detection and stripping, bad `seq_col`, bad input type, and CSV output file correctness.
- `DESCRIPTORS.md` and `descriptors.pdf` — new "Utility Functions" section documenting `all_descriptors` (updated parameters), `encode_dataset`, and `read_sequences_from_file` with parameter tables and code examples.

### Changed
- Version bumped to `1.4.1` across `pyproject.toml`, `protpy/__init__.py`, `docs/conf.py`, and `tests/test_protpy.py`.

---

## [1.4.0] - 2026-06-14

### Added
- `protpy/encode.py` — new module with a single public function `all_descriptors(sequence, params=None)` that computes all 30 descriptor groups and returns a single concatenated `DataFrame` of shape `1 × 10,757` with default parameters. Combined variants are used to avoid redundancy: `ctd_` with `all_ctd=True` covers all 7 CTD properties; `sequence_order_coupling_number_all` and `quasi_sequence_order_all` cover both physicochemical distance matrices.
  - Accepts an optional `params: dict[str, dict]` argument for per-descriptor keyword overrides. Keys are descriptor function names; values are dicts of kwargs forwarded directly to that descriptor. Descriptors absent from `params` use their defaults unchanged.
  - Validates `params` upfront: unknown or non-parameterisable descriptor names raise `ValueError` with the full list of valid keys; a non-dict `params` raises `TypeError`.
  - Parameterisable descriptors: `charge_distribution` (`ph`), `kmer_composition` (`k`), `reduced_alphabet_composition` (`alphabet_size`), `motif_composition` (`motifs`), `aggregation_propensity` (`window`, `hydrophobicity_threshold`, `charge_threshold`), `hydrophobic_moment` (`window`, `angle`), `pseudo_amino_acid_composition` (`lamda`, `weight`, `properties`), `amphiphilic_pseudo_amino_acid_composition` (`lamda`, `weight`, `properties`), `moreaubroto_autocorrelation` (`lag`, `properties`, `normalize`), `moran_autocorrelation` (`lag`, `properties`, `normalize`), `geary_autocorrelation` (`lag`, `properties`, `normalize`), `ctd_` (`all_ctd`), `sequence_order_coupling_number_all` (`lag`), `quasi_sequence_order_all` (`lag`, `weight`).
- `tests/test_encode.py` — new test module with four test cases: `test_all_descriptors_default` (output type, shape, column names, spot-check values for all four test sequences); `test_all_descriptors_custom_params` (shape and value changes for `k`, `ph`, `lag`, `alphabet_size`, and combined multi-param calls); `test_all_descriptors_params_validation` (unknown keys, non-parameterisable descriptors, wrong `params` types, misspelled names); `test_all_descriptors_invalid_sequence` (non-canonical amino acids, empty string, non-string types, combined with custom params).
- `all_descriptors` callable exported from top-level `protpy` namespace.
- `descriptors.pdf` — PDF version of `DESCRIPTORS.md`, providing a printable/downloadable reference for all 36 descriptors.

### Changed
- `all_descriptors` name-list in `protpy/__init__.py` renamed to `descriptor_names` to free the name for the new `all_descriptors()` function. References updated in `tests/test_protpy.py` and `tests/test_utils.py`. Users who accessed `protpy.all_descriptors` as a list should switch to `protpy.descriptor_names` or `protpy.list_descriptors()`.
- `aaindex` dependency version floor raised to `>=1.2.0` to pick up the new `aaindex2` and `aaindex3` modules, `pyproject.toml` migration, and all associated bug fixes in that release.
- Version bumped to `1.4.0` across `pyproject.toml`, `protpy/__init__.py`, `docs/conf.py`, and `tests/test_protpy.py`.

---

## [1.3.1] - 2026-04-04

### Added
- `protpy/utils.py` — new utility module with three public functions:
  - `get_descriptor_info(name)` — returns a metadata dict for any registered descriptor (description, module, category, abbreviation, parameters, output shape, speed rating, reference). Case-insensitive; includes close-match suggestions on unknown names.
  - `list_descriptors(category=None)` — returns a sorted list of all 36 descriptor names, optionally filtered by category.
  - `print_descriptor_info(name)` — pretty-prints descriptor metadata to stdout with speed symbols (✅ / 🟡 / 🔴).
  - `_DESCRIPTOR_REGISTRY` — internal dict containing structured metadata for all 36 descriptors.
- `examples/protpy_tutorial.ipynb` — interactive Jupyter notebook tutorial covering all 13 descriptor families with visualisations, performance benchmarks, and a multi-panel dashboard.
- `tests/test_utils.py` — unit tests for all three public utils functions and the `_DESCRIPTOR_REGISTRY` integrity (completeness, schema, controlled-vocabulary checks).
- `get_descriptor_info`, `list_descriptors`, and `print_descriptor_info` exported from top-level `protpy` package.

### Changed
- **Type hints** — all five descriptor modules (`autocorrelation.py`, `composition.py`, `ctd.py`, `sequence_order.py`, `_constants.py`) upgraded from bare `list`/`dict` annotations to subscripted generics (`list[str]`, `list[float]`, `dict[str, float]`, `dict[str, str]`, `tuple[float, float]`). Compatible with Python 3.9+ via `from __future__ import annotations`.
- **F-strings** — all `.format()` calls across all modules replaced with f-string equivalents (20 call sites updated).
- Version bumped to `1.3.1` across `pyproject.toml`, `protpy/__init__.py`, `docs/conf.py`, and `tests/test_protpy.py`.

---

## [1.3.0] - 2026-03-31

### Added
- 22 new composition descriptors: `gravy`, `aromaticity`, `instability_index`, `isoelectric_point`, `molecular_weight`, `charge_distribution`, `hydrophobic_polar_charged_composition`, `secondary_structure_propensity`, `kmer_composition`, `reduced_alphabet_composition`, `motif_composition`, `amino_acid_pair_composition`, `aliphatic_index`, `extinction_coefficient`, `boman_index`, `aggregation_propensity`, `hydrophobic_moment`, `shannon_entropy`, `pseudo_amino_acid_composition` (moved from autocorrelation group), `amphiphilic_pseudo_amino_acid_composition` (moved), `sequence_order_coupling_number_` (single-gap variant), `sequence_order_coupling_number_all`, `quasi_sequence_order_all`.
- `DESCRIPTORS.md` — comprehensive descriptor reference covering all 35 descriptors, including parameter tables, code examples, output shapes, speed/complexity ratings (✅ Fast / 🟡 Moderate / 🔴 Slow), and a full References section with original citations.
- `pyproject.toml` — modern PEP 517/518 build configuration replacing `setup.py` and `setup.cfg`.
- `images/protpy.png` — project logo added to repository and displayed in `README.md`.
- `CHANGELOG.md` — this file.

### Changed
- `all_descriptors` list in `protpy/__init__.py` expanded from 14 to 36 entries.
- `README.md` — Introduction descriptor list updated to all 35 descriptors grouped by category; Usage section expanded with code examples for all 17 previously undocumented composition descriptors.
- `tests/test_protpy.py` — `test_valid_descriptors` updated to assert `len(all_descriptors) == 36` and check all 36 descriptor names individually, grouped by category.
- `pyproject.toml` — version floor constraints added to core dependencies (`numpy>=1.16.0`, `pandas>=1.1.0`, `varname>=0.12.0`); `biopython` moved from `install_requires` to `optional-dependencies[test]`.
- Python version classifiers updated to include 3.11, 3.12, and 3.13 in both `pyproject.toml` and `docs/conf.py`.
- `requirements.txt` — `biopython` removed (test-only dependency).
- `.circleci/config.yml` CircleCI coverage job now actually runs `coverage run -m pytest tests/` and `coverage report` (previously installed tools without executing them).
- `kmer_composition` — added `ValueError` guard for `k > 4` to prevent accidental generation of 3.2 million+ column DataFrames.
- `docs/conf.py` — `release` updated to `1.3.0`.
- Version bumped to `1.3.0` across `pyproject.toml`, `protpy/__init__.py`, `docs/conf.py`, and `tests/test_protpy.py`.

### Removed
- `setup.py` and `setup.cfg` — replaced by `pyproject.toml`.
- `.circleci/` directory and `config.yml` — CircleCI pipeline removed; GitHub Actions is the sole CI/CD provider.
- CircleCI badge removed from `README.md`.
- CircleCI TODO item removed from `TODO.md`.

### Fixed
- `setup.cfg` typo `install_requies` → `install_requires` (previously caused all cfg-declared dependencies to be silently ignored).

---

## [1.2.1] - 2023-12-13

### Added
- `.readthedocs.yml` — ReadTheDocs configuration for hosted documentation.
- Renamed distance matrix JSON files from `physiochemical` to `physicochemical` (corrected spelling).

### Changed
- `README.md` — minor copy and link updates.
- `protpy/__init__.py` — metadata and keyword cleanup.
- `autocorrelation.py`, `composition.py`, `ctd.py`, `sequence_order.py` — comment and docstring corrections.
- `deploy_pypi.yml` and `deploy_test_pypi.yml` — workflow trigger and token reference updates.
- `tests/test_protpy.py` — assertion message updates.

### Fixed
- `protpy/data/` — distance matrix filenames corrected from `physiochemical` to `physicochemical` in both JSON filenames and all internal references.

---

## [1.2.0] - 2023-11-14

### Added
- Full GitHub Actions CI/CD pipeline:
  - `build_test.yml` — test matrix across Python 3.11, 3.12, 3.13 with pytest-cov and Codecov upload.
  - `deploy_test_pypi.yml` — automated deployment to TestPyPI after successful test run.
  - `deploy_pypi.yml` — automated deployment to production PyPI after successful TestPyPI deployment.
- `.github/workflows/README.md` — workflow documentation.
- `MANIFEST.in` — includes `protpy/data/` JSON files in source distributions.
- `protpy/data/README.md` — documentation for the bundled physicochemical distance matrices.
- `tests/README.md` — documentation for the test suite.
- `docs/` — Sphinx documentation structure with `conf.py`, `index.rst`, `api.rst`, `usage.rst`, `contributing.rst`, `Makefile`, and `make.bat`.
- `ToDo.md` — project task tracking file.

### Changed
- `README.md` — major expansion with full usage examples for all descriptors, requirements, installation, directory structure, and reference list.
- `setup.py` / `setup.cfg` — full metadata, classifiers, and dependency declarations.
- `protpy/__init__.py` — `all_descriptors` list formalised (14 entries).
- All source modules — comprehensive docstrings, inline comments, and reference citations added.
- `tests/` — unit test suites expanded and restructured across `test_composition.py`, `test_autocorrelation.py`, `test_conjoint_triad.py`, `test_ctd.py`, `test_sequence_order.py`, and `test_protpy.py`.

### Fixed
- Multiple descriptor calculation bugs across `composition.py`, `ctd.py`, `sequence_order.py`, and `autocorrelation.py`.
- CTD column naming corrected (zero-padded singular columns, e.g. `_01`).
- `secondary_struct` property key name corrected (was `sec_struct`).
- SOCN/QSO column naming standardised: `SOCNUm` → `SOCN`, `QSOrder` → `QSO`; `SW`/`Grant` suffixes appended.
- Sequences uppercased on input; leading/trailing whitespace stripped.
- Invalid amino acid characters now raise a `ValueError` with informative message.
- PseAAC hydrophobicity, hydrophilicity, and residue mass values hard-coded for reproducibility.
- `lag` and `weight` parameter validation added to sequence order functions.

---

## [1.1.0] - 2023-04-12

### Added
- Expanded unit test coverage across all descriptor modules.
- Additional inline comments and docstring improvements throughout all source files.

### Fixed
- CTD descriptor calculation corrections.
- Conjoint triad descriptor bug fixes.
- Pseudo amino acid composition function corrections.
- Sequence order module stability improvements.

---

## [1.0.0] - 2023-02-08

### Added
- Initial package release.
- Core descriptor modules:
  - `composition.py` — `amino_acid_composition`, `dipeptide_composition`, `tripeptide_composition`, `pseudo_amino_acid_composition`, `amphiphilic_pseudo_amino_acid_composition`.
  - `autocorrelation.py` — `moreaubroto_autocorrelation`, `moran_autocorrelation`, `geary_autocorrelation`.
  - `conjoint_triad.py` — `conjoint_triad`.
  - `ctd.py` — `ctd_composition`, `ctd_transition`, `ctd_distribution`, `ctd_`.
  - `sequence_order.py` — `sequence_order_coupling_number`, `quasi_sequence_order`.
- `protpy/__init__.py` — package entry point with metadata and `all_descriptors` list.
- `protpy/data/` — bundled Grantham and Schneider-Wrede physicochemical distance matrices (JSON).
- `tests/` — initial unit test suite with test FASTA data files.
- `requirements.txt` — `aaindex`, `numpy`, `pandas`, `varname`.
- `setup.py` / `setup.cfg` — package build and distribution configuration.
- `README.md` — initial project documentation.
- `LICENSE` — MIT licence.
- `.gitignore` — Python project ignore rules.

---

[Unreleased]: https://github.com/amckenna41/protPy/compare/v1.4.1...HEAD
[1.4.1]: https://github.com/amckenna41/protPy/compare/v1.4.0...v1.4.1
[1.4.0]: https://github.com/amckenna41/protPy/compare/v1.3.1...v1.4.0
[1.3.1]: https://github.com/amckenna41/protPy/compare/v1.3.0...v1.3.1
[1.3.0]: https://github.com/amckenna41/protPy/compare/v1.2.1...v1.3.0
[1.2.1]: https://github.com/amckenna41/protPy/compare/v1.2.0...v1.2.1
[1.2.0]: https://github.com/amckenna41/protPy/compare/v1.1.0...v1.2.0
[1.1.0]: https://github.com/amckenna41/protPy/compare/v1.0.0...v1.1.0
[1.0.0]: https://github.com/amckenna41/protPy/releases/tag/v1.0.0
