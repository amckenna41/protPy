# Changelog

All notable changes to **protpy** are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [Unreleased]

### Added
- `descriptors.pdf` — PDF version of `DESCRIPTORS.md`, providing a printable/downloadable reference for all 36 descriptors.

### Changed
- `aaindex` dependency version floor raised to `>=1.2.0` to pick up the new `aaindex2` and `aaindex3` modules, `pyproject.toml` migration, and all associated bug fixes in that release.

---

## [1.4.0] - 2026-04-04

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
- Version bumped to `1.4.0` across `pyproject.toml`, `protpy/__init__.py`, `docs/conf.py`, and `tests/test_protpy.py`.

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

[Unreleased]: https://github.com/amckenna41/protPy/compare/v1.3.0...HEAD
[1.3.0]: https://github.com/amckenna41/protPy/compare/v1.2.1...v1.3.0
[1.2.1]: https://github.com/amckenna41/protPy/compare/v1.2.0...v1.2.1
[1.2.0]: https://github.com/amckenna41/protPy/compare/v1.1.0...v1.2.0
[1.1.0]: https://github.com/amckenna41/protPy/compare/v1.0.0...v1.1.0
[1.0.0]: https://github.com/amckenna41/protPy/releases/tag/v1.0.0
