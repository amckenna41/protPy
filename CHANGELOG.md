# Changelog

All notable changes to **protpy** are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [1.5.0] - 2026-10-12

### Added
- **44 new registered descriptor functions** (plus the generic `aa_scale_descriptors`) with literature references (DOIs in the registry and `docs/references.bib`):
  - Composition: `cksaap`, `cksaagp`, `dde`, `asdc`, `distance_pair` (cp13/cp14/cp19/cp20), `natural_vector`, `atomic_composition`.
  - Scales (`protpy/scales.py`): `aa_scale_descriptors` with wrappers `z_scales_5`, `vhse`, `atchley_factors`, `t_scales`, `st_scales`, `protfp`, `fasgai`, `blosum_indices`, `ms_whim` (all with `domains=`); `auto_cross_covariance`; `physicochemical_distance_transformation`.
  - Physicochemical (`protpy/physicochemical.py`): `charge_fractions`, `ncpr`, `kappa` (reproduces the published Das & Pappu values), `sequence_charge_decoration`, `charge_density`, `hydrophobicity_profile`, `hydrophobic_ratio`, `half_life`, `ivywrel_index`, `cvp_bias`, `flexibility`.
  - Complexity (`protpy/complexity.py`): `wootton_federhen_complexity`, `linguistic_complexity`, `lempel_ziv_complexity` (`shannon_entropy` re-exported).
  - Signal (`protpy/signal.py`, numpy only): `fft_descriptors`, `dct_descriptors`.
  - Structural (`protpy/structural.py`): `aggrescan`, `uversky_charge_hydropathy`, `fold_index`, `disorder_propensity`, `transmembrane_features`, `signal_peptide_features`, `functional_site_counts`, `plaac` (port of the MIT-licensed PLAAC; reproduces its reference output), `antigenicity_surface`, `cysteine_features`.
- **Extensions to existing functions** (defaults unchanged): `reduced_alphabet_composition(alphabet=, k=)` (GAAC/GDPC/GTPC), `conjoint_triad(gap=)` (KSCTriad), `ctd_(property_set="ifeature13")`, `amino_acid_composition`/`dipeptide_composition(region=, n=)` (split composition), `isoelectric_point`/`charge_distribution(pka_scale=)` ("emboss", "bjellqvist", "ipc_protein"), `gravy(scale=)`, `hydrophobic_moment(scale=, angle="helix"/"sheet")`, `extinction_coefficient(absorbance=True)`.
- **`descriptor_set`, `include`, `exclude`** for `all_descriptors` and `encode_dataset`: `"classic"` (default, unchanged), `"extended"` (adds the new descriptors with ≤150 columns) and `"all"`. Set membership comes from the registry's new `descriptor_sets` field. Every parameterised v1.5.0 descriptor is accepted in `params`.
- **Data files** with provenance (`name`, `values`, `reference`, `doi`, `source_note`): `protpy/data/aa_scales.json`, `reduced_alphabets.json`, `residue_tables.json`, `structural_params.json`.
- Internal helpers in `protpy/_helpers.py`. The registry gains `doi`, `in_default_set` and `descriptor_sets`, and `print_descriptor_info` prints the DOI. New registry categories: Scales, Physicochemical, Complexity, Signal, Structural.
- `docs/references.bib` with BibTeX for every registry DOI.
- Tests: golden backwards-compatibility snapshot of every v1.4.2 default output (`tests/test_backwards_compat.py`), helper/registry tests, one test module per new module, descriptor-set tests, and an opt-in performance test (`PROTPY_PERF=1`).
- Plan items that could not be implemented without unverifiable source data (e.g. PseKRAAC, Kidera factors, GOR, COILS, Chou–Fasman) are listed with reasons in `planning/BLOCKED.md`.

### Changed
- `encode_dataset(verbose=False)` now also suppresses the descriptor-list line printed for the first sequence.
- `hydrophobic_moment` docstring now names the source of its default scale (Eisenberg, Schwarz, Komarony & Wall 1984).

### Fixed (post-1.4.2 audit)
- **Autocorrelation crash on fresh installs**: `aaindex1[code]` can return a new record on each access, so `.values.pop('-')` was lost and all three autocorrelation descriptors (and `all_descriptors`/`encode_dataset`) raised `IndexError`. Property values are now copied and filtered to the 20 canonical amino acids; aaindex records are no longer mutated.
- **Moran autocorrelation formula**: multiplied by the property variance instead of dividing by it.
- **CTD groups**: `polarity` groups 2/3 were copied from polarizability (A, D, G, S, T unassigned; 8 residues duplicated) — now `LIFWCMVY / PAGST / HQRKNED`. `normalized_vdwv` group 1 was missing C.
- **Overlapping counts**: dipeptide, tripeptide and conjoint triad used `str.count`, which skips overlapping matches (e.g. `"AA"` in `"AAA"`).
- **Short/repetitive sequences**: APAAComp raised `ZeroDivisionError` for sequences <= 30 residues; dipeptide composition, CTD transition, Moran and Geary crashed or returned NaN for single-residue or homopolymer sequences; lags with no residue pairs now return 0.0.
- **`all_descriptors(params={"ctd_": {"all_ctd": False}})`** raised `TypeError`.
- **Motif `dileucine_sorting`** never matched (literal `xxx` and stray `$`).
- **Boman index** now uses the published Radzicka-Wolfenden solubility scale and `BI = -sum/N`.
- **Isoelectric point** uses bisection over pH 0-14 (previously capped at 12.0 for very basic proteins).
- **PseAAC with AAIndex properties** included the `'-'` entry in normalization.
- **Glycine molecular weight** 75.032 -> 75.067.
- **`aggregation_propensity`**: `charge_threshold` is now an inclusive maximum, as documented; default changed 1 -> 0 so default results are unchanged.
- **`hydrophobic_moment`** registry entry now reports its `1 x 2` output shape.
- **MoreauBroto Autocorrelation (Bug 1)**: Two intertwined bugs in `moreaubroto_autocorrelation` — the lag index always read `sequence[j+1]` (hardcoded lag=1) regardless of the current lag `i`, and the accumulator `temp` was not reset between lags, causing values to accumulate across all lags. Both corrected: `sequence[j+1]` → `sequence[j+i]`, and `temp = 0` moved inside the lag loop.
- **PseAAC normalization denominator (Bug 2)**: `pseudo_amino_acid_composition` used `1 + weight + Σθ_k` where the Chou (2001) formula requires `1 + weight * Σθ_k`. Fixed operator `+` → `*`. Removed the redundant first normalization loop and unified output computation into a single pass.
- **APAAComp properties parameter ignored (Bug 3)**: `amphiphilic_pseudo_amino_acid_composition` accepted a `properties` parameter but the underlying helpers `amphiphilic_correllation_function` and `amphiphilic_sequence_order_correlation_factor` hardcoded `hydrophobicity_` and `hydrophilicity_`, silently ignoring any user-supplied properties. Fixed by threading `properties` through all three functions and pre-normalizing once per call.
- **`JSONDecodeError` constructor misuse (Bug 5)**: `sequence_order_coupling_number_` raised `json.JSONDecodeError(msg)` which fails at runtime with `TypeError` because `JSONDecodeError` requires `(msg, doc, pos)` positional arguments. Changed to `raise ValueError(...)`.

### Performance (post-1.4.2 audit)
- `kmer_composition`, dipeptide/tripeptide/pair composition count k-mers in one pass; QSO computes each SOCN once; PseAAC normalizes properties once per gap instead of per residue pair.
- **`sequence_order_coupling_number_` disk I/O (Issue 6)**: Added `@functools.lru_cache(maxsize=2)` loader `_load_distance_matrix` in `sequence_order.py`. The distance matrix JSON was read from disk on every call; for a sequence with `lag=30`, that was 30 file reads per descriptor invocation. Now both matrices are cached in memory after the first load.
- **`normalize_property` O(n²) (Issue 7)**: Refactored to compute `mean` and `std` once before the loop, reducing from O(n²) to O(n) for the 20-amino-acid property dictionaries used in PseAAC and APAAComp descriptors.

### Changed (post-1.4.2 audit)
- Invalid `lag`, `lamda`, `k`, `window` (non-integer or negative) and `weight` (non-numeric or negative) now raise `ValueError` instead of being silently replaced by defaults. `lag`/`lamda` larger than the sequence still fall back to 30.
- `all_descriptors` no longer writes `descriptors.csv` by default (`output_file=None`).
- `all_descriptors` no longer includes `amino_acid_pair_composition` (it duplicated dipeptide composition's 400 values).
- `encode_dataset`/`all_descriptors` refuse to overwrite an existing output file with different columns, append rows instead of rewriting the whole file per sequence, and `encode_dataset` replaces rows for re-encoded IDs instead of duplicating them.
- **`tripeptide_composition` returns percentages (Gap 8)**: Output values are now fractional percentages (`count / (N-2) * 100`, rounded to 3 dp) consistent with `dipeptide_composition`, instead of raw integer counts. Column dtype is now `float64` instead of `int64`.
- **`encode_dataset` verbose parameter (Gap 11)**: Added `verbose: bool = True` parameter to `encode_dataset`. When `False`, all per-sequence stdout status lines are suppressed. Default `True` preserves existing behaviour.

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
