# protpy Unit Tests

All tests in the project were ran using Python's unittest testing framework (https://docs.python.org/3/library/unittest.html).

**Run all unittests from main project directory:**
```
python3 -m unittest discover tests
```

You can add the flag *-b* to suppress some of the verbose output when running the unittests. The CI workflow runs the same tests with `python -m pytest`.

Unit tests
----------
* `test_protpy.py` - unit tests for overall `protpy` package (version, exported descriptor names).
* `test_autocorrelation.py` - unit tests for autocorrelation module and functionality.
* `test_composition.py` - unit tests for composition module and functionality, including pseudo composition and amphiphillic composition.
* `test_composition_extended.py` - v1.5.0 composition descriptors (CKSAAP, CKSAAGP, DDE, ASDC, distance pair, natural vector, atomic composition) and the GAAC/GDPC/GTPC, KSCTriad, 13-property CTD and split-composition extensions.
* `test_conjoint_triad.py` - unit tests for conjoint triad module and functionality.
* `test_ctd.py` - unit tests for CTD module and functionality.
* `test_sequence_order.py` - unit tests for sequence order module and functionality.
* `test_scales.py` - factor scales, auto/cross covariance and PDT (optional cross-check against `peptides`, skipped if not installed).
* `test_physicochemical.py` - charge patterning (including the published Das & Pappu kappa values), pKa scales, hydrophobicity profile, half-life, flexibility and the GRAVY / hydrophobic moment / extinction coefficient extensions.
* `test_complexity.py` - Shannon entropy, Wootton-Federhen, linguistic (deprecated) and Lempel-Ziv complexity.
* `test_spectral.py` - `spectral_descriptors` (FFT and DCT) and the deprecated `fft_descriptors` / `dct_descriptors` aliases.
* `test_structural.py` - AGGRESCAN, Uversky/FoldIndex, disorder propensity, TM and signal-peptide features, functional sites, PLAAC (against PLAAC's reference output), antigenicity and cysteine features.
* `test_encode.py` - `all_descriptors`, `encode_dataset` and `read_sequences_from_file`.
* `test_descriptor_sets.py` - `descriptor_set`, `include`, `exclude` and `decimals`, the 1,402-column classic set, short sequences (NaN-padded lags, same columns), `params`/`include` checks, and every registered descriptor on 1-8 residue sequences.
* `test_utils.py` - descriptor registry and metadata functions.
* `test_helpers.py` - internal helpers, data-table provenance, registry parameters matching function signatures, and registry DOIs present in `docs/references.bib`.
* `test_regressions.py` - regression tests for bugs fixed in the post-1.4.2 audit (updated for the v1.6.0 lag policy) and in the v1.6.0 final review (numeric validation, motif regex errors, LZ76 against Kaspar-Schuster, FASTA identifiers and stop codons, CSV escaping, scale-table copies).
* `test_backwards_compat.py` - golden snapshot: the default output of every registered descriptor and of the three `all_descriptors` sets must equal the v1.6.0 baseline (re-baselined from v1.4.2 in v1.6.0, whose audit fixes change values on purpose). Regenerate only for a deliberate, documented change: `python tests/test_backwards_compat.py --regen`.
* `test_reference_values.py` - checks against independent references rather than earlier protpy outputs: ProtParam values against Biopython, composition sums (QSO = 1, PseAAC = 100), the identities found by the v1.6.0 audit (Moreau-Broto = auto covariance, PDT = 2 s² Geary, DPC = k-mer k=2 = CKSAAP gap 0), ordered CTD quartiles and published kappa values.
* `_compare.py` - tolerance helper for expected values written from earlier, rounded outputs.
* `test_performance.py` - timing test for `descriptor_set="all"` on a 500-residue sequence; skipped unless `PROTPY_PERF=1` is set.

Test Files
----------
* `test_data/test_fasta1.fasta` - Spike glycoprotein (P59594 · SPIKE_SARS) [[1]](#references).
* `test_data/test_fasta2.fasta` - Nucleoprotein (P59595 · NCAP_SARS) [[2]](#references).
* `test_data/test_fasta3.fasta` - Envelope small membrane protein (P59637 · VEMP_SARS) [[3]](#references).
* `test_data/test_fasta4.fasta` - Nucleoprotein (P0DTC9 · NCAP_SARS2) [[4]](#references).
* `test_data/structural_reference.fasta` - MOT3 (from the MIT-licensed PLAAC repository example), glycophorin A (P02724) [[5]](#references) and amyloid-β 1-42 (P05067 residues 672-713) [[6]](#references), used by `test_structural.py`.
* `test_data/golden_v1_6_0.json.gz` - golden outputs (every registered descriptor and the three `all_descriptors` sets, for the four test FASTA sequences) used by `test_backwards_compat.py`.

References
----------
\[1\]: https://www.uniprot.org/uniprotkb/P59594/entry <br>
\[2\]: https://www.uniprot.org/uniprotkb/P59595/entry <br>
\[3\]: https://www.uniprot.org/uniprotkb/P59637/entry <br>
\[4\]: https://www.uniprot.org/uniprotkb/P0DTC9/entry <br>
\[5\]: https://www.uniprot.org/uniprotkb/P02724/entry <br>
\[6\]: https://www.uniprot.org/uniprotkb/P05067/entry 