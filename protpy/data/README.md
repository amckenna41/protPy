# Data files used by protpy descriptors

* `grantham-physicochemical-distance-matrix.json` - grantham distance score for each amino acid combination, used in the calculation of the sequence-order related protein descriptors [[1]](#references).
* `schneider-wrede-physicochemical-distance-matrix.json` - Schneider-Wrede distance score for each amino acid combination, used in the calculation of the sequence-order related protein descriptors [[2]](#references).

Added in v1.5.0. Each file maps a table name to `{"name", "values", "reference", "doi", "source_note"}`; `source_note` records where and when every table was transcribed.

* `aa_scales.json` - single-value hydrophobicity/property scales (Kyte-Doolittle, Hopp-Woods, two Eisenberg 1984 consensus scales — `eisenberg_jmb1984` (A = 0.62, the default `hydrophobic_moment` table) and `eisenberg_arb1984` (AAindex EISD840101, A = 0.25; called `eisenberg` before v1.6.0) — Wimley-White, Engelman, DomCut, EIIP, Vihinen flexibility, Parker, Janin, TOP-IDP) and the factor-scale families (z-scales 5, VHSE, Atchley, T-, ST-scales, ProtFP, FASGAI, BLOSUM indices, MS-WHIM). Tables with `"experimental": true` (z-scales 5, VHSE, ST-scales, FASGAI) have values that could not be checked against the primary paper; using them raises a `UserWarning`.
* `reduced_alphabets.json` - iFeature 5-group alphabet, cp13/cp14/cp19/cp20 cluster profiles and the six extra iFeature CTD partitions.
* `residue_tables.json` - sense-codon counts, free amino-acid formulas, pKa sets (protpy, EMBOSS, Bjellqvist, IPC_protein), the N-end rule half-life table and (v1.6.0) the Guruprasad DIWV instability table, copied from Biopython.
* `structural_params.json` - AGGRESCAN a3v, order/disorder-promoting residues and the PLAAC model parameters.

References
----------
\[1\]: Grantham R. Amino acid difference formula to help explain protein evolution. Science. 1974 Sep 6;185(4154):862-4. doi: 10.1126/science.185.4154.862. PMID: 4843792. <br>
\[2\]: Schneider G, Wrede P. The rational design of amino acid sequences by artificial neural networks and simulated molecular evolution: de novo design of an idealized leader peptidase cleavage site. Biophys J. 1994 Feb;66(2 Pt 1):335-44. doi: 10.1016/s0006-3495(94)80782-9. PMID: 8161687; PMCID: PMC1275700.