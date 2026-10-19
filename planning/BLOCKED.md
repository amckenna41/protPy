# v1.5.0 blocked items

Items from IMPLEMENTATION_PLAN.md that could not be implemented without inventing data.
Each entry names the table that is missing and where it would have to come from.

| # | Item | Missing data | Notes |
|---|---|---|---|
| 6 | `psekraac` | The 16 PseKRAAC reduced-alphabet type tables (types 1–16, every cluster size) | The Zuo et al. 2017 supplement (`bioinformatics_33_1_122_s1.doc`, fetched 2026-10-12) only lists the 16 types and their source papers; it does not contain the cluster tables. They would have to be transcribed from the 16 source papers (Liang 2015, Ogul 2007, Kosiol 2004, Zuo 2010, Melo 2006, Rakshit 2008, Susko 2007, Liu 2002, Li 2007, Peterson 2009, Cannata 2002, Li 2003) or obtained from the PseKRAAC server authors. The iFeature `PseKRAAC/` copies have no licence. |
| 9 | `z_scales` (Hellberg z1–z3) | Hellberg et al. 1987 Table (J Med Chem, 10.1021/jm00390a003) | Paywalled; no permissively licensed reproduction found. van Westen 2013's "Z-scales (3)" is a subset of Sandberg's 5, not Hellberg's, so it can't be substituted. `z_scales_5` is shipped. |
| 12 | `kidera_factors` | Kidera et al. 1985 10-factor table (J Protein Chem, 10.1007/BF01025492) | Paywalled; the only copy found was peptides.py (GPL, cross-check only). |
| 20 | `cruciani_properties` | Cruciani et al. 2004 3-property table (J Chemometrics, 10.1002/cem.856) | Paywalled; the only copy found was peptides.py (GPL, cross-check only). |
| 26 | `omega` | The exact Ω definition: binary alphabet encoding, blob size(s) and δ_max construction | Martin et al. 2016 only says Ω is "analogous to κ" and puts the calculation in its Supporting Information. The SI isn't open access (Europe PMC has no supplementary bundle; ACS and PMC downloads are behind browser challenges), so the definition and a reproducible reference value couldn't be checked. |
| 31 | `membrane_position` | The Eisenberg 1984 hydrophobic-moment-plot boundary lines (globular / surface / transmembrane) | Eisenberg, Schwarz et al. 1984 J Mol Biol (10.1016/0022-2836(84)90309-7) is paywalled (403). |
| 35 | `wilkinson_harrison_solubility` | λ1, λ2 (and the sign of λ2) and CV′ of the revised model | Davis et al. 1999 (Wiley) and Wilkinson & Harrison 1991 are paywalled. The Harrison lab calculator (biotech.ou.edu) fails TLS and the Internet Archive was offline. |
| Phase 4 scales | `hessa`, `moon_fleming` and `aasi` single-value scales (used by `hydrophobicity_profile`, E6 `gravy(scale=)` and E7 `hydrophobic_moment(scale=)`) | Hessa 2005 ΔG_app (Nature, paywalled); Moon & Fleming 2011 (values are in SI Table S1, not available); Juretić 2009 AASI (paywalled). Available scales: kyte_doolittle, hopp_woods, eisenberg, wimley_white, engelman, domcut, plus any AAindex1 code. |
| 47 (partial) | `signal_peptide_features`: `SP_vonHeijne_max` and `SP_small_m3_m1` | The von Heijne 1986 eukaryotic cleavage-site weight matrix | The NAR 1986 paper (PMC311474) is only a scanned PDF behind a browser challenge. `SP_small_m3_m1` depends on the matrix-selected site. Shipped: `SP_nregion_charge`, `SP_hregion_maxKD`. |
| 48 (partial) | `functional_site_counts`: monopartite/bipartite NLS and NES patterns | Verified regex consensus patterns | Kosugi 2009 (PMC2695404) defines NLS classes as activity profiles, not regexes; la Cour 2004 is behind a Cloudflare check. Shipped: N-glycosylation, ER retention and PTS1 as the PROSITE patterns PS00001, PS00014 and PS00342 (which are broader than the plan's narrower variants), plus tryptic sites. |
| 49 | `chou_fasman` | The full Chou–Fasman 1978 nucleation/extension/overlap rules and turn threshold | Chou & Fasman 1974/1978 are paywalled. Secondary sources are incomplete and disagree on the turn cutoff (7.5e-3 vs 7.5e-5). The Pα/Pβ/Pt and bend-frequency tables are in AAindex (CHOP780201–203, 212–215), but the rules can't be verified. |
| 50 | `gor` | GOR I 3×17×20 directional information tables | Garnier, Osguthorpe & Robson 1978 J Mol Biol is paywalled. |
| 51 | `coiled_coil` | Lupas 1991 MTK profile matrix and Gaussian parameters | Science 1991 is paywalled; the only machine-readable copy (ncoils) is GPL-2.0+. |
| 52 | `fold_amyloid` | Garbuzynskiy 2010 expected-contacts scale and threshold | Bioinformatics 2010 is behind a Cloudflare check and the FoldAmyloid server didn't respond. Campen 2008 lists a "FoldUnfold" contact scale, but it can't be confirmed to be the FoldAmyloid scale. |
| 54 (partial) | `antigenicity_surface`: Emini, Karplus–Schulz, Kolaskar–Tongaonkar | Emini surface-probability values/formula; Karplus–Schulz rigid-neighbour definition; Kolaskar–Tongaonkar scale | Emini 1985 (PMC255070) is a scanned PDF behind a browser challenge; Karplus & Schulz 1985 and Kolaskar & Tongaonkar 1990 are paywalled, and none of the three is on ExPASy ProtScale. Shipped: Parker (window 7) and Janin profiles. |
| 46 (partial) | `transmembrane_features`: `TM_min_dGapp` | Hessa 2005 position-independent ΔG_app scale | Nature 2005 is paywalled. The other four columns are shipped. |

## Shipped with caveats

- **Factor scales (Phase 2).** `z_scales_5`, `vhse`, `t_scales`, `st_scales`, `fasgai`, `blosum_indices` and `ms_whim` were transcribed from van Westen et al. 2013 Additional file 2 (CC-BY), because the primary papers are paywalled. Seven cells differ from peptides.py and are flagged `UNRESOLVED` in `protpy/data/aa_scales.json`:
  - FASGAI: V F1, P F6
  - ST: Q ST1, E ST2
  - VHSE: L VHSE3
  - z5: C Z3, L Z4

  Check these against the primary tables before release.
- **`distance_pair`.** The plan specified G + G²·(rng+1) outputs over distances 0..rng. Liu 2014 (Eqs. 5–12) defines d = 1..rng, giving G + G²·rng (602 for cp14, d = 3, as the paper states), and that is what is implemented. The paper's printed cp(14) lists E twice; the stray E in the first cluster is dropped, the only 14-cluster reading.
- **`kappa`.** δ_max is found by enumerating the segregated arrangements (positive block, negative block, neutral residues split over the three gaps), as the plan specifies, because the Das & Pappu SI is behind a browser challenge. It reproduces the published κ of the patterned (Glu-Lys)25 variants sv1, sv2, sv7, sv8, sv19 and sv30 to within 0.002. The irregular variants could not be read reliably from the Fig. 1 image and are not tested.
- **`isoelectric_point(pka_scale="emboss")`.** Uses EMBOSS 6.6.0 `Epk.dat`, whose N-terminal pKa differs from the 8.6 that Kozlowski 2016 Table 1 attributes to EMBOSS.
- **`flexibility`.** Uses Vihinen's average parameter set (AAindex VINM940101, the same table Biopython uses) with the window-9, 0.25-edge weighting. The neighbour-dependent sets (VINM940102–104) need Vihinen's definition of a rigid neighbour, which is paywalled, so they are not used.
- **`hydrophobic_moment(scale=None)`.** Keeps protpy's existing table, so the default output is unchanged. That table is the Eisenberg, Schwarz, Komarony & Wall 1984 J Mol Biol normalised consensus scale (A = 0.62; confirmed against ExPASy ProtScale). `scale="eisenberg"` is AAindex EISD840101 (Eisenberg 1984 Annu Rev Biochem, A = 0.25), a different published consensus scale.
- **`half_life`.** The E. coli value for an N-terminal Pro is "?" in ProtParam and is returned as NaN.
- **`linguistic_complexity`.** The Troyanskaya 2002 full text is behind a Cloudflare check. The ratio-of-sums form in the plan is used: Σ distinct k-mers / Σ min(20^k, L−k+1). This is consistent with the abstract, which says the method counts "the number of subwords present" with suffix trees. Some sources use a product-of-ratios form instead.
- **`wootton_federhen_complexity`.** Implements the SEG trigger stage only (K2 ≤ trigger windows), not SEG's extension threshold or interval optimisation, as the plan specifies.
- **`aggrescan`.** The a3v values and the hot-spot rules (HST −0.02, ≥5 residues, no Pro, trapezoidal areas) come from the paper's Additional files 1–2. The rule that picks the window size from sequence length, and the charge-corrected values for the two terminal residues, aren't in the published material. So `window` is a parameter (default 5) and the terminal residues take the edge-window value. Aβ42 gives the paper's two hot spots.
- **`fold_index`.** The default window of 51 comes from the plan; Prilusky 2005 and the FoldIndex server couldn't be accessed to confirm it. The 2.785 and 1.151 constants are confirmed via Das & Pappu 2013 (PMC).
- **`plaac`.** Ported from the MIT-licensed PLAAC Java source. It reproduces PLAAC 1.1.0's MOT3 reference output exactly (LLR 40.751, COREscore 40.751, PRDlen 296).
- **`atomic_composition`.** By default it matches Pfeature (no water removed); `subtract_water=True` gives the true peptide formula. Pfeature's W and Y formulas each have one H too many, and protpy uses the correct formulas.
