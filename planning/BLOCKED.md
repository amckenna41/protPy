# v1.5.0 blocked items

Items from IMPLEMENTATION_PLAN.md that could not be implemented without inventing data.
Each entry names the table that is missing and where it would have to come from.

| # | Item | Missing data | Notes |
|---|---|---|---|
| 6 | `psekraac` | The 16 PseKRAAC reduced-alphabet type tables (types 1–16, every cluster size) | The Zuo et al. 2017 supplement (`bioinformatics_33_1_122_s1.doc`, fetched 2026-10-12) only lists the 16 types and their source papers; it does not contain the cluster tables. They would have to be transcribed from the 16 source papers (Liang 2015, Ogul 2007, Kosiol 2004, Zuo 2010, Melo 2006, Rakshit 2008, Susko 2007, Liu 2002, Li 2007, Peterson 2009, Cannata 2002, Li 2003) or obtained from the PseKRAAC server authors. The iFeature `PseKRAAC/` copies have no licence. |
| 9 | `z_scales` (Hellberg z1–z3) | Hellberg et al. 1987 Table (J Med Chem, 10.1021/jm00390a003) | Paywalled; no permissively licensed reproduction found. van Westen 2013's "Z-scales (3)" is a subset of Sandberg's 5, not Hellberg's, so it can't be substituted. `z_scales_5` is shipped. |
| 12 | `kidera_factors` | Kidera et al. 1985 10-factor table (J Protein Chem, 10.1007/BF01025492) | Paywalled; the only copy found was peptides.py (GPL, cross-check only). |
| 20 | `cruciani_properties` | Cruciani et al. 2004 3-property table (J Chemometrics, 10.1002/cem.856) | Paywalled; the only copy found was peptides.py (GPL, cross-check only). |

## Shipped with caveats

- **Factor scales (Phase 2).** `z_scales_5`, `vhse`, `t_scales`, `st_scales`, `fasgai`, `blosum_indices` and `ms_whim` were transcribed from van Westen et al. 2013 Additional file 2 (CC-BY), because the primary papers are paywalled. Seven cells differ from peptides.py and are flagged `UNRESOLVED` in `protpy/data/aa_scales.json`:
  - FASGAI: V F1, P F6
  - ST: Q ST1, E ST2
  - VHSE: L VHSE3
  - z5: C Z3, L Z4

  Check these against the primary tables before release.
- **`distance_pair`.** The plan specified G + G²·(rng+1) outputs over distances 0..rng. Liu 2014 (Eqs. 5–12) defines d = 1..rng, giving G + G²·rng (602 for cp14, d = 3, as the paper states), and that is what is implemented. The paper's printed cp(14) lists E twice; the stray E in the first cluster is dropped, the only 14-cluster reading.
- **`atomic_composition`.** By default it matches Pfeature (no water removed); `subtract_water=True` gives the true peptide formula. Pfeature's W and Y formulas each have one H too many, and protpy uses the correct formulas.
