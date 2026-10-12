# v1.5.0 blocked items

Items from IMPLEMENTATION_PLAN.md that could not be implemented without inventing data.
Each entry names the table that is missing and where it would have to come from.

| # | Item | Missing data | Notes |
|---|---|---|---|
| 6 | `psekraac` | The 16 PseKRAAC reduced-alphabet type tables (types 1–16, every cluster size) | The Zuo et al. 2017 supplement (`bioinformatics_33_1_122_s1.doc`, fetched 2026-10-12) only lists the 16 types and their source papers; it does not contain the cluster tables. They would have to be transcribed from the 16 source papers (Liang 2015, Ogul 2007, Kosiol 2004, Zuo 2010, Melo 2006, Rakshit 2008, Susko 2007, Liu 2002, Li 2007, Peterson 2009, Cannata 2002, Li 2003) or obtained from the PseKRAAC server authors. The iFeature `PseKRAAC/` copies have no licence. |
