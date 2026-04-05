# protpy Descriptor Reference

A complete reference for all protein descriptors available in the **protpy** package, grouped by category.


## Composition Descriptors

### Amino Acid Composition (AAComp)

Fraction of each of the 20 standard amino acids in the sequence.

```python
result = protpy.amino_acid_composition(protein_seq)
# Shape: 1 x 20
# A      C      D      E      F  ...
# 6.693  3.108  5.817  3.347  6.614 ...
```

---

### Dipeptide Composition (DPComp)

Frequency of all 400 possible dipeptide (two-residue) combinations.

```python
result = protpy.dipeptide_composition(protein_seq)
# Shape: 1 x 400
# AA    AC    AD    AE    AF  ...
# 0.72  0.16  0.48  0.4   0.24 ...
```

---

### Tripeptide Composition (TPComp)

Frequency of all 8000 possible tripeptide (three-residue) combinations.

```python
result = protpy.tripeptide_composition(protein_seq)
# Shape: 1 x 8000
# AAA  AAC  AAD  AAE  AAF ...
# 1    0    0    2    0 ...
```

---

### Grand Average of Hydropathy (GRAVY)

Mean of the Kyte-Doolittle hydropathy values across all residues. A positive value indicates overall hydrophobicity; negative indicates overall hydrophilicity.

```python
result = protpy.gravy(protein_seq)
# Shape: 1 x 1
# GRAVY
# -0.045
```

---

### Aromaticity

Fraction of aromatic residues (F, W, Y, H) in the sequence.

```python
result = protpy.aromaticity(protein_seq)
# Shape: 1 x 1
# Aromaticity
# 0.118
```

---

### Instability Index

Stability classifier based on dipeptide instability weight values (DIWV). Values below 40 indicate a stable protein; 40 or above indicates instability.

```python
result = protpy.instability_index(protein_seq)
# Shape: 1 x 1
# InstabilityIndex
# 31.836
```

---

### Isoelectric Point

Estimated pH at which the protein carries no net charge, calculated iteratively using standard pKa values for ionisable residues.

```python
result = protpy.isoelectric_point(protein_seq)
# Shape: 1 x 1
# IsoelectricPoint
# 5.412
```

---

### Molecular Weight

Average molecular weight of the protein calculated from residue masses, corrected for water lost at each peptide bond.

```python
result = protpy.molecular_weight(protein_seq)
# Shape: 1 x 1
# MolecularWeight (Da)
# 139122.355
```

---

### Charge Distribution

Positive, negative, and net charge contributions of ionisable residues at a given pH using the Henderson-Hasselbalch equation.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `ph` | float | `7.4` | pH at which to calculate charge |

```python
# Default pH 7.4
result = protpy.charge_distribution(protein_seq)

# Custom pH
result = protpy.charge_distribution(protein_seq, ph=6.0)

# Shape: 1 x 3
# PositiveCharge  NegativeCharge  NetCharge
# 99.526          114.956         -15.43
```

---

### Hydrophobic/Polar/Charged Composition (HPC)

Percentage of residues belonging to each of three physicochemical groups: hydrophobic (A, C, F, I, L, M, V, W, Y), polar (G, N, Q, S, T), and charged (D, E, H, K, R).

```python
result = protpy.hydrophobic_polar_charged_composition(protein_seq)
# Shape: 1 x 3
# Hydrophobic  Polar   Charged
# 44.542       32.669  18.247
```

---

### Secondary Structure Propensity (SSP)

Average Chou-Fasman propensity values for alpha-helix, beta-sheet, and random coil conformations across all residues.

```python
result = protpy.secondary_structure_propensity(protein_seq)
# Shape: 1 x 3
# Helix  Sheet  Coil
# 0.983  1.05   1.043
```

---

### k-mer Composition

Frequency of all possible k-length residue subsequences, expressed as a percentage of total k-mers.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `k` | int | `2` | Length of each k-mer |

```python
# Default k=2 (dipeptide frequencies)
result = protpy.kmer_composition(protein_seq)

# Custom k
result = protpy.kmer_composition(protein_seq, k=3)

# Shape: 1 x 20^k  (e.g. 1 x 400 for k=2)
# AA     AC     AD  ...
# 0.797  0.159  ... ...
```

---

### Reduced Alphabet Composition

Amino acid composition after mapping residues to a reduced alphabet of physicochemical groups. Supported alphabet sizes: `2`, `3`, `4`, `6`.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `alphabet_size` | int | `6` | Number of reduced groups |

```python
# Default alphabet_size=6
result = protpy.reduced_alphabet_composition(protein_seq)

# Custom size
result = protpy.reduced_alphabet_composition(protein_seq, alphabet_size=4)

# Shape: 1 x alphabet_size
# Group_1  Group_2  Group_3  Group_4  Group_5  Group_6
# 25.339   34.741   9.163    9.084    10.837   10.837
```

---

### Motif Composition

Count of occurrences (including overlapping) of biological sequence motifs matched via regular expressions. Eight built-in motifs are used by default; a custom list can be supplied.

**Default motifs:**
| Column | Pattern | Biological meaning |
|--------|---------|-------------------|
| `NxST_glycosylation` | `N[^P][ST]` | N-linked glycosylation site |
| `RGD_integrin` | `RGD` | Integrin-binding RGD motif |
| `KDEL_retention` | `KDEL` | ER retention signal |
| `CxxC_zinc_finger` | `C..C` | Zinc-finger CxxC motif |
| `CAAX_prenylation` | `C[A-Z]{2}[CSIM]$` | CAAX prenylation box |
| `cAMP_PKA` | `[RK]{2}.[ST]` | cAMP/PKA phosphorylation site |
| `dileucine_sorting` | `[DE]xxxL[LI]` | Dileucine lysosomal sorting signal |
| `PEST_region` | `P.{1,10}[ED]` | PEST degradation signal |

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `motifs` | list or None | `None` | Custom regex patterns; uses built-in 8 if `None` |

```python
# Default built-in motifs
result = protpy.motif_composition(protein_seq)

# Custom motif list
result = protpy.motif_composition(protein_seq, motifs=[r'RGD', r'N[^P][ST]'])

# Shape: 1 x len(motifs)
# NxST_glycosylation  RGD_integrin  KDEL_retention  CxxC_zinc_finger  ...
# 23                  0             0               2                 ...
```

---

### Amino Acid Pair Composition

Frequency of all 400 residue-pair combinations with column names annotated by the physicochemical class of each residue (Hydrophobic, Polar, Charged, or Other).

```python
result = protpy.amino_acid_pair_composition(protein_seq)
# Shape: 1 x 400
# AA_Hydrophobic-Hydrophobic  AA_Hydrophobic-Polar  AA_Hydrophobic-Charged  ...
# 0.797                       0.159                 ...                     ...
```

---

### Aliphatic Index

A measure of the relative volume occupied by aliphatic side chains (Ala, Val, Ile, Leu). Higher values indicate greater thermostability. Formula: AI = Ala% + 2.9×Val% + 3.9×(Ile%+Leu%).

```python
result = protpy.aliphatic_index(protein_seq)
# Shape: 1 x 1
# AliphaticIndex
# 82.725
```

---

### Extinction Coefficient

The molar extinction coefficient at 280 nm, calculated from the number of Trp (W), Tyr (Y), and Cys (C) residues. Reported for both reduced (no disulfide bonds) and oxidized (all Cys paired) states.

```python
result = protpy.extinction_coefficient(protein_seq)
# Shape: 1 x 2
# ExtCoeff_Reduced  ExtCoeff_Oxidized
# 140960            143335
```

---

### Boman Index

Sum of solubility values for amino acids divided by sequence length, predicting potential for protein–protein interactions. Positive values suggest membrane-binding or interaction potential.

```python
result = protpy.boman_index(protein_seq)
# Shape: 1 x 1
# BomanIndex
# 0.119
```

---

### Aggregation Propensity

Estimates aggregation-prone regions via a sliding-window approach combining Kyte–Doolittle hydrophobicity and charge neutrality. Returns the count of qualifying windows and the fraction of the sequence covered.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `window` | int | `5` | Sliding window size |
| `hydrophobicity_threshold` | float | `2.0` | Minimum mean hydrophobicity |
| `charge_threshold` | int | `1` | Maximum charged residues per window |

```python
result = protpy.aggregation_propensity(protein_seq)
# Shape: 1 x 2
# AggregProneRegions  AggregProneFraction
# 58                  11.793
```

---

### Hydrophobic Moment

The mean and maximum hydrophobic moment across sliding windows, using the Eisenberg hydrophobicity scale and a helical-wheel projection. Captures amphipathicity of putative helix segments.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `window` | int | `11` | Sliding window size |
| `angle` | int | `100` | Residue rotation angle in degrees (100° for α-helix) |

```python
result = protpy.hydrophobic_moment(protein_seq)
# Shape: 1 x 2
# HydrophobicMoment_Mean  HydrophobicMoment_Max
# 0.272                   0.813
```

---

### Shannon Entropy

Information-theoretic measure of amino acid diversity in a sequence. Computed as:

$$H = -\sum_i p_i \log_2 p_i$$

where $p_i$ is the fractional frequency of each amino acid type present. A value of 0 indicates a completely repetitive (single amino acid) sequence; the theoretical maximum of $\log_2(20) \approx 4.322$ bits corresponds to a perfectly uniform distribution across all 20 canonical amino acids. Widely used as a low-complexity filter and diversity measure in ML feature pipelines.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `sequence` | str | — | Protein sequence |

```python
result = protpy.shannon_entropy(protein_seq)
# Shape: 1 x 1
# ShannonEntropy
# 4.163
```

---

### Pseudo Amino Acid Composition (PAAComp)

Augmented amino acid composition that incorporates sequence-order effects via correlation factors derived from physicochemical properties. Reduces the dimensionality problem of pure sequence-order information while retaining more sequence information than simple AAComp.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `lamda` | int | `30` | Number of sequence-order correlation factors to include |
| `weight` | float | `0.05` | Weighting factor for correlation layers |
| `properties` | list | `[]` | AAIndex accession numbers to use (uses built-in set if empty) |

```python
# Default parameters
result = protpy.pseudo_amino_acid_composition(protein_seq)

# Custom parameters
result = protpy.pseudo_amino_acid_composition(protein_seq, lamda=10, weight=0.1)

# Shape: 1 x (20 + lamda)  →  1 x 50 with defaults
# PAAC_1  PAAC_2  PAAC_3  ...
# 0.127   0.059   0.111   ...
```

---

### Amphiphilic Pseudo Amino Acid Composition (APAAComp)

Extension of PAAComp that uses both hydrophobicity and hydrophilicity properties to capture amphiphilic patterns (dual hydrophobic/hydrophilic character) along the sequence.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `lamda` | int | `30` | Number of sequence-order correlation factors |
| `weight` | float | `0.5` | Weighting factor for correlation layers |
| `properties` | list | `[]` | AAIndex accession numbers (defaults to hydrophobicity + hydrophilicity) |

```python
# Default parameters
result = protpy.amphiphilic_pseudo_amino_acid_composition(protein_seq)

# Custom parameters
result = protpy.amphiphilic_pseudo_amino_acid_composition(protein_seq, lamda=15, weight=0.3)

# Shape: 1 x (20 + 2*lamda)  →  1 x 80 with defaults
# APAAC_1  APAAC_2  APAAC_3  ...
# 6.624    3.076    5.757    ...
```

---

## Autocorrelation Descriptors

Autocorrelation descriptors measure the correlation between physicochemical property values of residues separated by a lag distance along the sequence. By default, 8 AAIndex properties are used, generating `lag × 8 = 240` features.

**Default properties:**
| AAIndex ID | Property |
|------------|----------|
| `CIDH920105` | Normalised average hydrophobicity |
| `BHAR880101` | Average flexibility indices |
| `CHAM820101` | Polarizability parameter |
| `CHAM820102` | Free energy of solution in water |
| `CHOC760101` | Residue accessible surface area in tripeptide |
| `BIGC670101` | Residue volume |
| `CHAM810101` | Steric parameter |
| `DAYM780201` | Relative mutability |

---

### Moreaubroto Autocorrelation (MBAuto)

Uses raw property values as the basis for correlation measurement.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `lag` | int | `30` | Maximum lag distance |
| `properties` | list | (8 defaults above) | AAIndex accession numbers |
| `normalize` | bool | `True` | Normalise output values |

```python
# Default parameters
result = protpy.moreaubroto_autocorrelation(protein_seq)

# Custom parameters
result = protpy.moreaubroto_autocorrelation(protein_seq, lag=15, properties=["CIDH920105"])

# Shape: 1 x (lag × len(properties))  →  1 x 240 with defaults
# MBAuto_CIDH920105_1  MBAuto_CIDH920105_2  ...
# -0.052               -0.104               ...
```

---

### Moran Autocorrelation (MAuto)

Uses the deviation from the mean property value, making it mean-centred and thereby less sensitive to the absolute property scale.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `lag` | int | `30` | Maximum lag distance |
| `properties` | list | (8 defaults above) | AAIndex accession numbers |
| `normalize` | bool | `True` | Normalise output values |

```python
# Default parameters
result = protpy.moran_autocorrelation(protein_seq)

# Custom parameters
result = protpy.moran_autocorrelation(protein_seq, lag=15)

# Shape: 1 x (lag × len(properties))  →  1 x 240 with defaults
# MAuto_CIDH920105_1  MAuto_CIDH920105_2  ...
# -0.07786            -0.07879            ...
```

---

### Geary Autocorrelation (GAuto)

Uses squared differences between property values at each lag, making it sensitive to local dissimilarities rather than global correlation.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `lag` | int | `30` | Maximum lag distance |
| `properties` | list | (8 defaults above) | AAIndex accession numbers |
| `normalize` | bool | `True` | Normalise output values |

```python
# Default parameters
result = protpy.geary_autocorrelation(protein_seq)

# Custom parameters
result = protpy.geary_autocorrelation(protein_seq, lag=10, normalize=False)

# Shape: 1 x (lag × len(properties))  →  1 x 240 with defaults
# GAuto_CIDH920105_1  GAuto_CIDH920105_2  ...
# 1.057               1.077               ...
```

---

## Conjoint Triad Descriptor

### Conjoint Triad (CTriad)

Encodes the sequence using a 7-class reduced amino acid alphabet and computes the frequency of all consecutive three-residue combinations (triads). The 7 classes are: (1) AGV, (2) ILFP, (3) YMTS, (4) HNQW, (5) RK, (6) DE, (7) C.

```python
result = protpy.conjoint_triad(protein_seq)
# Shape: 1 x 343  (7 × 7 × 7 class combinations)
# 111  112  113  114  ...
# 7    17   11   3    ...
```

---

## CTD Descriptors

CTD (Composition, Transition, Distribution) descriptors characterise the distribution of residues belonging to three physicochemical classes along the sequence. Seven physicochemical properties are supported.

**Supported properties:**

| Property key | Description | Classes |
|---|---|---|
| `hydrophobicity` | Hydrophobicity | Polar / Neutral / Hydrophobic |
| `normalized_vdwv` | Normalised van der Waals volume | 0–2.78 / 2.95–4.0 / 4.03–8.08 |
| `polarity` | Polarity | 4.9–6.2 / 8.0–9.2 / 10.4–13.0 |
| `charge` | Charge | Positive / Neutral / Negative |
| `secondary_struct` | Secondary structure | Helix / Strand / Coil |
| `solvent_accessibility` | Solvent accessibility | Buried / Exposed / Intermediate |
| `polarizability` | Polarizability | 0–0.108 / 0.128–0.186 / 0.219–0.409 |

---

### CTD Composition

Fraction of residues in each of the three physicochemical classes.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `property` | str | `"hydrophobicity"` | Physicochemical property to use |

```python
result = protpy.ctd_composition(protein_seq)
result = protpy.ctd_composition(protein_seq, property="charge")
# hydrophobicity_CTD_C_01  hydrophobicity_CTD_C_02  hydrophobicity_CTD_C_03
# 0.279                    0.386                    0.335
```

---

### CTD Transition

Fraction of transitions between each pair of the three physicochemical classes.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `property` | str | `"hydrophobicity"` | Physicochemical property to use |

```python
result = protpy.ctd_transition(protein_seq)
result = protpy.ctd_transition(protein_seq, property="polarity")
# hydrophobicity_CTD_T_12  hydrophobicity_CTD_T_13  hydrophobicity_CTD_T_23
# 0.181                    0.161                    0.179
```

---

### CTD Distribution

Position of the first, 25%, 50%, 75%, and last residue of each class within the sequence (as a percentage of total length).

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `property` | str | `"hydrophobicity"` | Physicochemical property to use |

```python
result = protpy.ctd_distribution(protein_seq)
result = protpy.ctd_distribution(protein_seq, property="secondary_struct")
# hydrophobicity_CTD_D_01_01  hydrophobicity_CTD_D_02_01  ...
# 0.0796                      0.557                       ...
```

---

### CTD Combined (`ctd_`)

Calculate Composition, Transition **and** Distribution for one or all supported properties in a single call.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `property` | str | `"hydrophobicity"` | Property to use when `all_ctd=False` |
| `all_ctd` | bool | `True` | If `True`, compute CTD for all 7 supported properties |

```python
# All 7 properties (default)
result = protpy.ctd_(protein_seq)

# Single property
result = protpy.ctd_(protein_seq, property="charge", all_ctd=False)

# Shape: 1 x (3 + 3 + 15) per property  →  1 x 147 for all 7 properties
# hydrophobicity_CTD_C_01  hydrophobicity_CTD_C_02  ...
# 0.279                    0.386                    ...
```

---

## Sequence Order Descriptors

Sequence order descriptors capture the effect of residue spacing along the sequence using physicochemical distance matrices. Two distance matrices are supported:

| Matrix | Description |
|--------|-------------|
| `schneider-wrede` | Physicochemical distance based on Schneider-Wrede scale (default) |
| `grantham` | Physicochemical distance based on Grantham's amino acid difference formula |

---

### Sequence Order Coupling Number — single (`sequence_order_coupling_number_`)

Computes the sum of squared physicochemical distances between all residue pairs separated by a gap of `d`. Returns a single float.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `d` | int | `1` | Gap between residue pairs |
| `distance_matrix` | str | `"schneider-wrede"` | Distance matrix to use |

```python
result = protpy.sequence_order_coupling_number_(protein_seq)
result = protpy.sequence_order_coupling_number_(protein_seq, d=5, distance_matrix="grantham")
# Returns: 401.387  (float)
```

---

### Sequence Order Coupling Number — series (`sequence_order_coupling_number`)

Calculates SOCN values across all gaps from 1 to `lag`.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `lag` | int | `30` | Maximum gap value |
| `distance_matrix` | str | `"schneider-wrede"` | Distance matrix to use |

```python
# Default parameters
result = protpy.sequence_order_coupling_number(protein_seq)

# Custom lag and matrix
result = protpy.sequence_order_coupling_number(protein_seq, lag=10, distance_matrix="grantham")

# Shape: 1 x lag  →  1 x 30 with defaults
# SOCN_SW1   SOCN_SW2   SOCN_SW3  ...
# 401.387    409.243    376.946   ...
```

---

### Sequence Order Coupling Number — all matrices (`sequence_order_coupling_number_all`)

Calculates SOCN across all lags using **both** the Schneider-Wrede and Grantham matrices and concatenates the results.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `lag` | int | `30` | Maximum gap value |

```python
result = protpy.sequence_order_coupling_number_all(protein_seq)
result = protpy.sequence_order_coupling_number_all(protein_seq, lag=15)
# Shape: 1 x (2 × lag)  →  1 x 60 with defaults
# SOCN_SW1  ...  SOCN_Grant1  ...
```

---

### Quasi Sequence Order (`quasi_sequence_order`)

Extends SOCN by combining standard amino acid composition with sequence-order coupling numbers, weighted by a factor `weight`. Captures both residue type and spatial distribution information.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `lag` | int | `30` | Maximum lag value |
| `weight` | float | `0.1` | Weighting factor for coupling terms |
| `distance_matrix` | str | `"schneider-wrede"` | Distance matrix to use |

```python
# Default parameters
result = protpy.quasi_sequence_order(protein_seq)

# Custom parameters
result = protpy.quasi_sequence_order(protein_seq, lag=10, weight=0.2, distance_matrix="grantham")

# Shape: 1 x (20 + lag)  →  1 x 50 with defaults
# QSO_SW1    QSO_SW2    QSO_SW3   ...
# 0.005692   0.002643   0.004947  ...
```

---

### Quasi Sequence Order — all matrices (`quasi_sequence_order_all`)

Calculates Quasi Sequence Order using **both** distance matrices and concatenates the results.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `lag` | int | `30` | Maximum lag value |
| `weight` | float | `0.1` | Weighting factor for coupling terms |

```python
result = protpy.quasi_sequence_order_all(protein_seq)
result = protpy.quasi_sequence_order_all(protein_seq, lag=15, weight=0.05)
# Shape: 1 x (2 × (20 + lag))  →  1 x 100 with defaults
# QSO_SW1  ...  QSO_Grant1  ...
```

---

## Descriptor Summary

Speed ratings reflect typical computation time for a single average-length protein (~500 residues) on a standard CPU:

| Rating | Meaning |
|--------|---------|
| ✅ Fast | < 1 ms — simple residue counting or scalar formula |
| 🟡 Moderate | 1–50 ms — sliding window or O(n²) pass |
| 🔴 Slow | > 50 ms — large feature space, iterative convergence, or many property lookups |

> Autocorrelation, PseAAC, and APAAComp scale with both sequence length **and** `lag` — reduce `lag` or the number of properties to speed them up. Tripeptide composition always produces 8000 columns regardless of sequence length.

| Descriptor | Function | Output shape | Category | Speed | Complexity |
|---|---|---|---|---|---|
| Amino Acid Composition | `amino_acid_composition(seq)` | 1 × 20 | Composition | ✅ Fast | O(n) |
| Dipeptide Composition | `dipeptide_composition(seq)` | 1 × 400 | Composition | ✅ Fast | O(n) |
| Tripeptide Composition | `tripeptide_composition(seq)` | 1 × 8000 | Composition | 🟡 Moderate | O(n) — large output |
| GRAVY | `gravy(seq)` | 1 × 1 | Composition | ✅ Fast | O(n) |
| Aromaticity | `aromaticity(seq)` | 1 × 1 | Composition | ✅ Fast | O(n) |
| Instability Index | `instability_index(seq)` | 1 × 1 | Composition | ✅ Fast | O(n) |
| Isoelectric Point | `isoelectric_point(seq)` | 1 × 1 | Composition | 🟡 Moderate | O(n · iterations) |
| Molecular Weight | `molecular_weight(seq)` | 1 × 1 | Composition | ✅ Fast | O(n) |
| Charge Distribution | `charge_distribution(seq, ph=7.4)` | 1 × 3 | Composition | ✅ Fast | O(n) |
| Hydrophobic/Polar/Charged | `hydrophobic_polar_charged_composition(seq)` | 1 × 3 | Composition | ✅ Fast | O(n) |
| Secondary Structure Propensity | `secondary_structure_propensity(seq)` | 1 × 3 | Composition | ✅ Fast | O(n) |
| k-mer Composition | `kmer_composition(seq, k=2)` | 1 × 20^k | Composition | 🟡 Moderate | O(n · 20^k) |
| Reduced Alphabet Composition | `reduced_alphabet_composition(seq, alphabet_size=6)` | 1 × alphabet_size | Composition | ✅ Fast | O(n) |
| Motif Composition | `motif_composition(seq, motifs=None)` | 1 × len(motifs) | Composition | 🟡 Moderate | O(n · m) per motif |
| Amino Acid Pair Composition | `amino_acid_pair_composition(seq)` | 1 × 400 | Composition | ✅ Fast | O(n) |
| Aliphatic Index | `aliphatic_index(seq)` | 1 × 1 | Composition | ✅ Fast | O(n) |
| Extinction Coefficient | `extinction_coefficient(seq)` | 1 × 2 | Composition | ✅ Fast | O(n) |
| Boman Index | `boman_index(seq)` | 1 × 1 | Composition | ✅ Fast | O(n) |
| Aggregation Propensity | `aggregation_propensity(seq, window=5)` | 1 × 2 | Composition | 🟡 Moderate | O(n · window) |
| Hydrophobic Moment | `hydrophobic_moment(seq, window=11, angle=100)` | 1 × 2 | Composition | 🟡 Moderate | O(n · window) |
| Shannon Entropy | `shannon_entropy(seq)` | 1 × 1 | Composition | ✅ Fast | O(n) |
| Pseudo AAComp | `pseudo_amino_acid_composition(seq, lamda=30, weight=0.05)` | 1 × (20 + lamda) | Composition | 🔴 Slow | O(n · lamda · props) |
| Amphiphilic Pseudo AAComp | `amphiphilic_pseudo_amino_acid_composition(seq, lamda=30, weight=0.5)` | 1 × (20 + 2×lamda) | Composition | 🔴 Slow | O(n · lamda · props) |
| Moreaubroto Autocorrelation | `moreaubroto_autocorrelation(seq, lag=30)` | 1 × (lag × props) | Autocorrelation | 🔴 Slow | O(n · lag · props) |
| Moran Autocorrelation | `moran_autocorrelation(seq, lag=30)` | 1 × (lag × props) | Autocorrelation | 🔴 Slow | O(n · lag · props) |
| Geary Autocorrelation | `geary_autocorrelation(seq, lag=30)` | 1 × (lag × props) | Autocorrelation | 🔴 Slow | O(n · lag · props) |
| Conjoint Triad | `conjoint_triad(seq)` | 1 × 343 | Conjoint Triad | ✅ Fast | O(n) |
| CTD Composition | `ctd_composition(seq, property="hydrophobicity")` | 1 × 3 | CTD | ✅ Fast | O(n) |
| CTD Transition | `ctd_transition(seq, property="hydrophobicity")` | 1 × 3 | CTD | ✅ Fast | O(n) |
| CTD Distribution | `ctd_distribution(seq, property="hydrophobicity")` | 1 × 15 | CTD | ✅ Fast | O(n) |
| CTD Combined | `ctd_(seq, property="hydrophobicity", all_ctd=True)` | 1 × 147 | CTD | 🟡 Moderate | O(n · props) |
| SOCN (single) | `sequence_order_coupling_number_(seq, d=1)` | float | Sequence Order | ✅ Fast | O(n) |
| SOCN (series) | `sequence_order_coupling_number(seq, lag=30)` | 1 × lag | Sequence Order | 🟡 Moderate | O(n · lag) |
| SOCN (all matrices) | `sequence_order_coupling_number_all(seq, lag=30)` | 1 × (2 × lag) | Sequence Order | 🟡 Moderate | O(n · lag) |
| Quasi Sequence Order | `quasi_sequence_order(seq, lag=30, weight=0.1)` | 1 × (20 + lag) | Sequence Order | 🟡 Moderate | O(n · lag) |
| Quasi Sequence Order (all) | `quasi_sequence_order_all(seq, lag=30, weight=0.1)` | 1 × (2 × (20 + lag)) | Sequence Order | 🟡 Moderate | O(n · lag) |

---

## References

The descriptors implemented in protpy are based on the following published methods:

**Composition**
- Amino acid, dipeptide, and tripeptide composition: Nakashima, H., Nishikawa, K., & Ooi, T. (1986). The folding type of a protein is relevant to the amino acid composition. *Journal of Biochemistry*, 99(1), 153–162.
- GRAVY: Kyte, J., & Doolittle, R. F. (1982). A simple method for displaying the hydropathic character of a protein. *Journal of Molecular Biology*, 157(1), 105–132.
- Aromaticity: Lobry, J. R., & Gautier, C. (1994). Hydrophobicity, expressivity and aromaticity are the major trends of amino-acid usage in 999 *Escherichia coli* chromosome-encoded genes. *Nucleic Acids Research*, 22(15), 3174–3180.
- Instability index: Guruprasad, K., Reddy, B. V. B., & Pandit, M. W. (1990). Correlation between stability of a protein and its dipeptide composition: a novel approach for predicting in vivo stability of a protein from its primary sequence. *Protein Engineering*, 4(2), 155–161.
- Isoelectric point: Bjellqvist, B., et al. (1994). The focusing positions of polypeptides in immobilized pH gradients can be predicted from their amino acid sequences. *Electrophoresis*, 14(1), 1023–1031.
- Molecular weight: Gasteiger, E., et al. (2005). Protein identification and analysis tools on the ExPASy server. In *The Proteomics Protocols Handbook*, Humana Press, 571–607.
- Isoelectric point, molecular weight, charge: Gasteiger, E., et al. (2005). *The Proteomics Protocols Handbook*, Humana Press, 571–607.
- Secondary structure propensity: Chou, P. Y., & Fasman, G. D. (1974). Prediction of protein conformation. *Biochemistry*, 13(2), 222–245.
- Aliphatic index: Ikai, A. J. (1980). Thermostability and aliphatic index of globular proteins. *Journal of Biochemistry*, 88(6), 1895–1898.
- Extinction coefficient: Pace, C. N., et al. (1995). How to measure and predict the molar absorption coefficient of a protein. *Protein Science*, 4(11), 2411–2423.
- Boman index: Boman, H. G. (2003). Antibacterial peptides: basic facts and emerging concepts. *Journal of Internal Medicine*, 254(3), 197–215.
- Hydrophobic moment: Eisenberg, D., Weiss, R. M., & Terwilliger, T. C. (1982). The helical hydrophobic moment: a measure of the amphiphilicity of a helix. *Nature*, 299, 371–374.
- Shannon entropy: Shannon, C. E. (1948). A mathematical theory of communication. *Bell System Technical Journal*, 27(3), 379–423.
- Pseudo amino acid composition (PseAAC): Chou, K.-C. (2001). Prediction of protein cellular attributes using pseudo-amino acid composition. *Proteins: Structure, Function, and Bioinformatics*, 43(3), 246–255.
- Amphiphilic PseAAC (APseAAC): Chou, K.-C. (2005). Using amphiphilic pseudo amino acid composition to predict enzyme subfamily classes. *Bioinformatics*, 21(1), 10–19.

**Autocorrelation**
- Moreau-Broto autocorrelation: Moreau, G., & Broto, P. (1980). The autocorrelation of a topological structure: A new molecular descriptor. *Nouveau Journal de Chimie*, 4, 359–360.
- Moran autocorrelation: Moran, P. A. P. (1950). Notes on continuous stochastic phenomena. *Biometrika*, 37(1–2), 17–23.
- Geary autocorrelation: Geary, R. C. (1954). The contiguity ratio and statistical mapping. *The Incorporated Statistician*, 5(3), 115–145.
- AAIndex properties: Kawashima, S., & Kanehisa, M. (2000). AAindex: amino acid index database. *Nucleic Acids Research*, 28(1), 374.

**Conjoint Triad**
- Liu, B., et al. (2008). Prediction of protein-protein interactions based on the naive Bayes classifier with amino acid composition features. *Biochemical and Biophysical Research Communications*, 368(2), 462–468. Doi: 10.1016/j.bbrc.2008.01.082.

**CTD**
- Dubchak, I., et al. (1995). Prediction of protein folding class using global description of amino acid sequence. *PNAS*, 92(19), 8700–8704.
- Dubchak, I., et al. (1999). Recognition of a protein fold in the context of the SCOP classification. *Proteins*, 35(4), 401–407.

**Sequence Order**
- Grantham, R. (1974). Amino acid difference formula to help explain protein evolution. *Science*, 185(4154), 862–864.
- Schneider, G., & Wrede, P. (1994). The rational design of amino acid sequences by artificial neural networks and simulated molecular evolution. *Biophysical Journal*, 66(2), 335–344.
- Chou, K.-C. (2000). Prediction of protein subcellular locations by incorporating quasi-sequence-order effect. *Biochemical and Biophysical Research Communications*, 278(2), 477–483.

[Back to top](#TOP)