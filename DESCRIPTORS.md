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
| `charge_threshold` | int | `0` | Maximum charged residues per window |

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

---

## v1.5.0 Descriptors

Descriptors added in v1.5.0. Each lists its formula, output shape, parameters and primary reference.
`descriptor_set="extended"` in `all_descriptors` / `encode_dataset` adds every descriptor below marked **extended**; `"all"` adds the rest as well.
Items from the implementation plan that could not be implemented without unverifiable data are listed in `planning/BLOCKED.md`.

### Composition

#### `cksaap` (CKSAAP)

Composition of k-spaced amino acid pairs: CKSAAP(XY,g) = N_g(XY) / (L-g-1) x 100 for g = 0..gap.

```text
Calculate the Composition of k-Spaced Amino Acid Pairs (CKSAAP) of a protein
sequence. For each spacing g = 0..gap, every residue pair (s_i, s_{i+g+1})
is counted and divided by the number of such pairs:

CKSAAP(XY, g) = N_g(XY) / (L - g - 1) x 100

g = 0 is the dipeptide composition. Values are percentages; a spacing with no
pairs (L - g - 1 <= 0) gives zeros.
```

```python
result = protpy.cksaap(protein_seq)
# Shape: 1 x 400·(gap+1)  (default 1 x 2400)
```

- **Columns:** `CKSAAP_AA_gap0`, `CKSAAP_AC_gap0`, `CKSAAP_AD_gap0`, `CKSAAP_AE_gap0`, `CKSAAP_AF_gap0`, `CKSAAP_AG_gap0`, …
- **Parameters:** `gap`: int (default=5) — maximum spacing
- **Descriptor set:** all
- **Reference:** Chen, Kurgan & Ruan (2007). BMC Struct. Biol., 7, 25. https://doi.org/10.1186/1472-6807-7-25

---

#### `cksaagp` (CKSAAGP)

Composition of k-spaced amino acid group pairs over the iFeature 5-group alphabet.

```text
Calculate the Composition of k-Spaced Amino Acid Group Pairs (CKSAAGP).
As CKSAAP, but residues are first mapped to the five iFeature groups:
aliphatic (GAVLMI), aromatic (FYW), positive (KRH), negative (DE) and
uncharged (STCPNQ).

CKSAAGP(G1G2, g) = N_g(G1G2) / (L - g - 1) x 100

Values are percentages and the 25 values for each spacing sum to 100
(0 when the spacing has no pairs).
```

```python
result = protpy.cksaagp(protein_seq)
# Shape: 1 x 25·(gap+1)  (default 1 x 150)
```

- **Columns:** `CKSAAGP_aliphatic.aliphatic_gap0`, `CKSAAGP_aliphatic.aromatic_gap0`, `CKSAAGP_aliphatic.positive_gap0`, `CKSAAGP_aliphatic.negative_gap0`, `CKSAAGP_aliphatic.uncharged_gap0`, `CKSAAGP_aromatic.aliphatic_gap0`, …
- **Parameters:** `gap`: int (default=5) — maximum spacing
- **Descriptor set:** extended, all
- **Reference:** Chen et al. (2018). Bioinformatics, 34(14), 2499-2502 (iFeature). https://doi.org/10.1093/bioinformatics/bty140

---

#### `dde` (DDE)

Dipeptide deviation from expected mean: (Dc - Tm) / sqrt(Tv), with Tm from codon counts.

```text
Calculate the Dipeptide Deviation from Expected Mean (DDE). The observed
dipeptide composition is standardised against the composition expected from
the number of codons coding each residue:

Dc(XY) = N(XY) / (L - 1)
Tm(XY) = (C_X / 61) x (C_Y / 61)
Tv(XY) = Tm(XY) x (1 - Tm(XY)) / (L - 1)
DDE(XY) = (Dc(XY) - Tm(XY)) / sqrt(Tv(XY))

where C_X is the number of sense codons for residue X in the standard genetic
code (61 in total). DDE is a standardised statistic, not a percentage.
```

```python
result = protpy.dde(protein_seq)
# Shape: 1 x 400
```

- **Columns:** `DDE_AA`, `DDE_AC`, `DDE_AD`, `DDE_AE`, `DDE_AF`, `DDE_AG`, …
- **Parameters:** sequence only
- **Descriptor set:** all
- **Reference:** Saravanan & Gautham (2015). OMICS, 19(10), 648-658. https://doi.org/10.1089/omi.2015.0095

---

#### `asdc` (ASDC)

Adaptive skip dipeptide composition: composition of all ordered pairs (s_i, s_j), i < j.

```text
Calculate the Adaptive Skip Dipeptide Composition (ASDC): the composition
of every ordered residue pair (s_i, s_j) with i < j, over all gaps:

ASDC(XY) = sum_g N_g(XY) / sum_g sum_XY N_g(XY) x 100

The denominator is the total number of ordered pairs, L(L-1)/2. Values are
percentages and sum to 100 (zeros for a single residue).
```

```python
result = protpy.asdc(protein_seq)
# Shape: 1 x 400
```

- **Columns:** `ASDC_AA`, `ASDC_AC`, `ASDC_AD`, `ASDC_AE`, `ASDC_AF`, `ASDC_AG`, …
- **Parameters:** sequence only
- **Descriptor set:** all
- **Reference:** Wei et al. (2018). Bioinformatics, 34(23), 4007-4016 (ACPred-FL). https://doi.org/10.1093/bioinformatics/bty451

---

#### `distance_pair` (DP)

PseAAC of distance-pairs over a reduced cluster profile (cp13/cp14/cp19/cp20).

```text
Calculate the PseAAC of distance-pairs over a reduced amino acid alphabet
(iDNA-Prot|dis). Residues are mapped to the clusters of a cluster
profile cp(n) and the vector is

[f^0(1..G), f^1(1..G^2), ..., f^d(1..G^2)]

where f^0 is the cluster composition (count / L) and f^d, d = 1..rng, is the
composition of cluster pairs (s_i, s_{i+d}) (count / (L - d)); d = 1 is the
nearest-neighbour pair. This follows Eqs. 5-12 of, giving G + G^2·rng
values (e.g. 14 + 14x14x3 = 602 for cp14, rng=3). Values are percentages.

Alphabets: cp13, cp14, cp19, cp20 (no reduction) from Eq. 10 of. The
printed cp(14) lists E twice; the stray E in the first cluster is dropped,
the only reading with 14 clusters.
```

```python
result = protpy.distance_pair(protein_seq)
# Shape: 1 x (G + G²·rng)  (default 1 x 602)
```

- **Columns:** `DP_IMV`, `DP_L`, `DP_F`, `DP_WY`, `DP_G`, `DP_P`, …
- **Parameters:** `rng`: int (default=3) — maximum pair distance; `alphabet`: str (default='cp14') — cluster profile
- **Descriptor set:** all
- **Reference:** Liu et al. (2014). PLoS ONE, 9(9), e106691 (iDNA-Prot|dis). https://doi.org/10.1371/journal.pone.0106691

---

#### `natural_vector` (NV)

Natural vector: count, mean position and normalised second central moment of each amino acid.

```text
Calculate the Natural Vector of a protein sequence. For each amino acid k:

n_k  = number of occurrences of k
mu_k = mean (1-based) position of k
D2_k = sum_i (pos_i - mu_k)^2 / (n_k x L)

mu_k and D2_k are 0 when k does not occur. Values are not percentages.
```

```python
result = protpy.natural_vector(protein_seq)
# Shape: 1 x 60
```

- **Columns:** `NV_n_A`, `NV_n_C`, `NV_n_D`, `NV_n_E`, `NV_n_F`, `NV_n_G`, …
- **Parameters:** sequence only
- **Descriptor set:** extended, all
- **Reference:** Yu et al. (2013). J. Theor. Biol., 318, 197-204. https://doi.org/10.1016/j.jtbi.2012.11.005

---

#### `atomic_composition` (ATC)

Percentage of C, H, N, O and S atoms over the residues' free amino-acid formulas.

```text
Calculate the Atomic Composition (ATC) of a protein sequence: the
percentage of C, H, N, O and S atoms, summed over the molecular formulas of
the residues' free amino acids.

ATC(e) = N_e / (N_C + N_H + N_N + N_O + N_S) x 100

By default free amino-acid formulas are summed without removing water,
matching Pfeature. With ``subtract_water=True``, (L - 1) H2O are removed
for the peptide bonds, giving the formula of the linear peptide (e.g. "GG"
-> C4H8N2O3).
```

```python
result = protpy.atomic_composition(protein_seq)
# Shape: 1 x 5
```

- **Columns:** `ATC_C`, `ATC_H`, `ATC_N`, `ATC_O`, `ATC_S`
- **Parameters:** `subtract_water`: bool (default=False) — remove (L-1) H2O
- **Descriptor set:** extended, all
- **Reference:** Pande et al. (2023). J. Comput. Biol., 30(2), 204-222 (Pfeature). https://doi.org/10.1089/cmb.2022.0241

---

### Scales

#### `z_scales_5` (Z5)

Mean of the Sandberg 5 z-scales over the sequence (optionally per domain).

```text
Calculate the mean z_scales_5 factor-scale descriptor of a protein sequence (Sandberg 5 z-scales).
See ``aa_scale_descriptors`` for the formula.
```

```python
result = protpy.z_scales_5(protein_seq)
# Shape: 1 x 5·domains  (default 1 x 5)
```

- **Columns:** `Z5_1`, `Z5_2`, `Z5_3`, `Z5_4`, `Z5_5`
- **Parameters:** `domains`: int (default=1) — number of near-equal segments
- **Descriptor set:** extended, all
- **Reference:** Sandberg et al. (1998). J. Med. Chem., 41(14), 2481-2491. https://doi.org/10.1021/jm9700575

---

#### `vhse` (VHSE)

Mean of the VHSE (hydrophobic, steric, electronic) scales over the sequence (optionally per domain).

```text
Calculate the mean vhse factor-scale descriptor of a protein sequence (principal components score Vectors of Hydrophobic, Steric and Electronic properties).
See ``aa_scale_descriptors`` for the formula.
```

```python
result = protpy.vhse(protein_seq)
# Shape: 1 x 8·domains  (default 1 x 8)
```

- **Columns:** `VHSE_1`, `VHSE_2`, `VHSE_3`, `VHSE_4`, `VHSE_5`, `VHSE_6`, …
- **Parameters:** `domains`: int (default=1) — number of near-equal segments
- **Descriptor set:** extended, all
- **Reference:** Mei et al. (2005). Biopolymers, 80(6), 775-786. https://doi.org/10.1002/bip.20296

---

#### `atchley_factors` (AF)

Mean of the Atchley factors I-V over the sequence (optionally per domain).

```text
Calculate the mean atchley_factors factor-scale descriptor of a protein sequence (Atchley factors I-V).
See ``aa_scale_descriptors`` for the formula.
```

```python
result = protpy.atchley_factors(protein_seq)
# Shape: 1 x 5·domains  (default 1 x 5)
```

- **Columns:** `AF_1`, `AF_2`, `AF_3`, `AF_4`, `AF_5`
- **Parameters:** `domains`: int (default=1) — number of near-equal segments
- **Descriptor set:** extended, all
- **Reference:** Atchley et al. (2005). PNAS, 102(18), 6395-6400. https://doi.org/10.1073/pnas.0408677102

---

#### `t_scales` (T)

Mean of the topological T-scales over the sequence (optionally per domain).

```text
Calculate the mean t_scales factor-scale descriptor of a protein sequence (topological T-scales).
See ``aa_scale_descriptors`` for the formula.
```

```python
result = protpy.t_scales(protein_seq)
# Shape: 1 x 5·domains  (default 1 x 5)
```

- **Columns:** `T_1`, `T_2`, `T_3`, `T_4`, `T_5`
- **Parameters:** `domains`: int (default=1) — number of near-equal segments
- **Descriptor set:** extended, all
- **Reference:** Tian, Zhou & Li (2007). J. Mol. Struct., 830, 106-115. https://doi.org/10.1016/j.molstruc.2006.07.004

---

#### `st_scales` (ST)

Mean of the structural topology ST-scales over the sequence (optionally per domain).

```text
Calculate the mean st_scales factor-scale descriptor of a protein sequence (structural topology ST-scales).
See ``aa_scale_descriptors`` for the formula.
```

```python
result = protpy.st_scales(protein_seq)
# Shape: 1 x 8·domains  (default 1 x 8)
```

- **Columns:** `ST_1`, `ST_2`, `ST_3`, `ST_4`, `ST_5`, `ST_6`, …
- **Parameters:** `domains`: int (default=1) — number of near-equal segments
- **Descriptor set:** extended, all
- **Reference:** Yang et al. (2010). Amino Acids, 38(3), 805-816. https://doi.org/10.1007/s00726-009-0287-y

---

#### `protfp` (ProtFP)

Mean of the ProtFP (PCA, 8 components) over the sequence (optionally per domain).

```text
Calculate the mean protfp factor-scale descriptor of a protein sequence (ProtFP PCA, 8 components).
See ``aa_scale_descriptors`` for the formula.
```

```python
result = protpy.protfp(protein_seq)
# Shape: 1 x 8·domains  (default 1 x 8)
```

- **Columns:** `ProtFP_1`, `ProtFP_2`, `ProtFP_3`, `ProtFP_4`, `ProtFP_5`, `ProtFP_6`, …
- **Parameters:** `domains`: int (default=1) — number of near-equal segments
- **Descriptor set:** extended, all
- **Reference:** van Westen et al. (2013). J. Cheminform., 5, 41. https://doi.org/10.1186/1758-2946-5-41

---

#### `fasgai` (FASGAI)

Mean of the FASGAI factor scales over the sequence (optionally per domain).

```text
Calculate the mean fasgai factor-scale descriptor of a protein sequence (factor analysis scales of generalized amino acid information).
See ``aa_scale_descriptors`` for the formula.
```

```python
result = protpy.fasgai(protein_seq)
# Shape: 1 x 6·domains  (default 1 x 6)
```

- **Columns:** `FASGAI_1`, `FASGAI_2`, `FASGAI_3`, `FASGAI_4`, `FASGAI_5`, `FASGAI_6`
- **Parameters:** `domains`: int (default=1) — number of near-equal segments
- **Descriptor set:** extended, all
- **Reference:** Liang & Li (2007). QSAR Comb. Sci., 26(6), 754-763. https://doi.org/10.1002/qsar.200630145

---

#### `blosum_indices` (BLOSUM)

Mean of the BLOSUM62-derived indices over the sequence (optionally per domain).

```text
Calculate the mean blosum_indices factor-scale descriptor of a protein sequence (BLOSUM62-derived indices).
See ``aa_scale_descriptors`` for the formula.
```

```python
result = protpy.blosum_indices(protein_seq)
# Shape: 1 x 10·domains  (default 1 x 10)
```

- **Columns:** `BLOSUM_1`, `BLOSUM_2`, `BLOSUM_3`, `BLOSUM_4`, `BLOSUM_5`, `BLOSUM_6`, …
- **Parameters:** `domains`: int (default=1) — number of near-equal segments
- **Descriptor set:** extended, all
- **Reference:** Georgiev (2009). J. Comput. Biol., 16(5), 703-723. https://doi.org/10.1089/cmb.2008.0173

---

#### `ms_whim` (MSWHIM)

Mean of the MS-WHIM scores over the sequence (optionally per domain).

```text
Calculate the mean ms_whim factor-scale descriptor of a protein sequence (MS-WHIM scores).
See ``aa_scale_descriptors`` for the formula.
```

```python
result = protpy.ms_whim(protein_seq)
# Shape: 1 x 3·domains  (default 1 x 3)
```

- **Columns:** `MSWHIM_1`, `MSWHIM_2`, `MSWHIM_3`
- **Parameters:** `domains`: int (default=1) — number of near-equal segments
- **Descriptor set:** extended, all
- **Reference:** Zaliani & Gancia (1999). J. Chem. Inf. Comput. Sci., 39(3), 525-533. https://doi.org/10.1021/ci980211b

---

#### `auto_cross_covariance` (ACC)

Auto/cross-covariance of an amino-acid scale: sum z_j(i) z_k(i+l) / (L-l).

```text
Calculate the auto- and cross-covariance (ACC) transform of an amino-acid
scale. For lag l and scale dimensions j, k:

ACC_{j,k}(l) = sum_{i=1}^{L-l} z_{j,i} x z_{k,i+l} / (L - l)

``mode="acc"`` returns every (j, k) pair, ``"ac"`` only j = k (auto
covariance) and ``"cc"`` only j != k (cross covariance). ``scale`` may be a
protpy factor-scale family (used as published) or an AAindex1 code / list
of codes, which are standardised over the 20 amino acids.
```

```python
result = protpy.auto_cross_covariance(protein_seq)
# Shape: 1 x d²·lag  (default 1 x 125)
```

- **Columns:** `ACC_1_1_lag1`, `ACC_1_2_lag1`, `ACC_1_3_lag1`, `ACC_1_4_lag1`, `ACC_1_5_lag1`, `ACC_2_1_lag1`, …
- **Parameters:** `scale`: str | list[str] (default='z_scales_5') — factor scale or AAindex code(s); `lag`: int (default=5) — maximum lag; `mode`: str (default='acc') — 'acc', 'ac' or 'cc'
- **Descriptor set:** all
- **Reference:** Wold et al. (1993). Anal. Chim. Acta, 277(2), 239-253. https://doi.org/10.1016/0003-2670(93)80437-P

---

#### `physicochemical_distance_transformation` (PDT)

Physicochemical distance transformation: mean squared difference of standardised AAindex properties at lag l.

```text
Calculate the Physicochemical Distance Transformation (PDT). For each
AAindex property P_j (standardised over the 20 amino acids) and lag l:

PDT_j(l) = sum_{i=1}^{L-l} (P_j(s_i) - P_j(s_{i+l}))^2 / (L - l)
```

```python
result = protpy.physicochemical_distance_transformation(protein_seq)
# Shape: 1 x len(properties)·lag  (default 1 x 24)
```

- **Columns:** `PDT_CIDH920105_1`, `PDT_CIDH920105_2`, `PDT_CIDH920105_3`, `PDT_BHAR880101_1`, `PDT_BHAR880101_2`, `PDT_BHAR880101_3`, …
- **Parameters:** `lag`: int (default=3) — maximum lag; `properties`: list[str] | None — AAindex codes (default: 8 autocorrelation properties)
- **Descriptor set:** extended, all
- **Reference:** Liu et al. (2012). PLoS ONE, 7(9), e46633. https://doi.org/10.1371/journal.pone.0046633

---

### Physicochemical

#### `charge_fractions` (FCR)

Fraction of charged residues (FCR), positive/negative fractions and charge asymmetry sigma = NCPR²/FCR.

```text
Calculate the charge fractions of a protein sequence:

FCR    = (n+ + n-) / L          fraction of charged residues
FPlus  = n+ / L,  FMinus = n- / L
Sigma  = NCPR^2 / FCR           overall charge asymmetry (0 if FCR = 0)

where n+ counts K/R and n- counts D/E (H is neutral). Values are fractions.
```

```python
result = protpy.charge_fractions(protein_seq)
# Shape: 1 x 4
```

- **Columns:** `Charge_FCR`, `Charge_FPlus`, `Charge_FMinus`, `Charge_Sigma`
- **Parameters:** sequence only
- **Descriptor set:** extended, all
- **Reference:** Das & Pappu (2013). PNAS, 110(33), 13392-13397. https://doi.org/10.1073/pnas.1304749110

---

#### `ncpr` (NCPR)

Net charge per residue: (n+ - n-) / L with K/R positive and D/E negative.

```text
Calculate the Net Charge Per Residue (NCPR) of a protein sequence:

NCPR = (n+ - n-) / L

with n+ counting K/R and n- counting D/E.
```

```python
result = protpy.ncpr(protein_seq)
# Shape: 1 x 1
```

- **Columns:** `NCPR`
- **Parameters:** sequence only
- **Descriptor set:** extended, all
- **Reference:** Mao et al. (2010). PNAS, 107(18), 8183-8188. https://doi.org/10.1073/pnas.0911107107

---

#### `kappa` (Kappa)

Das-Pappu kappa: normalised blob charge-asymmetry deviation measuring segregation of opposite charges (0 mixed, 1 segregated).

```text
Calculate the kappa charge-patterning parameter of Das & Pappu. For blob
size g, the sequence is split into the L - g + 1 overlapping blobs and

sigma_b = (f+ - f-)^2 / (f+ + f-)        per blob (0 if uncharged)
delta_g = mean over blobs of (sigma_b - sigma)^2
kappa_g = delta_g / delta_max,g

where sigma = NCPR^2 / FCR of the whole sequence and delta_max,g is the
largest delta over the maximally segregated arrangements of the same
composition (a positive block and a negative block with the neutral residues
distributed over the three surrounding gaps; every distinct arrangement is
evaluated). kappa is the mean of kappa_g
over ``blob_sizes`` (5 and 6 in); 0 means well mixed and 1 fully
segregated. Charges: K, R = +1; D, E = -1; H = 0.

kappa is reported as 0 (not NaN) when the sequence has no charges or only one
charge type, where oppositely charged residues cannot be segregated. A blob
size longer than the sequence is reduced to the sequence length.
```

```python
result = protpy.kappa(protein_seq)
# Shape: 1 x (1 + len(blob_sizes))  (default 1 x 3)
```

- **Columns:** `Kappa`, `Kappa_5`, `Kappa_6`
- **Parameters:** `blob_sizes`: tuple[int] (default=(5, 6)) — blob sizes averaged
- **Descriptor set:** extended, all
- **Reference:** Das & Pappu (2013). PNAS, 110(33), 13392-13397. https://doi.org/10.1073/pnas.1304749110

---

#### `sequence_charge_decoration` (SCD)

Sequence charge decoration: (1/L) sum_{m>n} q_m q_n sqrt(m-n).

```text
Calculate the Sequence Charge Decoration (SCD) of Sawle & Ghosh:

SCD = (1/L) sum_{m=2}^{L} sum_{n=1}^{m-1} q_m q_n (m - n)^{1/2}

with q = +1 for K/R, -1 for D/E and 0 otherwise. More negative values
indicate more segregated opposite charges; uncharged sequences give 0.
```

```python
result = protpy.sequence_charge_decoration(protein_seq)
# Shape: 1 x 1
```

- **Columns:** `Charge_SCD`
- **Parameters:** sequence only
- **Descriptor set:** extended, all
- **Reference:** Sawle & Ghosh (2015). J. Chem. Phys., 143(8), 085101. https://doi.org/10.1063/1.4929391

---

#### `charge_density` (ChargeDensity)

Net charge at the given pH divided by molecular weight (charges per Da).

```text
Calculate the charge density of a protein sequence: the net charge at
``ph`` (from ``charge_distribution``) divided by the molecular weight (from
``molecular_weight``), in charges per Da.

ChargeDensity = NetCharge(pH) / MW
```

```python
result = protpy.charge_density(protein_seq)
# Shape: 1 x 1
```

- **Columns:** `ChargeDensity`
- **Parameters:** `ph`: float (default=7.4) — pH
- **Descriptor set:** extended, all
- **Reference:** Müller et al. (2017). Bioinformatics, 33(17), 2753-2755 (modlAMP). https://doi.org/10.1093/bioinformatics/btx285

---

#### `hydrophobicity_profile` (HP)

Summary statistics of a sliding-window hydrophobicity profile (ProtScale-style edge weighting).

```text
Summarise a sliding-window hydrophobicity profile. Each residue is
replaced by its scale value and a window of length ``window`` is slid along
the sequence; ``edge < 1`` weights the window ends linearly down to ``edge``
(ExPASy ProtScale semantics). The profile is reduced to its mean, max, min
and standard deviation, plus the fraction of windows above ``threshold``
when one is given (e.g. 1.6 for Kyte-Doolittle transmembrane segments).
Windows longer than the sequence use the whole sequence.
```

```python
result = protpy.hydrophobicity_profile(protein_seq)
# Shape: 1 x 4 (1 x 5 with threshold)
```

- **Columns:** `HP_kyte_doolittle_w9_mean`, `HP_kyte_doolittle_w9_max`, `HP_kyte_doolittle_w9_min`, `HP_kyte_doolittle_w9_std`
- **Parameters:** `scale`: str (default='kyte_doolittle') — single-value scale or AAindex code; `window`: int (default=9) — window length; `edge`: float (default=1.0) — relative end weight; `threshold`: float | None (default=None) — adds fraction of windows above threshold
- **Descriptor set:** extended, all
- **Reference:** Kyte & Doolittle (1982). J. Mol. Biol., 157(1), 105-132. https://doi.org/10.1016/0022-2836(82)90515-0

---

#### `hydrophobic_ratio` (HydrophobicRatio)

APD hydrophobic ratio: percentage of L, I, V, F, M, A, C, W residues.

```text
Calculate the hydrophobic ratio of the Antimicrobial Peptide Database (APD)
: the percentage of residues that are L, I, V, F, M, A, C or W.

HydrophobicRatio = N(LIVFMACW) / L x 100
```

```python
result = protpy.hydrophobic_ratio(protein_seq)
# Shape: 1 x 1
```

- **Columns:** `HydrophobicRatio`
- **Parameters:** sequence only
- **Descriptor set:** extended, all
- **Reference:** Wang (2020). Protein Sci., 29(1), 8-18 (APD); Wang, Li & Wang (2016). Nucleic Acids Res., 44, D1087 (APD3). https://doi.org/10.1002/pro.3702

---

#### `half_life` (HalfLife)

N-end rule half-life (hours) in mammalian, yeast and E. coli systems from the N-terminal residue.

```text
Estimate the in vivo half-life of a protein from its N-terminal residue by
the N-end rule, as tabulated by ExPASy ProtParam: mammalian
reticulocytes (in vitro), yeast and E. coli (in vivo).

Values are in hours. Lower bounds in the source (">20 h" in yeast, ">10 h"
in E. coli, ">20 h" for mammalian Pro) are reported as the bound. The E. coli
value for an N-terminal Pro is unknown in the source and returned as NaN.
Not applicable to N-terminally modified proteins.
```

```python
result = protpy.half_life(protein_seq)
# Shape: 1 x 3
```

- **Columns:** `HalfLife_Mammal_h`, `HalfLife_Yeast_h`, `HalfLife_Ecoli_h`
- **Parameters:** sequence only
- **Descriptor set:** extended, all
- **Reference:** Bachmair, Finley & Varshavsky (1986). Science, 234(4773), 179-186. https://doi.org/10.1126/science.3018930

---

#### `ivywrel_index` (IVYWREL)

Fraction of I, V, Y, W, R, E, L residues (thermophilic adaptation index).

```text
Calculate the IVYWREL index of Zeldovich et al.: the fraction of residues
that are I, V, Y, W, R, E or L, which correlates with optimal growth
temperature.

IVYWREL = N(IVYWREL) / L
```

```python
result = protpy.ivywrel_index(protein_seq)
# Shape: 1 x 1
```

- **Columns:** `IVYWREL`
- **Parameters:** sequence only
- **Descriptor set:** extended, all
- **Reference:** Zeldovich, Berezovsky & Shakhnovich (2007). PLoS Comput. Biol., 3(1), e5. https://doi.org/10.1371/journal.pcbi.0030005

---

#### `cvp_bias` (CvP)

Charged versus polar bias: %(D+E+K+R) - %(N+Q+S+T).

```text
Calculate the charged-versus-polar (CvP) bias of Suhre & Claverie, an
indicator of hyperthermostability:

CvP = %(D + E + K + R) - %(N + Q + S + T)
```

```python
result = protpy.cvp_bias(protein_seq)
# Shape: 1 x 1
```

- **Columns:** `CvP_Bias`
- **Parameters:** sequence only
- **Descriptor set:** extended, all
- **Reference:** Suhre & Claverie (2003). J. Biol. Chem., 278(19), 17198-17202. https://doi.org/10.1074/jbc.M301327200

---

#### `flexibility` (Flex)

Summary statistics of the Vihinen average-B-value flexibility profile (window 9, end weight 0.25).

```text
Summarise the flexibility profile of Vihinen et al.. Average normalised
B-values (Vihinen's "average" parameter set, AAindex VINM940101) are averaged
over a window of 9 whose end weights fall linearly to 0.25 (weights 0.25,
0.4375, 0.625, 0.8125, 1, ..., 0.25; sum 5.25), and the profile is reduced
to its mean, max, min and standard deviation.

The neighbour-dependent parameter sets of (residues with 0, 1 or 2 rigid
neighbours) are not used. Compared with Biopython's
``ProteinAnalysis.flexibility``, which uses the same values and weights,
protpy centres each window on its middle residue and includes the last
window; Biopython takes the residue after the centre and omits the last window.
```

```python
result = protpy.flexibility(protein_seq)
# Shape: 1 x 4
```

- **Columns:** `Flex_mean`, `Flex_max`, `Flex_min`, `Flex_std`
- **Parameters:** `window`: int (default=9) — window length
- **Descriptor set:** extended, all
- **Reference:** Vihinen, Torkkila & Riikonen (1994). Proteins, 19(2), 141-149. https://doi.org/10.1002/prot.340190207

---

### Complexity

#### `wootton_federhen_complexity` (WF)

Wootton-Federhen local complexity (K1, K2) over sliding windows and the SEG-trigger low-complexity fraction.

```text
Calculate Wootton-Federhen local compositional complexity. For every
window of length W with residue counts n_i:

K1 = (1/W) log_20( W! / prod_i n_i! )             (complexity state, 0-1)
K2 = -sum_i (n_i/W) log2(n_i/W)                   (Shannon entropy, bits)

The low-complexity fraction is the fraction of residues covered by at least
one window with K2 <= ``k2_trigger``. This is the SEG trigger stage only:
SEG's extension to the K2 extension threshold and its interval merging and
optimisation are not performed. A window longer than the sequence uses the
whole sequence.
```

```python
result = protpy.wootton_federhen_complexity(protein_seq)
# Shape: 1 x 4
```

- **Columns:** `WF_K1_mean`, `WF_K2_mean`, `WF_K2_min`, `WF_LowComplexity_frac`
- **Parameters:** `window`: int (default=12) — window length; `k2_trigger`: float (default=2.2) — K2 trigger (bits)
- **Descriptor set:** extended, all
- **Reference:** Wootton & Federhen (1993). Comput. Chem., 17(2), 149-163. https://doi.org/10.1016/0097-8485(93)85006-X

---

#### `linguistic_complexity` (LC)

Linguistic complexity: distinct subwords observed / maximum possible (ratio of sums over word lengths).

```text
Calculate the linguistic complexity (LC) of a protein sequence: the
number of distinct subwords (k-mers) observed, over the maximum number
possible for a sequence of that length over the 20-letter alphabet:

LC = sum_{k=1}^{K} U_k / sum_{k=1}^{K} min(20^k, L - k + 1)

where U_k is the number of distinct k-mers and K = ``max_k`` (default L, i.e.
all subwords). This is the ratio-of-sums form, counting all subwords as in
; repetitive sequences have low LC.
```

```python
result = protpy.linguistic_complexity(protein_seq)
# Shape: 1 x 1
```

- **Columns:** `LinguisticComplexity`
- **Parameters:** `max_k`: int | None (default=None) — longest word length (default L)
- **Descriptor set:** extended, all
- **Reference:** Troyanskaya et al. (2002). Bioinformatics, 18(5), 679-688. https://doi.org/10.1093/bioinformatics/18.5.679

---

#### `lempel_ziv_complexity` (LZ)

Lempel-Ziv (LZ76) complexity c(n) and its normalised value c(n)·log20(n)/n.

```text
Calculate the Lempel-Ziv (LZ76) complexity of a protein sequence: the
number c(n) of new components found when the sequence is parsed from left to
right, each component being the shortest substring not producible by
copying from the preceding history. The normalised value for a 20-letter
alphabet is

LZ_Normalized = c(n) x log_20(n) / n

(0 for a single residue).
```

```python
result = protpy.lempel_ziv_complexity(protein_seq)
# Shape: 1 x 2
```

- **Columns:** `LZ_Complexity`, `LZ_Normalized`
- **Parameters:** sequence only
- **Descriptor set:** extended, all
- **Reference:** Lempel & Ziv (1976). IEEE Trans. Inf. Theory, 22(1), 75-81. https://doi.org/10.1109/TIT.1976.1055501

---

### Signal

#### `fft_descriptors` (FFT)

Mean power of the FFT spectrum of a mean-centred property signal (default EIIP) in equal-width frequency bins.

```text
Calculate a length-invariant binned power spectrum of a property signal, as
in the resonant recognition model of Cosic. Each residue is replaced by
its scale value (default: electron-ion interaction potential, EIIP), the
signal is mean-centred, and its power spectrum is

P(f) = |FFT(x)(f)|^2 / L,   f = 0 .. 0.5 cycles/residue

The normalised frequency range 0-0.5 is divided into ``n_bins`` equal-width
bins (the last bin includes 0.5) and the mean power in each bin is returned;
bins containing no frequency are 0.
```

```python
result = protpy.fft_descriptors(protein_seq)
# Shape: 1 x n_bins  (default 1 x 10)
```

- **Columns:** `FFT_eiip_bin1`, `FFT_eiip_bin2`, `FFT_eiip_bin3`, `FFT_eiip_bin4`, `FFT_eiip_bin5`, `FFT_eiip_bin6`, …
- **Parameters:** `scale`: str | dict (default='eiip') — single-value scale, AAindex code or dict; `n_bins`: int (default=10) — number of frequency bins
- **Descriptor set:** extended, all
- **Reference:** Cosic (1994). IEEE Trans. Biomed. Eng., 41(12), 1101-1114. https://doi.org/10.1109/10.335859

---

#### `dct_descriptors` (DCT)

First n orthonormal DCT-II coefficients of a mean-centred property signal (default Kyte-Doolittle).

```text
Calculate the first ``n_coeffs`` coefficients of the orthonormal type-II
discrete cosine transform (DCT-II) of a mean-centred property signal:

X_k = s_k sum_{n=0}^{L-1} x_n cos(pi k (2n + 1) / (2L)),
s_0 = sqrt(1/L), s_k = sqrt(2/L) for k > 0

X_0 is 0 because the signal is mean-centred.
```

```python
result = protpy.dct_descriptors(protein_seq)
# Shape: 1 x n_coeffs  (default 1 x 10)
```

- **Columns:** `DCT_kyte_doolittle_0`, `DCT_kyte_doolittle_1`, `DCT_kyte_doolittle_2`, `DCT_kyte_doolittle_3`, `DCT_kyte_doolittle_4`, `DCT_kyte_doolittle_5`, …
- **Parameters:** `scale`: str | dict (default='kyte_doolittle') — single-value scale, AAindex code or dict; `n_coeffs`: int (default=10) — number of coefficients (<= L)
- **Descriptor set:** extended, all
- **Reference:** Wang et al. (2017). PLoS ONE, 12(9), e0185587. https://doi.org/10.1371/journal.pone.0185587

---

### Structural

#### `aggrescan` (AGG)

AGGRESCAN aggregation hot-spot descriptors from the a3v propensity profile (hot spots: >=5 residues above -0.02, no Pro).

```text
Calculate AGGRESCAN aggregation "hot spot" descriptors. Each residue is
given its experimental aggregation propensity a3v; the a4v profile (AP) is the
a3v mean over a sliding window, assigned to the window's central residue.
Off-centre N- and C-terminal residues take the first and last window values.
A hot spot (HS) is a run of >= 5 residues with a4v above the hot-spot
threshold HST = -0.02 and containing no proline. Areas use trapezoidal
integration of the AP with HST as the zero axis.

Outputs:
    a3vSA   mean a3v over the sequence
    nHS     number of hot spots;     NnHS = nHS / L x 100
    AAT     area of the AP above HST over the whole sequence;  AATr = AAT / L
    THSA    total hot-spot area (sum of HSA);                   THSAr = THSA / L
    TA      total area of the AP (HST as zero axis, signed)
    a4vSS   sum of a4v;              Na4vSS = a4vSS / L x 100
    a4vAHS  mean a4v inside hot spots (mean over hot spots; 0 if none)
    HSA_max largest single hot-spot area (0 if none)

AGGRESCAN chooses the window (5, 7, 9 or 11) from the sequence length and
gives the two extreme residues specific charge-corrected values; neither
rule is given in the published supplementary material, so ``window`` is a
parameter (default 5) and the extreme residues are treated like the other
off-centre residues.
```

```python
result = protpy.aggrescan(protein_seq)
# Shape: 1 x 12
```

- **Columns:** `AGG_a3vSA`, `AGG_nHS`, `AGG_NnHS`, `AGG_AAT`, `AGG_AATr`, `AGG_THSA`, …
- **Parameters:** `window`: int (default=5) — odd window length (AGGRESCAN uses 5, 7, 9, 11)
- **Descriptor set:** extended, all
- **Reference:** Conchillo-Solé et al. (2007). BMC Bioinformatics, 8, 65. https://doi.org/10.1186/1471-2105-8-65

---

#### `uversky_charge_hydropathy` (UV)

Uversky mean hydropathy <H>, mean net charge <R> and boundary distance HQ = 2.785<H> - <R> - 1.151.

```text
Calculate the charge-hydropathy descriptors of Uversky et al.:

<H> = mean of (KD + 4.5) / 9            (Kyte-Doolittle rescaled to 0-1)
<R> = |n+ - n-| / L                     (K/R positive, D/E negative)
HQ  = 2.785 <H> - <R> - 1.151

The boundary <R> = 2.785 <H> - 1.151 separates natively unfolded from folded
proteins (as restated in); ``UV_Unfolded`` is 1 when HQ < 0.
```

```python
result = protpy.uversky_charge_hydropathy(protein_seq)
# Shape: 1 x 4
```

- **Columns:** `UV_MeanH`, `UV_MeanR`, `UV_HQ`, `UV_Unfolded`
- **Parameters:** sequence only
- **Descriptor set:** extended, all
- **Reference:** Uversky, Gillespie & Fink (2000). Proteins, 41(3), 415-427. https://doi.org/10.1002/1097-0134(20001115)41:3<415::AID-PROT130>3.0.CO;2-7

---

#### `fold_index` (FI)

FoldIndex profile summary: windowed 2.785<H> - |<R>| - 1.151 and disordered-residue statistics.

```text
Calculate FoldIndex profile descriptors. For each window of length
``window``:

FI = 2.785 <H> - |<R>| - 1.151

with <H> and <R> as in ``uversky_charge_hydropathy`` but computed over the
window. Residues covered by at least one window with FI < 0 are counted as
disordered. A window longer than the sequence uses the whole sequence, so
FI then equals ``UV_HQ``.
```

```python
result = protpy.fold_index(protein_seq)
# Shape: 1 x 4
```

- **Columns:** `FI_mean`, `FI_min`, `FI_frac_disordered`, `FI_longest_disordered`
- **Parameters:** `window`: int (default=51) — window length
- **Descriptor set:** extended, all
- **Reference:** Prilusky et al. (2005). Bioinformatics, 21(16), 3435-3438. https://doi.org/10.1093/bioinformatics/bti537

---

#### `disorder_propensity` (DP)

Fractions of disorder- and order-promoting residues and mean TOP-IDP value.

```text
Calculate disorder-propensity descriptors: the fractions of
disorder-promoting (A, R, G, Q, S, P, E, K) and order-promoting
(W, C, F, I, Y, V, L, N) residues, and the mean TOP-IDP value
(higher = more disorder-prone).
```

```python
result = protpy.disorder_propensity(protein_seq)
# Shape: 1 x 3
```

- **Columns:** `DP_frac_disorder_promoting`, `DP_frac_order_promoting`, `DP_TOP_IDP_mean`
- **Parameters:** sequence only
- **Descriptor set:** extended, all
- **Reference:** Campen et al. (2008). Protein Pept. Lett., 15(9), 956-963. https://doi.org/10.2174/092986608785849164

---

#### `transmembrane_features` (TM)

Kyte-Doolittle transmembrane-segment heuristics (segment count, max hydropathy, TM fraction, positive-inside bias).

```text
Calculate heuristic transmembrane-segment features from the Kyte-Doolittle
hydropathy profile. Windows of length ``window`` with mean hydropathy
above ``kd_threshold`` (1.6 for 19-residue windows in) mark putative
transmembrane segments; overlapping windows are merged into segments.

TM_n_segments            number of merged segments
TM_max_KD                maximum window hydropathy
TM_frac_in_TM            fraction of residues inside a segment
TM_positive_inside_bias  (K+R in the 15 residues before the first segment)
                         - (K+R in the 15 residues after it); 0 if no segment

This is a hydropathy heuristic, not a topology predictor.
```

```python
result = protpy.transmembrane_features(protein_seq)
# Shape: 1 x 4
```

- **Columns:** `TM_n_segments`, `TM_max_KD`, `TM_frac_in_TM`, `TM_positive_inside_bias`
- **Parameters:** `window`: int (default=19) — hydropathy window; `kd_threshold`: float (default=1.6) — segment threshold
- **Descriptor set:** extended, all
- **Reference:** Kyte & Doolittle (1982). J. Mol. Biol., 157(1), 105-132. https://doi.org/10.1016/0022-2836(82)90515-0

---

#### `signal_peptide_features` (SP)

N-terminal signal-peptide heuristics: n-region net charge and maximum 8-residue hydropathy.

```text
Calculate simple signal-peptide features of the N-terminus:

SP_nregion_charge  net charge of the first 5 residues (K/R +1, D/E -1)
SP_hregion_maxKD   maximum mean Kyte-Doolittle hydropathy over 8-residue
                   windows within the first ``n_term`` residues

The von Heijne cleavage-site weight matrix is not included.
```

```python
result = protpy.signal_peptide_features(protein_seq)
# Shape: 1 x 2
```

- **Columns:** `SP_nregion_charge`, `SP_hregion_maxKD`
- **Parameters:** `n_term`: int (default=30) — N-terminal region length
- **Descriptor set:** extended, all
- **Reference:** von Heijne (1983). Eur. J. Biochem., 133(1), 17-21. https://doi.org/10.1111/j.1432-1033.1983.tb07424.x

---

#### `functional_site_counts` (FS)

Counts and per-100-residue densities of PROSITE N-glycosylation, ER-retention and PTS1 motifs and tryptic sites.

```text
Count functional sequence motifs and their density per 100 residues:
N-glycosylation sites (PROSITE PS00001), the C-terminal ER-retention signal
(PS00014) and peroxisomal targeting signal 1 (PS00342), and tryptic
cleavage sites (after K/R, not before P). Overlapping matches are
counted. These complement the default motifs of ``motif_composition``.
```

```python
result = protpy.functional_site_counts(protein_seq)
# Shape: 1 x 8
```

- **Columns:** `FS_N_glycosylation`, `FS_ER_retention`, `FS_PTS1`, `FS_tryptic_site`, `FS_N_glycosylation_per100`, `FS_ER_retention_per100`, …
- **Parameters:** sequence only
- **Descriptor set:** extended, all
- **Reference:** Sigrist et al. (2013). Nucleic Acids Res., 41(D1), D344-D347 (PROSITE). https://doi.org/10.1093/nar/gks1067

---

#### `plaac` (PLAAC)

PLAAC prion-like domain scores: max window LLR, HMM core score, PrD length and presence.

```text
Calculate prion-like domain descriptors with the PLAAC model, ported from
the PLAAC source (MIT License; S. cerevisiae parameters, alpha = 1). Each
residue has the log-likelihood ratio LLR = ln(f_prion / f_background).

PLAAC_LLR        maximum LLR sum over any window of ``core_length`` residues
PLAAC_COREscore  maximum LLR sum over a ``core_length`` window lying inside a
                 prion-like run of the two-state HMM Viterbi parse
PLAAC_PrD_length length of the prion-like run containing the best core
PLAAC_has_PrD    1 if a prion-like run of >= ``core_length`` residues exists

As in PLAAC, LLR and COREscore are NaN when undefined (sequence shorter than
``core_length`` / no prion-like run of that length); PrD_length is then 0.
```

```python
result = protpy.plaac(protein_seq)
# Shape: 1 x 4
```

- **Columns:** `PLAAC_LLR`, `PLAAC_COREscore`, `PLAAC_PrD_length`, `PLAAC_has_PrD`
- **Parameters:** `core_length`: int (default=60) — minimal core length
- **Descriptor set:** extended, all
- **Reference:** Lancaster et al. (2014). Bioinformatics, 30(17), 2501-2502. https://doi.org/10.1093/bioinformatics/btu310

---

#### `antigenicity_surface` (Antigenicity)

Profile statistics of the Parker hydrophilicity and Janin surface scales over 7-residue windows.

```text
Summarise hydrophilicity/surface profiles associated with antigenic sites:
the Parker HPLC hydrophilicity scale averaged over 7-residue windows and
the Janin free energy of transfer from inside to outside of a globular
protein over the same window. Each profile is reduced to its mean, max,
min and standard deviation.
```

```python
result = protpy.antigenicity_surface(protein_seq)
# Shape: 1 x 8
```

- **Columns:** `PARKER_mean`, `PARKER_max`, `PARKER_min`, `PARKER_std`, `JANIN_mean`, `JANIN_max`, …
- **Parameters:** `window`: int (default=7) — window length
- **Descriptor set:** extended, all
- **Reference:** Parker, Guo & Hodges (1986). Biochemistry, 25(19), 5425-5432. https://doi.org/10.1021/bi00367a013

---

#### `cysteine_features` (Cys)

Cysteine count, fraction, parity, maximum disulfides, CxxC motifs and mean spacing.

```text
Calculate cysteine features:

Cys_count          number of C
Cys_frac           C / L
Cys_even           1 if the count is even
Cys_max_disulfides floor(count / 2)
Cys_CxxC           overlapping C-x-x-C motifs
Cys_mean_spacing   mean gap between consecutive cysteines (0 if fewer than 2)
```

```python
result = protpy.cysteine_features(protein_seq)
# Shape: 1 x 6
```

- **Columns:** `Cys_count`, `Cys_frac`, `Cys_even`, `Cys_max_disulfides`, `Cys_CxxC`, `Cys_mean_spacing`
- **Parameters:** sequence only
- **Descriptor set:** extended, all
- **Reference:** Fomenko & Gladyshev (2003). Biochemistry, 42(38), 11214-11225. https://doi.org/10.1021/bi034459s

---

### Extensions to existing descriptors (v1.5.0)

Defaults are unchanged; the new behaviour is opt-in.

| Function | New parameter(s) | Effect | Reference |
|---|---|---|---|
| `reduced_alphabet_composition` | `alphabet`, `k` | `alphabet="ifeature5"` with `k=1/2/3` gives GAAC (1 x 5), GDPC (1 x 25), GTPC (1 x 125); cp13/14/19/20 also accepted | Chen et al. 2018, 10.1093/bioinformatics/bty140 |
| `conjoint_triad` | `gap` | `gap>0` gives KSCTriad, triads (s_i, s_{i+g+1}, s_{i+2g+2}) for g = 0..gap, columns suffixed `_gap{g}` | Chen et al. 2018; Shen et al. 2007, 10.1073/pnas.0607879104 |
| `ctd_` | `property_set` | `"ifeature13"` adds 6 AAindex hydrophobicity partitions (13 properties, 1 x 273) | Dubchak et al. 1995, 10.1073/pnas.92.19.8700 |
| `amino_acid_composition`, `dipeptide_composition` | `region`, `n` | `"nterm"`, `"cterm"` or `"split"` (N-terminal / middle / C-terminal blocks, prefixes `NT_`/`MID_`/`CT_`) | Kumar et al. 2006, 10.1074/jbc.M511061200 |
| `isoelectric_point`, `charge_distribution` | `pka_scale` | `"emboss"`, `"bjellqvist"` (with residue-specific terminal pKa) or `"ipc_protein"` | Bjellqvist et al. 1993, 10.1002/elps.11501401163; Kozlowski 2016, 10.1186/s13062-016-0159-9 |
| `gravy` | `scale` | average any single-value scale or AAindex code (column `GRAVY_{scale}`) | per scale |
| `hydrophobic_moment` | `scale`, `angle="helix"/"sheet"` | alternative hydrophobicity scale; named angles 100° / 160° | Eisenberg et al. 1984, 10.1073/pnas.81.1.140 |
| `extinction_coefficient` | `absorbance` | adds `Abs01_Reduced`, `Abs01_Oxidized` (A280 of a 1 g/L solution, ε / MW) | Pace et al. 1995, 10.1002/pro.5560041120 |

---

## Utility Functions

### `all_descriptors`

Compute all available descriptors for a single protein sequence and return them concatenated into a single DataFrame.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `sequence` | str | — | Protein sequence |
| `params` | dict or None | `None` | Per-descriptor keyword arguments keyed by descriptor name |
| `show_progress` | bool | `False` | Display a `tqdm` progress bar naming each descriptor as it is computed |
| `output_file` | str or None | `None` | Path to save results; `None` disables file output |
| `save_each_sequence` | bool | `True` | Append this row to the output CSV immediately after computation; useful when calling in a loop |
| `descriptor_set` | str | `"classic"` | `"classic"` (v1.4 set, unchanged), `"extended"` (adds the v1.5 descriptors with ≤150 columns) or `"all"` (adds every v1.5 descriptor) |
| `include` | list[str] or None | `None` | Extra registered descriptors to compute |
| `exclude` | list[str] or None | `None` | Descriptors to skip (classic ones included, e.g. `"tripeptide_composition"`) |

```python
import protpy

seq = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDISIQEFKD"

# v1.5 descriptor sets
df = protpy.all_descriptors(seq, descriptor_set="extended")
df = protpy.all_descriptors(seq, include=["kappa", "plaac"], exclude=["tripeptide_composition"],
                            params={"kappa": {"blob_sizes": (5,)}})

# Basic usage
df = protpy.all_descriptors(seq)
# Shape: 1 x 10757 (all descriptor features concatenated)

# With progress bar
df = protpy.all_descriptors(seq, show_progress=True)
# [Composition]    Amino acid composition  10%|███     |  3/30 [00:00<00:02]
# ...

# Save to file
df = protpy.all_descriptors(seq, output_file="out.csv")

# No file output
df = protpy.all_descriptors(seq, output_file=None)
```

---

### `encode_dataset`

Compute all available descriptors for every sequence in a dataset (list or DataFrame) and return a DataFrame whose first column is a unique sequence identifier.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `sequences` | list[str] or pd.DataFrame | — | Sequences as a plain list or a two-column DataFrame from `read_sequences_from_file(..., id_col=...)` |
| `seq_col` | str | `"sequence"` | Column name containing sequences when input is a DataFrame |
| `params` | dict or None | `None` | Per-descriptor keyword arguments forwarded to each `all_descriptors` call |
| `show_progress` | bool | `False` | Show a per-descriptor `tqdm` bar for each individual sequence |
| `output_file` | str or None | `"descriptors.csv"` | CSV path; each completed row is flushed immediately so interrupted runs preserve progress |
| `save_each_sequence` | bool | `True` | Flush the CSV after each sequence; set `False` to skip all file output within this call |
| `strip_gaps` | bool | `False` | Remove gap characters (`-`) before encoding; useful for MSA-derived sequences |

```python
import protpy

# Plain list — identifiers auto-generated as seq1, seq2, …
seqs = [
    "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGD",
    "ACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWY",
]
df = protpy.encode_dataset(seqs)
# df.columns[0] == "sequence_id"
# df["sequence_id"].tolist() == ["seq1", "seq2"]
# Shape: 2 x (1 + 10757)

# Named identifiers from a multi-column file
data = protpy.read_sequences_from_file(
    "thermostability.txt", seq_col="sequence", id_col="name"
)
df = protpy.encode_dataset(data, seq_col="sequence")
# df.columns[0] == "name"
# Shape: n_sequences x (1 + 10757)

# MSA sequences with gap characters
df = protpy.encode_dataset(seqs, strip_gaps=True)
```

---

### `read_sequences_from_file`

Parse protein sequences from a FASTA, TXT, CSV, or TSV file. The format is inferred from the file extension.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `filepath` | str or os.PathLike | — | Path to the input file |
| `seq_col` | str, int, or None | `None` | Column containing sequences (name or zero-based index); ignored for plain FASTA and plain TXT |
| `id_col` | str, int, or None | `None` | Identifier column; when provided returns a DataFrame instead of a list |
| `has_header` | bool | `True` | Whether the first row of a delimited file is a header; columns must be selected by integer index when `False` |

**Supported formats:**

| Extension | Behaviour |
|-----------|-----------|
| `.fasta` / `.fa` | Standard FASTA; multi-line sequences joined; returns `list[str]` |
| `.txt` (no `seq_col`) | One raw sequence per non-blank line; returns `list[str]` |
| `.txt` (with `seq_col`) | Treated as delimited; separator auto-detected from first line |
| `.csv` | Comma-separated values |
| `.tsv` | Tab-separated values |

When `id_col` is provided the function returns a `pd.DataFrame` with only the identifier and sequence columns; all other columns are dropped and rows with blank sequences are excluded. When `id_col` is omitted the function returns a plain `list[str]`.

**Raises:** `FileNotFoundError` if the file does not exist; `ValueError` for unsupported extensions or missing columns; `TypeError` if `seq_col` or `id_col` is not a str, int, or None.

```python
import protpy

# FASTA
seqs = protpy.read_sequences_from_file("proteins.fasta")
# ['MKTAYIAKQRQISFVKSHFS...', 'ACDEFGHIKLM...', ...]

# Plain TXT (one sequence per line)
seqs = protpy.read_sequences_from_file("sequences.txt")

# CSV / TSV — by column name
seqs = protpy.read_sequences_from_file("data.csv", seq_col="sequence")

# CSV / TSV — by column index
seqs = protpy.read_sequences_from_file("data.tsv", seq_col=1)

# Multi-column file — keep id + sequence, drop everything else
df = protpy.read_sequences_from_file(
    "thermostability.txt", seq_col="sequence", id_col="name"
)
# Returns DataFrame:
#        name                           sequence
# 0  00000000  MTIKEMPQPKTFGELKNLPLLNTD...
# 1  11111111  KETSPIPQPKTFGPLGNLPLIDKD...

# Headerless file — select by index
seqs = protpy.read_sequences_from_file("seqs.csv", has_header=False, seq_col=0)
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

**v1.5.0 descriptors**
- Each v1.5.0 descriptor above lists its primary reference and DOI. BibTeX for every reference is in [`docs/references.bib`](docs/references.bib).

**Sequence Order**
- Grantham, R. (1974). Amino acid difference formula to help explain protein evolution. *Science*, 185(4154), 862–864.
- Schneider, G., & Wrede, P. (1994). The rational design of amino acid sequences by artificial neural networks and simulated molecular evolution. *Biophysical Journal*, 66(2), 335–344.
- Chou, K.-C. (2000). Prediction of protein subcellular locations by incorporating quasi-sequence-order effect. *Biochemical and Biophysical Research Communications*, 278(2), 477–483.

[Back to top](#TOP)