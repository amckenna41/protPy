################################################################################
#############                      Utilities                      ##############
################################################################################

from __future__ import annotations

"""
Utility and auxiliary functions for the protpy package. Includes descriptor
metadata lookup, shared helpers, and other convenience tools.
"""

#registry of all available descriptors with metadata
_DESCRIPTOR_REGISTRY: dict[str, dict] = {
    # ── Composition ────────────────────────────────────────────────────────────
    "amino_acid_composition": {
        "description": (
            "Fraction of each of the 20 canonical amino acids in the sequence. "
            "AA_Comp(s) = AA(t) / N(s)."
        ),
        "module": "composition",
        "category": "Composition",
        "abbreviation": "AAComp",
        "parameters": {
            "sequence": "str — protein sequence",
        },
        "output_shape": "1 x 20",
        "speed": "Fast",
        "reference": "[1] Gromiha, M. M. (2010). Protein Bioinformatics. Elsevier.",
    },
    "dipeptide_composition": {
        "description": (
            "Fraction of each of the 400 possible dipeptide types (20²) in the sequence. "
            "DPComp(s,t) = AA(s,t) / (N - 1)."
        ),
        "module": "composition",
        "category": "Composition",
        "abbreviation": "DPComp",
        "parameters": {
            "sequence": "str — protein sequence",
        },
        "output_shape": "1 x 400",
        "speed": "Fast",
        "reference": "[1] Gromiha, M. M. (2010). Protein Bioinformatics. Elsevier.",
    },
    "tripeptide_composition": {
        "description": (
            "Frequency of each of the 8000 possible tripeptide types (20³) in the sequence."
        ),
        "module": "composition",
        "category": "Composition",
        "abbreviation": "TPComp",
        "parameters": {
            "sequence": "str — protein sequence",
        },
        "output_shape": "1 x 8000",
        "speed": "Moderate",
        "reference": "[1] Gromiha, M. M. (2010). Protein Bioinformatics. Elsevier.",
    },
    "gravy": {
        "description": (
            "Grand Average of Hydropathy (GRAVY) — sum of Kyte-Doolittle hydropathy "
            "values divided by sequence length. Positive = hydrophobic, negative = hydrophilic."
        ),
        "module": "composition",
        "category": "Composition",
        "abbreviation": "GRAVY",
        "parameters": {
            "sequence": "str — protein sequence",
        },
        "output_shape": "1 x 1",
        "speed": "Fast",
        "reference": "[8] Kyte & Doolittle (1982). J. Mol. Biol., 157(1), 105-132.",
    },
    "aromaticity": {
        "description": (
            "Fraction of aromatic amino acids (Phe, Trp, Tyr) in the sequence."
        ),
        "module": "composition",
        "category": "Composition",
        "abbreviation": "Aromaticity",
        "parameters": {
            "sequence": "str — protein sequence",
        },
        "output_shape": "1 x 1",
        "speed": "Fast",
        "reference": "[9] Lobry & Gautier (1994). Nucleic Acids Research, 22(15), 3174-3180.",
    },
    "instability_index": {
        "description": (
            "Instability Index — weighted sum of dipeptide instability weight values. "
            "Values < 40 predict in vivo stability; >= 40 predicts instability."
        ),
        "module": "composition",
        "category": "Composition",
        "abbreviation": "II",
        "parameters": {
            "sequence": "str — protein sequence",
        },
        "output_shape": "1 x 1",
        "speed": "Fast",
        "reference": "[10] Guruprasad et al. (1990). Protein Engineering, 4(2), 155-161.",
    },
    "isoelectric_point": {
        "description": (
            "Theoretical isoelectric point (pI) — the pH at which the net protein "
            "charge is zero, estimated via iterative charge-balance."
        ),
        "module": "composition",
        "category": "Composition",
        "abbreviation": "pI",
        "parameters": {
            "sequence": "str — protein sequence",
        },
        "output_shape": "1 x 1",
        "speed": "Fast",
        "reference": "[11] Bjellqvist et al. (1993). Electrophoresis, 14(1), 1023-1031.",
    },
    "molecular_weight": {
        "description": (
            "Molecular Weight (MW) — sum of average residue masses minus "
            "one water molecule per peptide bond (18.015 Da)."
        ),
        "module": "composition",
        "category": "Composition",
        "abbreviation": "MW",
        "parameters": {
            "sequence": "str — protein sequence",
        },
        "output_shape": "1 x 1",
        "speed": "Fast",
        "reference": "[12] Gasteiger et al. (2005). The Proteomics Protocols Handbook.",
    },
    "charge_distribution": {
        "description": (
            "Positive, negative and net charge at a given pH, calculated using "
            "Henderson-Hasselbalch equations for K, R, H, D, E residues."
        ),
        "module": "composition",
        "category": "Composition",
        "abbreviation": "ChargeDist",
        "parameters": {
            "sequence": "str — protein sequence",
            "ph": "float (default=7.4) — pH for charge calculation",
        },
        "output_shape": "1 x 3",
        "speed": "Fast",
        "reference": "[13] Cameselle-Teijeiro (1979). Biochemical Education, 7(3), 69-70.",
    },
    "hydrophobic_polar_charged_composition": {
        "description": (
            "Percentage composition of residues grouped into three physicochemical "
            "classes: Hydrophobic (A,C,F,I,L,M,V,W,Y), Polar (G,N,Q,S,T), "
            "and Charged (D,E,R,H,K)."
        ),
        "module": "composition",
        "category": "Composition",
        "abbreviation": "HPC",
        "parameters": {
            "sequence": "str — protein sequence",
        },
        "output_shape": "1 x 3",
        "speed": "Fast",
        "reference": "[1] Gromiha, M. M. (2010). Protein Bioinformatics. Elsevier.",
    },
    "secondary_structure_propensity": {
        "description": (
            "Mean Chou-Fasman propensity values for alpha-helix, beta-sheet and "
            "coil conformations across the entire sequence."
        ),
        "module": "composition",
        "category": "Composition",
        "abbreviation": "SSP",
        "parameters": {
            "sequence": "str — protein sequence",
        },
        "output_shape": "1 x 3",
        "speed": "Fast",
        "reference": "[14] Chou & Fasman (1974). Biochemistry, 13(2), 211-222.",
    },
    "kmer_composition": {
        "description": (
            "Frequency of all 20^k possible k-length amino acid subsequences in the "
            "sequence, expressed as a percentage of total k-mers. k <= 4."
        ),
        "module": "composition",
        "category": "Composition",
        "abbreviation": "kmer",
        "parameters": {
            "sequence": "str — protein sequence",
            "k": "int (default=2) — k-mer length (1 ≤ k ≤ 4)",
        },
        "output_shape": "1 x 20^k  (e.g. 1 x 400 for k=2)",
        "speed": "Fast (k≤3) / Slow (k=4)",
        "reference": "[1] Gromiha, M. M. (2010). Protein Bioinformatics. Elsevier.",
    },
    "reduced_alphabet_composition": {
        "description": (
            "Fractional composition after mapping 20 amino acids to a reduced "
            "alphabet of 2, 3, 4 or 6 physicochemically similar groups."
        ),
        "module": "composition",
        "category": "Composition",
        "abbreviation": "ReducedAlphabet",
        "parameters": {
            "sequence": "str — protein sequence",
            "alphabet_size": "int (default=6) — number of groups (2, 3, 4 or 6)",
        },
        "output_shape": "1 x alphabet_size",
        "speed": "Fast",
        "reference": "[1] Gromiha, M. M. (2010). Protein Bioinformatics. Elsevier.",
    },
    "motif_composition": {
        "description": (
            "Count of regex-matched occurrences for each motif in the sequence. "
            "Uses a built-in set of 8 biologically relevant motifs by default "
            "(N-glycosylation, RGD, KDEL, CxxC, CAAX, PKA, dileucine, PEST)."
        ),
        "module": "composition",
        "category": "Composition",
        "abbreviation": "MotifComp",
        "parameters": {
            "sequence": "str — protein sequence",
            "motifs": "dict[str, str] | None (default=None) — {name: regex_pattern}; None uses built-in defaults",
        },
        "output_shape": "1 x M  (M = number of motifs, default M=8)",
        "speed": "Fast",
        "reference": "[15] PROSITE: https://prosite.expasy.org/",
    },
    "amino_acid_pair_composition": {
        "description": (
            "Fractional frequency of all 400 consecutive residue pairs, annotated "
            "with the physicochemical class (Hydrophobic/Polar/Charged) of each residue."
        ),
        "module": "composition",
        "category": "Composition",
        "abbreviation": "PairComp",
        "parameters": {
            "sequence": "str — protein sequence",
        },
        "output_shape": "1 x 400",
        "speed": "Moderate",
        "reference": "[1] Gromiha, M. M. (2010). Protein Bioinformatics. Elsevier.",
    },
    "aliphatic_index": {
        "description": (
            "Aliphatic Index (AI) — relative volume of aliphatic side chains "
            "(Ala, Val, Ile, Leu). Higher values indicate thermostability. "
            "AI = Xala + 2.9·Xval + 3.9·(Xile + Xleu)."
        ),
        "module": "composition",
        "category": "Composition",
        "abbreviation": "AI",
        "parameters": {
            "sequence": "str — protein sequence",
        },
        "output_shape": "1 x 1",
        "speed": "Fast",
        "reference": "[16] Ikai (1980). J. Biochemistry, 88(6), 1895-1898.",
    },
    "extinction_coefficient": {
        "description": (
            "Molar extinction coefficient at 280 nm for both reduced "
            "(no disulfide bonds) and oxidized (SS bonds) forms, using Trp, Tyr, Cys counts."
        ),
        "module": "composition",
        "category": "Composition",
        "abbreviation": "ExtCoeff",
        "parameters": {
            "sequence": "str — protein sequence",
        },
        "output_shape": "1 x 2  [ExtCoeff_Reduced, ExtCoeff_Oxidized]",
        "speed": "Fast",
        "reference": "[17] Gasteiger et al. (2005). The Proteomics Protocols Handbook.",
    },
    "boman_index": {
        "description": (
            "Boman Index (potential protein interaction index) — sum of residue "
            "solubility values divided by sequence length. Values > 2.48 suggest "
            "high protein-binding potential."
        ),
        "module": "composition",
        "category": "Composition",
        "abbreviation": "BI",
        "parameters": {
            "sequence": "str — protein sequence",
        },
        "output_shape": "1 x 1",
        "speed": "Fast",
        "reference": "[18] Boman (2003). J. Internal Medicine, 254(3), 197-215.",
    },
    "aggregation_propensity": {
        "description": (
            "Number of aggregation-prone regions (APRs) and the fraction of "
            "residues covered by any APR, identified by high local hydrophobicity "
            "and low charge within a sliding window."
        ),
        "module": "composition",
        "category": "Composition",
        "abbreviation": "AggregProp",
        "parameters": {
            "sequence": "str — protein sequence",
            "window": "int (default=5) — sliding window length",
            "hydrophobicity_threshold": "float (default=2.0) — minimum mean Kyte-Doolittle hydrophobicity",
            "charge_threshold": "int (default=1) — max charged residues allowed in an APR window",
        },
        "output_shape": "1 x 2  [AggregProneRegions, AggregProneFraction]",
        "speed": "Fast",
        "reference": "[8] Kyte & Doolittle (1982). J. Mol. Biol., 157(1), 105-132.",
    },
    "hydrophobic_moment": {
        "description": (
            "Eisenberg hydrophobic moment (μH) — measure of amphiphilicity calculated "
            "as the vector sum of residue hydrophobicities projected at a given helical "
            "angle. Averaged over all windows of the given size."
        ),
        "module": "composition",
        "category": "Composition",
        "abbreviation": "μH",
        "parameters": {
            "sequence": "str — protein sequence",
            "window": "int (default=11) — helical window size",
            "angle": "float (default=100) — residue rotation angle in degrees (100° for alpha-helix)",
        },
        "output_shape": "1 x 1",
        "speed": "Fast",
        "reference": "[19] Eisenberg et al. (1982). Nature, 299, 371-374.",
    },
    "shannon_entropy": {
        "description": (
            "Shannon entropy of the amino acid composition — measure of sequence "
            "complexity/diversity. H = -sum(p_i * log2(p_i)) over all 20 amino acids."
        ),
        "module": "composition",
        "category": "Composition",
        "abbreviation": "ShannonH",
        "parameters": {
            "sequence": "str — protein sequence",
        },
        "output_shape": "1 x 1",
        "speed": "Fast",
        "reference": "[1] Gromiha, M. M. (2010). Protein Bioinformatics. Elsevier.",
    },
    "pseudo_amino_acid_composition": {
        "description": (
            "Pseudo Amino Acid Composition (PAAComp) — augments standard AAComp with "
            "sequence-order correlation factors derived from physicochemical properties. "
            "First 20 values are weighted AAComp; remaining lamda values are correlation factors."
        ),
        "module": "composition",
        "category": "Composition",
        "abbreviation": "PAAComp",
        "parameters": {
            "sequence": "str — protein sequence",
            "lamda": "int (default=30) — rank of correlation (output = 20 + lamda features)",
            "weight": "float (default=0.05) — weighting factor for correlation terms",
            "properties": "list[str] (default=[]) — AAIndex1 accession codes; empty = use built-in hydrophobicity/hydrophilicity/mass",
        },
        "output_shape": "1 x (20 + lamda)  (default 1 x 50)",
        "speed": "Moderate",
        "reference": "[4] Chou, K-C. (2001). PROTEINS, 43, 246-255.",
    },
    "amphiphilic_pseudo_amino_acid_composition": {
        "description": (
            "Amphiphilic Pseudo Amino Acid Composition (APAAComp) — extends PAAComp "
            "with hydrophobicity and hydrophilicity distribution patterns. First 20 "
            "values are weighted AAComp; remaining 2·lamda values encode amphiphilic correlations."
        ),
        "module": "composition",
        "category": "Composition",
        "abbreviation": "APAAComp",
        "parameters": {
            "sequence": "str — protein sequence",
            "lamda": "int (default=30) — rank of correlation (output = 20 + 2·lamda features)",
            "weight": "float (default=0.5) — weighting factor",
            "properties": "list[dict] (default=[hydrophobicity_, hydrophilicity_]) — property value dicts",
        },
        "output_shape": "1 x (20 + 2·lamda)  (default 1 x 80)",
        "speed": "Moderate",
        "reference": "[5] Chou, K-C. (2005). Bioinformatics, 21(1), 10-19.",
    },
    # ── Autocorrelation ────────────────────────────────────────────────────────
    "moreaubroto_autocorrelation": {
        "description": (
            "MoreauBroto Autocorrelation (MBAuto) — topological descriptor based on "
            "the product of property values for residue pairs separated by lag positions. "
            "Generates lag × n_properties features."
        ),
        "module": "autocorrelation",
        "category": "Autocorrelation",
        "abbreviation": "MBAuto",
        "parameters": {
            "sequence": "str — protein sequence",
            "lag": "int (default=30) — max lag between residue pairs",
            "properties": "list[str] (default=8 AAIndex codes) — AAIndex1 accession numbers",
            "normalize": "bool (default=True) — normalise property values before calculation",
        },
        "output_shape": "1 x (lag × n_properties)  (default 1 x 240)",
        "speed": "Moderate",
        "reference": "[1] Hollas (2003). J. Math. Chem., 33(2), 91-101.",
    },
    "moran_autocorrelation": {
        "description": (
            "Moran Autocorrelation (MAuto) — uses property deviations from the mean "
            "rather than raw values to compute autocorrelation. "
            "Generates lag × n_properties features."
        ),
        "module": "autocorrelation",
        "category": "Autocorrelation",
        "abbreviation": "MAuto",
        "parameters": {
            "sequence": "str — protein sequence",
            "lag": "int (default=30) — max lag between residue pairs",
            "properties": "list[str] (default=8 AAIndex codes) — AAIndex1 accession numbers",
            "normalize": "bool (default=True) — normalise property values",
        },
        "output_shape": "1 x (lag × n_properties)  (default 1 x 240)",
        "speed": "Moderate",
        "reference": "[2] Ong et al. (2007). BMC Bioinformatics, 8, 300.",
    },
    "geary_autocorrelation": {
        "description": (
            "Geary Autocorrelation (GAuto) — uses the squared difference of property "
            "values for residue pairs separated by lag, instead of vector products. "
            "Generates lag × n_properties features."
        ),
        "module": "autocorrelation",
        "category": "Autocorrelation",
        "abbreviation": "GAuto",
        "parameters": {
            "sequence": "str — protein sequence",
            "lag": "int (default=30) — max lag between residue pairs",
            "properties": "list[str] (default=8 AAIndex codes) — AAIndex1 accession numbers",
            "normalize": "bool (default=True) — normalise property values",
        },
        "output_shape": "1 x (lag × n_properties)  (default 1 x 240)",
        "speed": "Moderate",
        "reference": "[2] Ong et al. (2007). BMC Bioinformatics, 8, 300.",
    },
    # ── Conjoint Triad ─────────────────────────────────────────────────────────
    "conjoint_triad": {
        "description": (
            "Conjoint Triad (CTriad) — encodes each sequence using a 7-letter reduced "
            "alphabet and computes the frequency distribution of all 7³=343 triads "
            "(consecutive amino acid triplets)."
        ),
        "module": "conjoint_triad",
        "category": "Conjoint Triad",
        "abbreviation": "CTriad",
        "parameters": {
            "sequence": "str — protein sequence",
        },
        "output_shape": "1 x 343",
        "speed": "Fast",
        "reference": "[1] Shen et al. (2007). PNAS, 104(11), 4337-4341.",
    },
    # ── CTD ────────────────────────────────────────────────────────────────────
    "ctd_composition": {
        "description": (
            "CTD Composition (CTD-C) — fraction of residues belonging to each of three "
            "property-defined classes (class 1, 2 and 3). 3 features per property."
        ),
        "module": "ctd",
        "category": "CTD",
        "abbreviation": "CTD-C",
        "parameters": {
            "sequence": "str — protein sequence",
            "property": "str (default='hydrophobicity') — one of: hydrophobicity, normalized_vdwv, polarity, charge, secondary_struct, solvent_accessibility, polarizability",
        },
        "output_shape": "1 x 3  (per property)",
        "speed": "Fast",
        "reference": "[1] Dubchak et al. (1995). PNAS, 92, 8700-8704.",
    },
    "ctd_transition": {
        "description": (
            "CTD Transition (CTD-T) — percentage frequency of transitions between "
            "different property classes for consecutive residues. 3 features per property."
        ),
        "module": "ctd",
        "category": "CTD",
        "abbreviation": "CTD-T",
        "parameters": {
            "sequence": "str — protein sequence",
            "property": "str (default='hydrophobicity') — CTD physicochemical property",
        },
        "output_shape": "1 x 3  (per property)",
        "speed": "Fast",
        "reference": "[1] Dubchak et al. (1995). PNAS, 92, 8700-8704.",
    },
    "ctd_distribution": {
        "description": (
            "CTD Distribution (CTD-D) — chain positions where the 1st, 25%, 50%, 75% "
            "and 100% of residues belonging to a property class are located. "
            "15 features per property."
        ),
        "module": "ctd",
        "category": "CTD",
        "abbreviation": "CTD-D",
        "parameters": {
            "sequence": "str — protein sequence",
            "property": "str (default='hydrophobicity') — CTD physicochemical property",
        },
        "output_shape": "1 x 15  (per property)",
        "speed": "Fast",
        "reference": "[1] Dubchak et al. (1995). PNAS, 92, 8700-8704.",
    },
    "ctd_": {
        "description": (
            "Combined CTD descriptor (CTD-C + CTD-T + CTD-D) — concatenation of "
            "Composition, Transition and Distribution for one or all 7 properties. "
            "Single property = 21 features; all properties = 147 features."
        ),
        "module": "ctd",
        "category": "CTD",
        "abbreviation": "CTD",
        "parameters": {
            "sequence": "str — protein sequence",
            "property": "str (default='hydrophobicity') — CTD property (used when all_ctd=False)",
            "all_ctd": "bool (default=True) — if True, compute over all 7 properties",
        },
        "output_shape": "1 x 21  (single property) / 1 x 147  (all 7 properties)",
        "speed": "Fast (single) / Moderate (all)",
        "reference": "[1] Dubchak et al. (1995). PNAS, 92, 8700-8704.",
    },
    # ── Sequence Order ─────────────────────────────────────────────────────────
    "sequence_order_coupling_number_": {
        "description": (
            "Single Sequence Order Coupling Number (SOCN) — dissimilarity between "
            "all amino acid pairs separated by gap d in the sequence, using a "
            "selected physicochemical distance matrix. Returns a single float."
        ),
        "module": "sequence_order",
        "category": "Sequence Order",
        "abbreviation": "SOCN_d",
        "parameters": {
            "sequence": "str — protein sequence",
            "d": "int (default=1) — gap between amino acid pairs",
            "distance_matrix": "str (default='schneider-wrede') — 'schneider-wrede' or 'grantham'",
        },
        "output_shape": "scalar float",
        "speed": "Fast",
        "reference": "[1] Chou, K-C. (2000). Biochem. Biophys. Res. Commun., 278, 477-483.",
    },
    "sequence_order_coupling_number": {
        "description": (
            "Sequence Order Coupling Number (SOCN) — SOCN values for gaps 1 to lag, "
            "using one physicochemical distance matrix. N = lag features."
        ),
        "module": "sequence_order",
        "category": "Sequence Order",
        "abbreviation": "SOCN",
        "parameters": {
            "sequence": "str — protein sequence",
            "lag": "int (default=30) — maximum gap; output features = lag",
            "distance_matrix": "str (default='schneider-wrede') — 'schneider-wrede' or 'grantham'",
        },
        "output_shape": "1 x lag  (default 1 x 30)",
        "speed": "Moderate",
        "reference": "[1] Chou, K-C. (2000). Biochem. Biophys. Res. Commun., 278, 477-483.",
    },
    "sequence_order_coupling_number_all": {
        "description": (
            "SOCN using both distance matrices (schneider-wrede + grantham) concatenated "
            "into a single DataFrame. Output features = lag × 2."
        ),
        "module": "sequence_order",
        "category": "Sequence Order",
        "abbreviation": "SOCN_all",
        "parameters": {
            "sequence": "str — protein sequence",
            "lag": "int (default=30) — maximum gap; output features = lag × 2",
        },
        "output_shape": "1 x (lag × 2)  (default 1 x 60)",
        "speed": "Moderate",
        "reference": "[1] Chou, K-C. (2000). Biochem. Biophys. Res. Commun., 278, 477-483.",
    },
    "quasi_sequence_order": {
        "description": (
            "Quasi Sequence Order (QSO) — combines amino acid composition with SOCN "
            "values to encode sequence-order effects. First 20 values are weighted "
            "AAComp; remaining lag values are SOCN contributions. N + 20 features."
        ),
        "module": "sequence_order",
        "category": "Sequence Order",
        "abbreviation": "QSO",
        "parameters": {
            "sequence": "str — protein sequence",
            "lag": "int (default=30) — maximum gap; output features = lag + 20",
            "weight": "float (default=0.1) — weighting factor for SOCN terms",
            "distance_matrix": "str (default='schneider-wrede') — 'schneider-wrede' or 'grantham'",
        },
        "output_shape": "1 x (lag + 20)  (default 1 x 50)",
        "speed": "Moderate",
        "reference": "[1] Chou, K-C. (2000). Biochem. Biophys. Res. Commun., 278, 477-483.",
    },
    "quasi_sequence_order_all": {
        "description": (
            "QSO using both distance matrices concatenated. Output = (lag + 20) × 2 features."
        ),
        "module": "sequence_order",
        "category": "Sequence Order",
        "abbreviation": "QSO_all",
        "parameters": {
            "sequence": "str — protein sequence",
            "lag": "int (default=30) — maximum gap",
            "weight": "float (default=0.1) — weighting factor",
        },
        "output_shape": "1 x ((lag + 20) × 2)  (default 1 x 100)",
        "speed": "Moderate",
        "reference": "[1] Chou, K-C. (2000). Biochem. Biophys. Res. Commun., 278, 477-483.",
    },
}

#speed rating display symbols for pretty-printing
_SPEED_SYMBOL: dict[str, str] = {
    "Fast": "✅ Fast",
    "Moderate": "🟡 Moderate",
    "Slow": "🔴 Slow",
}


def get_descriptor_info(name: str) -> dict:
    """
    Return metadata for a protpy descriptor by name.

    Retrieves a copy of the descriptor's registry entry, including its
    description, module, category, abbreviation, parameters, output shape,
    speed rating and primary reference.

    Parameters
    ==========
    :name: str
        Name of the descriptor function (e.g. "amino_acid_composition").
        Case-insensitive; leading/trailing whitespace is stripped.

    Returns
    =======
    :info: dict
        Dictionary with keys:
        - ``name``            : canonical function name
        - ``description``     : human-readable description
        - ``module``          : source module within protpy
        - ``category``        : descriptor family (e.g. "Composition")
        - ``abbreviation``    : short abbreviation used in the literature
        - ``parameters``      : dict of {param_name: description_str}
        - ``output_shape``    : string describing output DataFrame shape
        - ``speed``           : speed rating (Fast / Moderate / Slow)
        - ``reference``       : primary citation

    Raises
    ======
    KeyError
        If ``name`` does not match any registered descriptor.

    Examples
    ========
    >>> import protpy
    >>> info = protpy.get_descriptor_info("amino_acid_composition")
    >>> print(info["output_shape"])
    1 x 20
    """
    #normalise input name
    name = name.strip().lower()

    if name not in _DESCRIPTOR_REGISTRY:
        #suggest close matches to aid discoverability
        from difflib import get_close_matches
        suggestions = get_close_matches(name, _DESCRIPTOR_REGISTRY.keys(), n=3, cutoff=0.5)
        hint = f"  Did you mean: {suggestions}?" if suggestions else ""
        raise KeyError(f"Descriptor '{name}' not found in the protpy registry.{hint}")

    #return a copy so callers cannot mutate the registry
    info = {"name": name}
    info.update(_DESCRIPTOR_REGISTRY[name])
    return info


def list_descriptors(category: str | None = None) -> list[str]:
    """
    Return a list of all registered descriptor names, optionally filtered by category.

    Parameters
    ==========
    :category: str | None (default=None)
        If provided, only descriptors belonging to this category are returned.
        Available categories: "Composition", "Autocorrelation", "Conjoint Triad",
        "CTD", "Sequence Order". Case-insensitive.

    Returns
    =======
    :names: list[str]
        Sorted list of matching descriptor names.

    Examples
    ========
    >>> import protpy
    >>> protpy.list_descriptors("CTD")
    ['ctd_', 'ctd_composition', 'ctd_distribution', 'ctd_transition']
    """
    if category is None:
        return sorted(_DESCRIPTOR_REGISTRY.keys())

    #case-insensitive category filter
    cat_lower = category.strip().lower()
    return sorted(
        name for name, meta in _DESCRIPTOR_REGISTRY.items()
        if meta["category"].lower() == cat_lower
    )


def print_descriptor_info(name: str) -> None:
    """
    Pretty-print the metadata for a descriptor to stdout.

    Parameters
    ==========
    :name: str
        Descriptor function name (case-insensitive).

    Examples
    ========
    >>> import protpy
    >>> protpy.print_descriptor_info("gravy")
    """
    info = get_descriptor_info(name)

    #speed display with symbol
    speed_label = _SPEED_SYMBOL.get(info["speed"], info["speed"])

    print(f"\n{'='*60}")
    print(f"  Descriptor : {info['name']}  [{info['abbreviation']}]")
    print(f"  Category   : {info['category']}  |  Module: protpy.{info['module']}")
    print(f"  Speed      : {speed_label}")
    print(f"  Output     : {info['output_shape']}")
    print(f"{'='*60}")
    print(f"\n  {info['description']}\n")
    print("  Parameters:")
    for param, desc in info["parameters"].items():
        print(f"    - {param}: {desc}")
    print(f"\n  Reference: {info['reference']}")
    print(f"{'='*60}\n")
