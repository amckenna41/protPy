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
            "region": "str | None (default=None) — None, 'nterm', 'cterm' or 'split' (NT_/MID_/CT_ blocks; split AAC, Kumar et al. 2006)",
            "n": "int (default=25) — terminal block length",
        },
        "output_shape": "1 x 20",
        "speed": "Fast",
        "reference": "[1] Gromiha, M. M. (2010). Protein Bioinformatics. Elsevier.",
        "doi": None,
        "in_default_set": True,
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
            "region": "str | None (default=None) — None, 'nterm', 'cterm' or 'split' (NT_/MID_/CT_ blocks)",
            "n": "int (default=25) — terminal block length",
        },
        "output_shape": "1 x 400",
        "speed": "Fast",
        "reference": "[1] Gromiha, M. M. (2010). Protein Bioinformatics. Elsevier.",
        "doi": None,
        "in_default_set": True,
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
        "doi": None,
        "in_default_set": True,
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
        "doi": "10.1016/0022-2836(82)90515-0",
        "in_default_set": True,
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
        "doi": "10.1093/nar/22.15.3174",
        "in_default_set": True,
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
        "doi": "10.1093/protein/4.2.155",
        "in_default_set": True,
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
        "doi": "10.1002/elps.11501401163",
        "in_default_set": True,
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
        "doi": "10.1385/1-59259-890-0:571",
        "in_default_set": True,
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
        "doi": None,
        "in_default_set": True,
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
        "doi": None,
        "in_default_set": True,
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
        "doi": "10.1021/bi00699a001",
        "in_default_set": True,
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
        "doi": None,
        "in_default_set": True,
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
            "alphabet": "str | None (default=None) — 'ifeature5' (GAAC/GDPC/GTPC), 'cp13', 'cp14', 'cp19' or 'cp20'",
            "k": "int (default=1) — group k-mer length, 1-3",
        },
        "output_shape": "1 x alphabet_size",
        "speed": "Fast",
        "reference": "[1] Gromiha, M. M. (2010). Protein Bioinformatics. Elsevier.",
        "doi": None,
        "in_default_set": True,
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
        "doi": "10.1093/nar/gks1067",
        "in_default_set": True,
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
        "doi": None,
        "in_default_set": True,
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
        "doi": "10.1093/oxfordjournals.jbchem.a133168",
        "in_default_set": True,
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
        "doi": "10.1385/1-59259-890-0:571",
        "in_default_set": True,
    },
    "boman_index": {
        "description": (
            "Boman Index (potential protein interaction index) — negated sum of residue "
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
        "doi": "10.1046/j.1365-2796.2003.01228.x",
        "in_default_set": True,
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
            "charge_threshold": "int (default=0) — max charged residues allowed in an APR window",
        },
        "output_shape": "1 x 2  [AggregProneRegions, AggregProneFraction]",
        "speed": "Fast",
        "reference": "[8] Kyte & Doolittle (1982). J. Mol. Biol., 157(1), 105-132.",
        "doi": "10.1016/0022-2836(82)90515-0",
        "in_default_set": True,
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
        "output_shape": "1 x 2  [HydrophobicMoment_Mean, HydrophobicMoment_Max]",
        "speed": "Fast",
        "reference": "[19] Eisenberg et al. (1982). Nature, 299, 371-374.",
        "doi": "10.1038/299371a0",
        "in_default_set": True,
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
        "doi": None,
        "in_default_set": True,
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
        "doi": "10.1002/prot.1035",
        "in_default_set": True,
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
        "doi": "10.1093/bioinformatics/bth466",
        "in_default_set": True,
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
        "doi": "10.1023/A:1023247831238",
        "in_default_set": True,
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
        "doi": "10.1186/1471-2105-8-300",
        "in_default_set": True,
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
        "doi": "10.1186/1471-2105-8-300",
        "in_default_set": True,
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
            "gap": "int (default=0) — >0 gives KSCTriad (columns suffixed _gap{g})",
        },
        "output_shape": "1 x 343",
        "speed": "Fast",
        "reference": "[1] Shen et al. (2007). PNAS, 104(11), 4337-4341.",
        "doi": "10.1073/pnas.0607879104",
        "in_default_set": True,
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
        "doi": "10.1073/pnas.92.19.8700",
        "in_default_set": True,
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
        "doi": "10.1073/pnas.92.19.8700",
        "in_default_set": True,
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
        "doi": "10.1073/pnas.92.19.8700",
        "in_default_set": True,
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
            "property_set": "str (default='default') — 'default' (7 properties) or 'ifeature13' (13 properties, 1 x 273)",
        },
        "output_shape": "1 x 21  (single property) / 1 x 147  (all 7 properties)",
        "speed": "Fast (single) / Moderate (all)",
        "reference": "[1] Dubchak et al. (1995). PNAS, 92, 8700-8704.",
        "doi": "10.1073/pnas.92.19.8700",
        "in_default_set": True,
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
        "doi": "10.1006/bbrc.2000.3815",
        "in_default_set": True,
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
        "doi": "10.1006/bbrc.2000.3815",
        "in_default_set": True,
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
        "doi": "10.1006/bbrc.2000.3815",
        "in_default_set": True,
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
        "doi": "10.1006/bbrc.2000.3815",
        "in_default_set": True,
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
        "doi": "10.1006/bbrc.2000.3815",
        "in_default_set": True,
    },
    "cksaap": {
        "description": "Composition of k-spaced amino acid pairs: CKSAAP(XY,g) = N_g(XY) / (L-g-1) x 100 for g = 0..gap.",
        "module": "composition",
        "category": "Composition",
        "abbreviation": "CKSAAP",
        "parameters": {
            "sequence": "str — protein sequence",
            "gap": "int (default=5) — maximum spacing",
        },
        "output_shape": "1 x 400·(gap+1)  (default 1 x 2400)",
        "speed": "Moderate",
        "reference": "Chen, Kurgan & Ruan (2007). BMC Struct. Biol., 7, 25.",
        "doi": "10.1186/1472-6807-7-25",
        "descriptor_sets": ["all"],
        "in_default_set": False,
    },
    "cksaagp": {
        "description": "Composition of k-spaced amino acid group pairs over the iFeature 5-group alphabet.",
        "module": "composition",
        "category": "Composition",
        "abbreviation": "CKSAAGP",
        "parameters": {
            "sequence": "str — protein sequence",
            "gap": "int (default=5) — maximum spacing",
        },
        "output_shape": "1 x 25·(gap+1)  (default 1 x 150)",
        "speed": "Fast",
        "reference": "Chen et al. (2018). Bioinformatics, 34(14), 2499-2502 (iFeature).",
        "doi": "10.1093/bioinformatics/bty140",
        "descriptor_sets": ["extended", "all"],
        "in_default_set": False,
    },
    "dde": {
        "description": "Dipeptide deviation from expected mean: (Dc - Tm) / sqrt(Tv), with Tm from codon counts.",
        "module": "composition",
        "category": "Composition",
        "abbreviation": "DDE",
        "parameters": {
            "sequence": "str — protein sequence (length >= 2)",
        },
        "output_shape": "1 x 400",
        "speed": "Fast",
        "reference": "Saravanan & Gautham (2015). OMICS, 19(10), 648-658.",
        "doi": "10.1089/omi.2015.0095",
        "descriptor_sets": ["all"],
        "in_default_set": False,
    },
    "asdc": {
        "description": "Adaptive skip dipeptide composition: composition of all ordered pairs (s_i, s_j), i < j.",
        "module": "composition",
        "category": "Composition",
        "abbreviation": "ASDC",
        "parameters": {
            "sequence": "str — protein sequence",
        },
        "output_shape": "1 x 400",
        "speed": "Fast",
        "reference": "Wei et al. (2018). Bioinformatics, 34(23), 4007-4016 (ACPred-FL).",
        "doi": "10.1093/bioinformatics/bty451",
        "descriptor_sets": ["all"],
        "in_default_set": False,
    },
    "distance_pair": {
        "description": "PseAAC of distance-pairs over a reduced cluster profile (cp13/cp14/cp19/cp20).",
        "module": "composition",
        "category": "Composition",
        "abbreviation": "DP",
        "parameters": {
            "sequence": "str — protein sequence",
            "rng": "int (default=3) — maximum pair distance",
            "alphabet": "str (default='cp14') — cluster profile",
        },
        "output_shape": "1 x (G + G²·rng)  (default 1 x 602)",
        "speed": "Fast",
        "reference": "Liu et al. (2014). PLoS ONE, 9(9), e106691 (iDNA-Prot|dis).",
        "doi": "10.1371/journal.pone.0106691",
        "descriptor_sets": ["all"],
        "in_default_set": False,
    },
    "natural_vector": {
        "description": "Natural vector: count, mean position and normalised second central moment of each amino acid.",
        "module": "composition",
        "category": "Composition",
        "abbreviation": "NV",
        "parameters": {
            "sequence": "str — protein sequence",
        },
        "output_shape": "1 x 60",
        "speed": "Fast",
        "reference": "Yu et al. (2013). J. Theor. Biol., 318, 197-204.",
        "doi": "10.1016/j.jtbi.2012.11.005",
        "descriptor_sets": ["extended", "all"],
        "in_default_set": False,
    },
    "atomic_composition": {
        "description": "Percentage of C, H, N, O and S atoms over the residues' free amino-acid formulas.",
        "module": "composition",
        "category": "Composition",
        "abbreviation": "ATC",
        "parameters": {
            "sequence": "str — protein sequence",
            "subtract_water": "bool (default=False) — remove (L-1) H2O",
        },
        "output_shape": "1 x 5",
        "speed": "Fast",
        "reference": "Pande et al. (2023). J. Comput. Biol., 30(2), 204-222 (Pfeature).",
        "doi": "10.1089/cmb.2022.0241",
        "descriptor_sets": ["extended", "all"],
        "in_default_set": False,
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
        - ``doi``             : DOI of the primary citation (None for textbooks / sources without one)
        - ``in_default_set``  : whether the descriptor is part of the classic ``all_descriptors`` set

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
        "CTD", "Sequence Order", "Scales", "Physicochemical", "Complexity",
        "Signal", "Structural". Case-insensitive.

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
    if info.get("doi"):
        print(f"  DOI      : https://doi.org/{info['doi']}")
    print(f"{'='*60}\n")
