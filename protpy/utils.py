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
            "scale": "str (default='kyte_doolittle') — single-value scale or AAindex code (column GRAVY_{scale})",
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
            "pka_scale": "str (default='protpy') — pKa set: 'protpy', 'emboss', 'bjellqvist' or 'ipc_protein'",
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
            "pka_scale": "str (default='protpy') — pKa set: 'protpy', 'emboss', 'bjellqvist' or 'ipc_protein'",
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
            "absorbance": "bool (default=False) — add Abs01_Reduced / Abs01_Oxidized (A280 of 1 g/L)",
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
            "angle": "float | str (default=100) — rotation angle in degrees, or 'helix' (100) / 'sheet' (160)",
            "scale": "str | None (default=None) — hydrophobicity scale; None keeps protpy's default table",
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
    "z_scales_5": {
        "description": "Mean of the Sandberg 5 z-scales over the sequence (optionally per domain).",
        "module": "scales",
        "category": "Scales",
        "abbreviation": "Z5",
        "parameters": {
            "sequence": "str — protein sequence",
            "domains": "int (default=1) — number of near-equal segments",
        },
        "output_shape": "1 x 5·domains  (default 1 x 5)",
        "speed": "Fast",
        "reference": "Sandberg et al. (1998). J. Med. Chem., 41(14), 2481-2491.",
        "doi": "10.1021/jm9700575",
        "descriptor_sets": ["extended", "all"],
        "in_default_set": False,
    },
    "vhse": {
        "description": "Mean of the VHSE (hydrophobic, steric, electronic) scales over the sequence (optionally per domain).",
        "module": "scales",
        "category": "Scales",
        "abbreviation": "VHSE",
        "parameters": {
            "sequence": "str — protein sequence",
            "domains": "int (default=1) — number of near-equal segments",
        },
        "output_shape": "1 x 8·domains  (default 1 x 8)",
        "speed": "Fast",
        "reference": "Mei et al. (2005). Biopolymers, 80(6), 775-786.",
        "doi": "10.1002/bip.20296",
        "descriptor_sets": ["extended", "all"],
        "in_default_set": False,
    },
    "atchley_factors": {
        "description": "Mean of the Atchley factors I-V over the sequence (optionally per domain).",
        "module": "scales",
        "category": "Scales",
        "abbreviation": "AF",
        "parameters": {
            "sequence": "str — protein sequence",
            "domains": "int (default=1) — number of near-equal segments",
        },
        "output_shape": "1 x 5·domains  (default 1 x 5)",
        "speed": "Fast",
        "reference": "Atchley et al. (2005). PNAS, 102(18), 6395-6400.",
        "doi": "10.1073/pnas.0408677102",
        "descriptor_sets": ["extended", "all"],
        "in_default_set": False,
    },
    "t_scales": {
        "description": "Mean of the topological T-scales over the sequence (optionally per domain).",
        "module": "scales",
        "category": "Scales",
        "abbreviation": "T",
        "parameters": {
            "sequence": "str — protein sequence",
            "domains": "int (default=1) — number of near-equal segments",
        },
        "output_shape": "1 x 5·domains  (default 1 x 5)",
        "speed": "Fast",
        "reference": "Tian, Zhou & Li (2007). J. Mol. Struct., 830, 106-115.",
        "doi": "10.1016/j.molstruc.2006.07.004",
        "descriptor_sets": ["extended", "all"],
        "in_default_set": False,
    },
    "st_scales": {
        "description": "Mean of the structural topology ST-scales over the sequence (optionally per domain).",
        "module": "scales",
        "category": "Scales",
        "abbreviation": "ST",
        "parameters": {
            "sequence": "str — protein sequence",
            "domains": "int (default=1) — number of near-equal segments",
        },
        "output_shape": "1 x 8·domains  (default 1 x 8)",
        "speed": "Fast",
        "reference": "Yang et al. (2010). Amino Acids, 38(3), 805-816.",
        "doi": "10.1007/s00726-009-0287-y",
        "descriptor_sets": ["extended", "all"],
        "in_default_set": False,
    },
    "protfp": {
        "description": "Mean of the ProtFP (PCA, 8 components) over the sequence (optionally per domain).",
        "module": "scales",
        "category": "Scales",
        "abbreviation": "ProtFP",
        "parameters": {
            "sequence": "str — protein sequence",
            "domains": "int (default=1) — number of near-equal segments",
        },
        "output_shape": "1 x 8·domains  (default 1 x 8)",
        "speed": "Fast",
        "reference": "van Westen et al. (2013). J. Cheminform., 5, 41.",
        "doi": "10.1186/1758-2946-5-41",
        "descriptor_sets": ["extended", "all"],
        "in_default_set": False,
    },
    "fasgai": {
        "description": "Mean of the FASGAI factor scales over the sequence (optionally per domain).",
        "module": "scales",
        "category": "Scales",
        "abbreviation": "FASGAI",
        "parameters": {
            "sequence": "str — protein sequence",
            "domains": "int (default=1) — number of near-equal segments",
        },
        "output_shape": "1 x 6·domains  (default 1 x 6)",
        "speed": "Fast",
        "reference": "Liang & Li (2007). QSAR Comb. Sci., 26(6), 754-763.",
        "doi": "10.1002/qsar.200630145",
        "descriptor_sets": ["extended", "all"],
        "in_default_set": False,
    },
    "blosum_indices": {
        "description": "Mean of the BLOSUM62-derived indices over the sequence (optionally per domain).",
        "module": "scales",
        "category": "Scales",
        "abbreviation": "BLOSUM",
        "parameters": {
            "sequence": "str — protein sequence",
            "domains": "int (default=1) — number of near-equal segments",
        },
        "output_shape": "1 x 10·domains  (default 1 x 10)",
        "speed": "Fast",
        "reference": "Georgiev (2009). J. Comput. Biol., 16(5), 703-723.",
        "doi": "10.1089/cmb.2008.0173",
        "descriptor_sets": ["extended", "all"],
        "in_default_set": False,
    },
    "ms_whim": {
        "description": "Mean of the MS-WHIM scores over the sequence (optionally per domain).",
        "module": "scales",
        "category": "Scales",
        "abbreviation": "MSWHIM",
        "parameters": {
            "sequence": "str — protein sequence",
            "domains": "int (default=1) — number of near-equal segments",
        },
        "output_shape": "1 x 3·domains  (default 1 x 3)",
        "speed": "Fast",
        "reference": "Zaliani & Gancia (1999). J. Chem. Inf. Comput. Sci., 39(3), 525-533.",
        "doi": "10.1021/ci980211b",
        "descriptor_sets": ["extended", "all"],
        "in_default_set": False,
    },
    "auto_cross_covariance": {
        "description": "Auto/cross-covariance of an amino-acid scale: sum z_j(i) z_k(i+l) / (L-l).",
        "module": "scales",
        "category": "Scales",
        "abbreviation": "ACC",
        "parameters": {
            "sequence": "str — protein sequence",
            "scale": "str | list[str] (default='z_scales_5') — factor scale or AAindex code(s)",
            "lag": "int (default=5) — maximum lag",
            "mode": "str (default='acc') — 'acc', 'ac' or 'cc'",
        },
        "output_shape": "1 x d²·lag  (default 1 x 125)",
        "speed": "Fast",
        "reference": "Wold et al. (1993). Anal. Chim. Acta, 277(2), 239-253.",
        "doi": "10.1016/0003-2670(93)80437-P",
        "descriptor_sets": ["all"],
        "in_default_set": False,
    },
    "physicochemical_distance_transformation": {
        "description": "Physicochemical distance transformation: mean squared difference of standardised AAindex properties at lag l.",
        "module": "scales",
        "category": "Scales",
        "abbreviation": "PDT",
        "parameters": {
            "sequence": "str — protein sequence",
            "lag": "int (default=3) — maximum lag",
            "properties": "list[str] | None — AAindex codes (default: 8 autocorrelation properties)",
        },
        "output_shape": "1 x len(properties)·lag  (default 1 x 24)",
        "speed": "Fast",
        "reference": "Liu et al. (2012). PLoS ONE, 7(9), e46633.",
        "doi": "10.1371/journal.pone.0046633",
        "descriptor_sets": ["extended", "all"],
        "in_default_set": False,
    },
    "charge_fractions": {
        "description": "Fraction of charged residues (FCR), positive/negative fractions and charge asymmetry sigma = NCPR²/FCR.",
        "module": "physicochemical",
        "category": "Physicochemical",
        "abbreviation": "FCR",
        "parameters": {
            "sequence": "str — protein sequence",
        },
        "output_shape": "1 x 4",
        "speed": "Fast",
        "reference": "Das & Pappu (2013). PNAS, 110(33), 13392-13397.",
        "doi": "10.1073/pnas.1304749110",
        "descriptor_sets": ["extended", "all"],
        "in_default_set": False,
    },
    "ncpr": {
        "description": "Net charge per residue: (n+ - n-) / L with K/R positive and D/E negative.",
        "module": "physicochemical",
        "category": "Physicochemical",
        "abbreviation": "NCPR",
        "parameters": {
            "sequence": "str — protein sequence",
        },
        "output_shape": "1 x 1",
        "speed": "Fast",
        "reference": "Mao et al. (2010). PNAS, 107(18), 8183-8188.",
        "doi": "10.1073/pnas.0911107107",
        "descriptor_sets": ["extended", "all"],
        "in_default_set": False,
    },
    "kappa": {
        "description": "Das-Pappu kappa: normalised blob charge-asymmetry deviation measuring segregation of opposite charges (0 mixed, 1 segregated).",
        "module": "physicochemical",
        "category": "Physicochemical",
        "abbreviation": "Kappa",
        "parameters": {
            "sequence": "str — protein sequence",
            "blob_sizes": "tuple[int] (default=(5, 6)) — blob sizes averaged",
        },
        "output_shape": "1 x (1 + len(blob_sizes))  (default 1 x 3)",
        "speed": "Moderate",
        "reference": "Das & Pappu (2013). PNAS, 110(33), 13392-13397.",
        "doi": "10.1073/pnas.1304749110",
        "descriptor_sets": ["extended", "all"],
        "in_default_set": False,
    },
    "sequence_charge_decoration": {
        "description": "Sequence charge decoration: (1/L) sum_{m>n} q_m q_n sqrt(m-n).",
        "module": "physicochemical",
        "category": "Physicochemical",
        "abbreviation": "SCD",
        "parameters": {
            "sequence": "str — protein sequence",
        },
        "output_shape": "1 x 1",
        "speed": "Fast",
        "reference": "Sawle & Ghosh (2015). J. Chem. Phys., 143(8), 085101.",
        "doi": "10.1063/1.4929391",
        "descriptor_sets": ["extended", "all"],
        "in_default_set": False,
    },
    "charge_density": {
        "description": "Net charge at the given pH divided by molecular weight (charges per Da).",
        "module": "physicochemical",
        "category": "Physicochemical",
        "abbreviation": "ChargeDensity",
        "parameters": {
            "sequence": "str — protein sequence",
            "ph": "float (default=7.4) — pH",
        },
        "output_shape": "1 x 1",
        "speed": "Fast",
        "reference": "Müller et al. (2017). Bioinformatics, 33(17), 2753-2755 (modlAMP).",
        "doi": "10.1093/bioinformatics/btx285",
        "descriptor_sets": ["extended", "all"],
        "in_default_set": False,
    },
    "hydrophobicity_profile": {
        "description": "Summary statistics of a sliding-window hydrophobicity profile (ProtScale-style edge weighting).",
        "module": "physicochemical",
        "category": "Physicochemical",
        "abbreviation": "HP",
        "parameters": {
            "sequence": "str — protein sequence",
            "scale": "str (default='kyte_doolittle') — single-value scale or AAindex code",
            "window": "int (default=9) — window length",
            "edge": "float (default=1.0) — relative end weight",
            "threshold": "float | None (default=None) — adds fraction of windows above threshold",
        },
        "output_shape": "1 x 4 (1 x 5 with threshold)",
        "speed": "Fast",
        "reference": "Kyte & Doolittle (1982). J. Mol. Biol., 157(1), 105-132.",
        "doi": "10.1016/0022-2836(82)90515-0",
        "descriptor_sets": ["extended", "all"],
        "in_default_set": False,
    },
    "hydrophobic_ratio": {
        "description": "APD hydrophobic ratio: percentage of L, I, V, F, M, A, C, W residues.",
        "module": "physicochemical",
        "category": "Physicochemical",
        "abbreviation": "HydrophobicRatio",
        "parameters": {
            "sequence": "str — protein sequence",
        },
        "output_shape": "1 x 1",
        "speed": "Fast",
        "reference": "Wang (2020). Protein Sci., 29(1), 8-18 (APD); Wang, Li & Wang (2016). Nucleic Acids Res., 44, D1087 (APD3).",
        "doi": "10.1002/pro.3702",
        "descriptor_sets": ["extended", "all"],
        "in_default_set": False,
    },
    "half_life": {
        "description": "N-end rule half-life (hours) in mammalian, yeast and E. coli systems from the N-terminal residue.",
        "module": "physicochemical",
        "category": "Physicochemical",
        "abbreviation": "HalfLife",
        "parameters": {
            "sequence": "str — protein sequence",
        },
        "output_shape": "1 x 3",
        "speed": "Fast",
        "reference": "Bachmair, Finley & Varshavsky (1986). Science, 234(4773), 179-186.",
        "doi": "10.1126/science.3018930",
        "descriptor_sets": ["extended", "all"],
        "in_default_set": False,
    },
    "ivywrel_index": {
        "description": "Fraction of I, V, Y, W, R, E, L residues (thermophilic adaptation index).",
        "module": "physicochemical",
        "category": "Physicochemical",
        "abbreviation": "IVYWREL",
        "parameters": {
            "sequence": "str — protein sequence",
        },
        "output_shape": "1 x 1",
        "speed": "Fast",
        "reference": "Zeldovich, Berezovsky & Shakhnovich (2007). PLoS Comput. Biol., 3(1), e5.",
        "doi": "10.1371/journal.pcbi.0030005",
        "descriptor_sets": ["extended", "all"],
        "in_default_set": False,
    },
    "cvp_bias": {
        "description": "Charged versus polar bias: %(D+E+K+R) - %(N+Q+S+T).",
        "module": "physicochemical",
        "category": "Physicochemical",
        "abbreviation": "CvP",
        "parameters": {
            "sequence": "str — protein sequence",
        },
        "output_shape": "1 x 1",
        "speed": "Fast",
        "reference": "Suhre & Claverie (2003). J. Biol. Chem., 278(19), 17198-17202.",
        "doi": "10.1074/jbc.M301327200",
        "descriptor_sets": ["extended", "all"],
        "in_default_set": False,
    },
    "flexibility": {
        "description": "Summary statistics of the Vihinen average-B-value flexibility profile (window 9, end weight 0.25).",
        "module": "physicochemical",
        "category": "Physicochemical",
        "abbreviation": "Flex",
        "parameters": {
            "sequence": "str — protein sequence",
            "window": "int (default=9) — window length",
        },
        "output_shape": "1 x 4",
        "speed": "Fast",
        "reference": "Vihinen, Torkkila & Riikonen (1994). Proteins, 19(2), 141-149.",
        "doi": "10.1002/prot.340190207",
        "descriptor_sets": ["extended", "all"],
        "in_default_set": False,
    },
    "wootton_federhen_complexity": {
        "description": "Wootton-Federhen local complexity (K1, K2) over sliding windows and the SEG-trigger low-complexity fraction.",
        "module": "complexity",
        "category": "Complexity",
        "abbreviation": "WF",
        "parameters": {
            "sequence": "str — protein sequence",
            "window": "int (default=12) — window length",
            "k2_trigger": "float (default=2.2) — K2 trigger (bits)",
        },
        "output_shape": "1 x 4",
        "speed": "Fast",
        "reference": "Wootton & Federhen (1993). Comput. Chem., 17(2), 149-163.",
        "doi": "10.1016/0097-8485(93)85006-X",
        "descriptor_sets": ["extended", "all"],
        "in_default_set": False,
    },
    "linguistic_complexity": {
        "description": "Linguistic complexity: distinct subwords observed / maximum possible (ratio of sums over word lengths).",
        "module": "complexity",
        "category": "Complexity",
        "abbreviation": "LC",
        "parameters": {
            "sequence": "str — protein sequence",
            "max_k": "int | None (default=None) — longest word length (default L)",
        },
        "output_shape": "1 x 1",
        "speed": "Fast",
        "reference": "Troyanskaya et al. (2002). Bioinformatics, 18(5), 679-688.",
        "doi": "10.1093/bioinformatics/18.5.679",
        "descriptor_sets": ["extended", "all"],
        "in_default_set": False,
    },
    "lempel_ziv_complexity": {
        "description": "Lempel-Ziv (LZ76) complexity c(n) and its normalised value c(n)·log20(n)/n.",
        "module": "complexity",
        "category": "Complexity",
        "abbreviation": "LZ",
        "parameters": {
            "sequence": "str — protein sequence",
        },
        "output_shape": "1 x 2",
        "speed": "Fast",
        "reference": "Lempel & Ziv (1976). IEEE Trans. Inf. Theory, 22(1), 75-81.",
        "doi": "10.1109/TIT.1976.1055501",
        "descriptor_sets": ["extended", "all"],
        "in_default_set": False,
    },
    "fft_descriptors": {
        "description": "Mean power of the FFT spectrum of a mean-centred property signal (default EIIP) in equal-width frequency bins.",
        "module": "signal",
        "category": "Signal",
        "abbreviation": "FFT",
        "parameters": {
            "sequence": "str — protein sequence",
            "scale": "str | dict (default='eiip') — single-value scale, AAindex code or dict",
            "n_bins": "int (default=10) — number of frequency bins",
        },
        "output_shape": "1 x n_bins  (default 1 x 10)",
        "speed": "Fast",
        "reference": "Cosic (1994). IEEE Trans. Biomed. Eng., 41(12), 1101-1114.",
        "doi": "10.1109/10.335859",
        "descriptor_sets": ["extended", "all"],
        "in_default_set": False,
    },
    "dct_descriptors": {
        "description": "First n orthonormal DCT-II coefficients of a mean-centred property signal (default Kyte-Doolittle).",
        "module": "signal",
        "category": "Signal",
        "abbreviation": "DCT",
        "parameters": {
            "sequence": "str — protein sequence",
            "scale": "str | dict (default='kyte_doolittle') — single-value scale, AAindex code or dict",
            "n_coeffs": "int (default=10) — number of coefficients (<= L)",
        },
        "output_shape": "1 x n_coeffs  (default 1 x 10)",
        "speed": "Fast",
        "reference": "Wang et al. (2017). PLoS ONE, 12(9), e0185587.",
        "doi": "10.1371/journal.pone.0185587",
        "descriptor_sets": ["extended", "all"],
        "in_default_set": False,
    },
    "aggrescan": {
        "description": "AGGRESCAN aggregation hot-spot descriptors from the a3v propensity profile (hot spots: >=5 residues above -0.02, no Pro).",
        "module": "structural",
        "category": "Structural",
        "abbreviation": "AGG",
        "parameters": {
            "sequence": "str — protein sequence",
            "window": "int (default=5) — odd window length (AGGRESCAN uses 5, 7, 9, 11)",
        },
        "output_shape": "1 x 12",
        "speed": "Fast",
        "reference": "Conchillo-Solé et al. (2007). BMC Bioinformatics, 8, 65.",
        "doi": "10.1186/1471-2105-8-65",
        "descriptor_sets": ["extended", "all"],
        "in_default_set": False,
    },
    "uversky_charge_hydropathy": {
        "description": "Uversky mean hydropathy <H>, mean net charge <R> and boundary distance HQ = 2.785<H> - <R> - 1.151.",
        "module": "structural",
        "category": "Structural",
        "abbreviation": "UV",
        "parameters": {
            "sequence": "str — protein sequence",
        },
        "output_shape": "1 x 4",
        "speed": "Fast",
        "reference": "Uversky, Gillespie & Fink (2000). Proteins, 41(3), 415-427.",
        "doi": "10.1002/1097-0134(20001115)41:3<415::AID-PROT130>3.0.CO;2-7",
        "descriptor_sets": ["extended", "all"],
        "in_default_set": False,
    },
    "fold_index": {
        "description": "FoldIndex profile summary: windowed 2.785<H> - abs(<R>) - 1.151 and disordered-residue statistics.",
        "module": "structural",
        "category": "Structural",
        "abbreviation": "FI",
        "parameters": {
            "sequence": "str — protein sequence",
            "window": "int (default=51) — window length",
        },
        "output_shape": "1 x 4",
        "speed": "Fast",
        "reference": "Prilusky et al. (2005). Bioinformatics, 21(16), 3435-3438.",
        "doi": "10.1093/bioinformatics/bti537",
        "descriptor_sets": ["extended", "all"],
        "in_default_set": False,
    },
    "disorder_propensity": {
        "description": "Fractions of disorder- and order-promoting residues and mean TOP-IDP value.",
        "module": "structural",
        "category": "Structural",
        "abbreviation": "DP",
        "parameters": {
            "sequence": "str — protein sequence",
        },
        "output_shape": "1 x 3",
        "speed": "Fast",
        "reference": "Campen et al. (2008). Protein Pept. Lett., 15(9), 956-963.",
        "doi": "10.2174/092986608785849164",
        "descriptor_sets": ["extended", "all"],
        "in_default_set": False,
    },
    "transmembrane_features": {
        "description": "Kyte-Doolittle transmembrane-segment heuristics (segment count, max hydropathy, TM fraction, positive-inside bias).",
        "module": "structural",
        "category": "Structural",
        "abbreviation": "TM",
        "parameters": {
            "sequence": "str — protein sequence",
            "window": "int (default=19) — hydropathy window",
            "kd_threshold": "float (default=1.6) — segment threshold",
        },
        "output_shape": "1 x 4",
        "speed": "Fast",
        "reference": "Kyte & Doolittle (1982). J. Mol. Biol., 157(1), 105-132.",
        "doi": "10.1016/0022-2836(82)90515-0",
        "descriptor_sets": ["extended", "all"],
        "in_default_set": False,
    },
    "signal_peptide_features": {
        "description": "N-terminal signal-peptide heuristics: n-region net charge and maximum 8-residue hydropathy.",
        "module": "structural",
        "category": "Structural",
        "abbreviation": "SP",
        "parameters": {
            "sequence": "str — protein sequence",
            "n_term": "int (default=30) — N-terminal region length",
        },
        "output_shape": "1 x 2",
        "speed": "Fast",
        "reference": "von Heijne (1983). Eur. J. Biochem., 133(1), 17-21.",
        "doi": "10.1111/j.1432-1033.1983.tb07424.x",
        "descriptor_sets": ["extended", "all"],
        "in_default_set": False,
    },
    "functional_site_counts": {
        "description": "Counts and per-100-residue densities of PROSITE N-glycosylation, ER-retention and PTS1 motifs and tryptic sites.",
        "module": "structural",
        "category": "Structural",
        "abbreviation": "FS",
        "parameters": {
            "sequence": "str — protein sequence",
        },
        "output_shape": "1 x 8",
        "speed": "Fast",
        "reference": "Sigrist et al. (2013). Nucleic Acids Res., 41(D1), D344-D347 (PROSITE).",
        "doi": "10.1093/nar/gks1067",
        "descriptor_sets": ["extended", "all"],
        "in_default_set": False,
    },
    "plaac": {
        "description": "PLAAC prion-like domain scores: max window LLR, HMM core score, PrD length and presence.",
        "module": "structural",
        "category": "Structural",
        "abbreviation": "PLAAC",
        "parameters": {
            "sequence": "str — protein sequence",
            "core_length": "int (default=60) — minimal core length",
        },
        "output_shape": "1 x 4",
        "speed": "Fast",
        "reference": "Lancaster et al. (2014). Bioinformatics, 30(17), 2501-2502.",
        "doi": "10.1093/bioinformatics/btu310",
        "descriptor_sets": ["extended", "all"],
        "in_default_set": False,
    },
    "antigenicity_surface": {
        "description": "Profile statistics of the Parker hydrophilicity and Janin surface scales over 7-residue windows.",
        "module": "structural",
        "category": "Structural",
        "abbreviation": "Antigenicity",
        "parameters": {
            "sequence": "str — protein sequence",
            "window": "int (default=7) — window length",
        },
        "output_shape": "1 x 8",
        "speed": "Fast",
        "reference": "Parker, Guo & Hodges (1986). Biochemistry, 25(19), 5425-5432.",
        "doi": "10.1021/bi00367a013",
        "descriptor_sets": ["extended", "all"],
        "in_default_set": False,
    },
    "cysteine_features": {
        "description": "Cysteine count, fraction, parity, maximum disulfides, CxxC motifs and mean spacing.",
        "module": "structural",
        "category": "Structural",
        "abbreviation": "Cys",
        "parameters": {
            "sequence": "str — protein sequence",
        },
        "output_shape": "1 x 6",
        "speed": "Fast",
        "reference": "Fomenko & Gladyshev (2003). Biochemistry, 42(38), 11214-11225.",
        "doi": "10.1021/bi034459s",
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
