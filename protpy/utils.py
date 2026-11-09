################################################################################
#############                      Utilities                      ##############
################################################################################

from __future__ import annotations

"""
Utility and auxiliary functions for the protpy package. Includes descriptor
metadata lookup, shared helpers, and other convenience tools.
"""

#registry of the descriptors with metadata. Deprecated functions are not registered.
#descriptor_sets: the all_descriptors sets the descriptor belongs to ([] = opt-in via include);
#alias_of: the family function a permanent literature alias calls.
_DESCRIPTOR_REGISTRY: dict[str, dict] = {
    'amino_acid_composition': {
        'description': 'Percentage of each of the 20 canonical amino acids in the sequence. AA_Comp(t) = N(t) / L x 100.',
        'module': 'composition',
        'category': 'Composition',
        'abbreviation': 'AAComp',
        "parameters": {
            'sequence': 'str — protein sequence',
            'region': "str | None (default=None) — None, 'nterm', 'cterm' or 'split' (NT_/MID_/CT_ blocks; split AAC, Kumar et al. 2006)",
            'n': 'int (default=25) — terminal block length',
        },
        'output_shape': '1 x 20',
        'speed': 'Fast',
        'reference': '[1] Gromiha, M. M. (2010). Protein Bioinformatics. Elsevier.',
        'doi': None,
        'descriptor_sets': ['classic', 'extended', 'all'],
        'in_default_set': True,
    },
    'dipeptide_composition': {
        'description': "Percentage of each of the 400 dipeptides: DPComp(s,t) = N(s,t) / (L - 1) x 100; normalize='dde' gives the dipeptide deviation from expected mean.",
        'module': 'composition',
        'category': 'Composition',
        'abbreviation': 'DPComp',
        "parameters": {
            'sequence': 'str — protein sequence',
            'region': "str | None (default=None) — None, 'nterm', 'cterm' or 'split' (NT_/MID_/CT_ blocks)",
            'n': 'int (default=25) — terminal block length',
            'normalize': "str | None (default=None) — None (percentages) or 'dde' (dipeptide deviation from expected mean)",
        },
        'output_shape': '1 x 400',
        'speed': 'Fast',
        'reference': '[1] Gromiha, M. M. (2010). Protein Bioinformatics. Elsevier.',
        'doi': None,
        'descriptor_sets': ['classic', 'extended', 'all'],
        'in_default_set': True,
    },
    'gravy': {
        'description': 'Grand Average of Hydropathy (GRAVY) — sum of Kyte-Doolittle hydropathy values divided by sequence length. Positive = hydrophobic, negative = hydrophilic.',
        'module': 'composition',
        'category': 'Composition',
        'abbreviation': 'GRAVY',
        "parameters": {
            'sequence': 'str — protein sequence',
            'scale': "str (default='kyte_doolittle') — single-value scale or AAindex code (column GRAVY_{scale})",
        },
        'output_shape': '1 x 1',
        'speed': 'Fast',
        'reference': '[8] Kyte & Doolittle (1982). J. Mol. Biol., 157(1), 105-132.',
        'doi': '10.1016/0022-2836(82)90515-0',
        'descriptor_sets': ['classic', 'extended', 'all'],
        'in_default_set': True,
    },
    'aromaticity': {
        'description': 'Fraction of aromatic amino acids (Phe, Trp, Tyr) in the sequence.',
        'module': 'composition',
        'category': 'Composition',
        'abbreviation': 'Aromaticity',
        "parameters": {
            'sequence': 'str — protein sequence',
        },
        'output_shape': '1 x 1',
        'speed': 'Fast',
        'reference': '[9] Lobry & Gautier (1994). Nucleic Acids Research, 22(15), 3174-3180.',
        'doi': '10.1093/nar/22.15.3174',
        'descriptor_sets': ['classic', 'extended', 'all'],
        'in_default_set': True,
    },
    'instability_index': {
        'description': 'Instability Index — 10/L times the sum of the Guruprasad DIWV values of all dipeptides (Biopython table). Values < 40 predict in vivo stability; >= 40 predicts instability.',
        'module': 'composition',
        'category': 'Composition',
        'abbreviation': 'II',
        "parameters": {
            'sequence': 'str — protein sequence',
        },
        'output_shape': '1 x 1',
        'speed': 'Fast',
        'reference': '[10] Guruprasad et al. (1990). Protein Engineering, 4(2), 155-161.',
        'doi': '10.1093/protein/4.2.155',
        'descriptor_sets': ['classic', 'extended', 'all'],
        'in_default_set': True,
    },
    'isoelectric_point': {
        'description': 'Theoretical isoelectric point (pI) — the pH at which the net charge is zero, by bisection; Bjellqvist pKa values by default (as ExPASy ProtParam / Biopython).',
        'module': 'composition',
        'category': 'Composition',
        'abbreviation': 'pI',
        "parameters": {
            'sequence': 'str — protein sequence',
            'pka_scale': "str (default='bjellqvist') — pKa set: 'bjellqvist', 'emboss', 'ipc_protein' or 'protpy' (v1.0-v1.5 default)",
        },
        'output_shape': '1 x 1',
        'speed': 'Fast',
        'reference': '[11] Bjellqvist et al. (1993). Electrophoresis, 14(1), 1023-1031.',
        'doi': '10.1002/elps.11501401163',
        'descriptor_sets': ['classic', 'extended', 'all'],
        'in_default_set': True,
    },
    'molecular_weight': {
        'description': 'Molecular Weight (MW) — sum of average residue masses minus one water molecule per peptide bond (18.015 Da).',
        'module': 'composition',
        'category': 'Composition',
        'abbreviation': 'MW',
        "parameters": {
            'sequence': 'str — protein sequence',
        },
        'output_shape': '1 x 1',
        'speed': 'Fast',
        'reference': '[12] Gasteiger et al. (2005). The Proteomics Protocols Handbook.',
        'doi': '10.1385/1-59259-890-0:571',
        'descriptor_sets': ['classic', 'extended', 'all'],
        'in_default_set': True,
    },
    'aliphatic_index': {
        'description': 'Aliphatic Index (AI) — relative volume of aliphatic side chains (Ala, Val, Ile, Leu). Higher values indicate thermostability. AI = Xala + 2.9·Xval + 3.9·(Xile + Xleu).',
        'module': 'composition',
        'category': 'Composition',
        'abbreviation': 'AI',
        "parameters": {
            'sequence': 'str — protein sequence',
        },
        'output_shape': '1 x 1',
        'speed': 'Fast',
        'reference': '[16] Ikai (1980). J. Biochemistry, 88(6), 1895-1898.',
        'doi': '10.1093/oxfordjournals.jbchem.a133168',
        'descriptor_sets': ['classic', 'extended', 'all'],
        'in_default_set': True,
    },
    'extinction_coefficient': {
        'description': 'Molar extinction coefficient at 280 nm from W, Y and cystines (reduced and oxidized forms; all_descriptors keeps only ExtCoeff_Oxidized).',
        'module': 'composition',
        'category': 'Composition',
        'abbreviation': 'ExtCoeff',
        "parameters": {
            'sequence': 'str — protein sequence',
            'absorbance': 'bool (default=False) — add Abs01_Reduced / Abs01_Oxidized (A280 of 1 g/L)',
        },
        'output_shape': '1 x 2  [ExtCoeff_Reduced, ExtCoeff_Oxidized]',
        'speed': 'Fast',
        'reference': '[17] Gasteiger et al. (2005). The Proteomics Protocols Handbook.',
        'doi': '10.1385/1-59259-890-0:571',
        'descriptor_sets': ['classic', 'extended', 'all'],
        'in_default_set': True,
    },
    'boman_index': {
        'description': "Boman (protein-interaction) index: minus the mean Radzicka-Wolfenden solubility value; equals -gravy(scale='RADA880101').",
        'module': 'composition',
        'category': 'Composition',
        'abbreviation': 'BI',
        "parameters": {
            'sequence': 'str — protein sequence',
        },
        'output_shape': '1 x 1',
        'speed': 'Fast',
        'reference': '[18] Boman (2003). J. Internal Medicine, 254(3), 197-215.',
        'doi': '10.1046/j.1365-2796.2003.01228.x',
        'alias_of': 'gravy',
        'descriptor_sets': ['classic', 'extended', 'all'],
        'in_default_set': True,
    },
    'charge_distribution': {
        'description': "Henderson-Hasselbalch positive, negative and net side-chain charge at a given pH (Bjellqvist pKa); normalize='length' or 'mass' gives charge per residue or per Da.",
        'module': 'composition',
        'category': 'Composition',
        'abbreviation': 'ChargeDist',
        "parameters": {
            'sequence': 'str — protein sequence',
            'ph': 'float (default=7.4) — pH for charge calculation',
            'pka_scale': "str (default='bjellqvist') — pKa set: 'bjellqvist', 'emboss', 'ipc_protein' or 'protpy'",
            'normalize': "str (default='none') — 'none', 'length' (per residue) or 'mass' (per Da)",
        },
        'output_shape': '1 x 3',
        'speed': 'Fast',
        'reference': '[13] Cameselle-Teijeiro (1979). Biochemical Education, 7(3), 69-70.',
        'doi': None,
        'descriptor_sets': ['classic', 'extended', 'all'],
        'in_default_set': True,
    },
    'hydrophobic_moment': {
        'description': 'Eisenberg hydrophobic moment (μH) — measure of amphiphilicity calculated as the vector sum of residue hydrophobicities projected at a given helical angle. Averaged over all windows of the given size.',
        'module': 'composition',
        'category': 'Composition',
        'abbreviation': 'μH',
        "parameters": {
            'sequence': 'str — protein sequence',
            'window': 'int (default=11) — helical window size',
            'angle': "float | str (default=100) — rotation angle in degrees, or 'helix' (100) / 'sheet' (160)",
            'scale': "str | None (default=None) — hydrophobicity scale; None keeps protpy's default table",
        },
        'output_shape': '1 x 2  [HydrophobicMoment_Mean, HydrophobicMoment_Max]',
        'speed': 'Fast',
        'reference': '[19] Eisenberg et al. (1982). Nature, 299, 371-374.',
        'doi': '10.1038/299371a0',
        'descriptor_sets': ['classic', 'extended', 'all'],
        'in_default_set': True,
    },
    'shannon_entropy': {
        'description': 'Shannon entropy of the amino acid composition — measure of sequence complexity/diversity. H = -sum(p_i * log2(p_i)) over all 20 amino acids.',
        'module': 'complexity',
        'category': 'Complexity',
        'abbreviation': 'ShannonH',
        "parameters": {
            'sequence': 'str — protein sequence',
        },
        'output_shape': '1 x 1',
        'speed': 'Fast',
        'reference': '[1] Gromiha, M. M. (2010). Protein Bioinformatics. Elsevier.',
        'doi': None,
        'descriptor_sets': ['classic', 'extended', 'all'],
        'in_default_set': True,
    },
    'motif_composition': {
        'description': 'Counts of regex sequence motifs: 8 default motifs, plus the PROSITE N-glycosylation, ER-retention and PTS1 patterns and tryptic sites by name; optional per-100-residue densities.',
        'module': 'composition',
        'category': 'Composition',
        'abbreviation': 'MotifComp',
        "parameters": {
            'sequence': 'str — protein sequence',
            'motifs': 'dict | list | None (default=None) — {name: regex}, or a list of built-in motif names / regex patterns; None uses the 8 defaults',
            'per_100': 'bool (default=False) — also return counts per 100 residues',
        },
        'output_shape': '1 x M  (default 1 x 8; 1 x 2M with per_100)',
        'speed': 'Fast',
        'reference': '[15] PROSITE: https://prosite.expasy.org/',
        'doi': '10.1093/nar/gks1067',
        'descriptor_sets': ['classic', 'extended', 'all'],
        'in_default_set': True,
    },
    'pseudo_amino_acid_composition': {
        'description': 'Pseudo Amino Acid Composition (PAAComp) — augments standard AAComp with sequence-order correlation factors derived from physicochemical properties. First 20 values are weighted AAComp; remaining lamda values are correlation factors.',
        'module': 'composition',
        'category': 'Composition',
        'abbreviation': 'PAAComp',
        "parameters": {
            'sequence': 'str — protein sequence',
            'lamda': 'int (default=30) — rank of correlation (output = 20 + lamda features)',
            'weight': 'float (default=0.05) — weighting factor for correlation terms',
            'properties': 'list[str] (default=[]) — AAIndex1 accession codes; empty = use built-in hydrophobicity/hydrophilicity/mass',
        },
        'output_shape': '1 x (20 + lamda)  (default 1 x 50)',
        'speed': 'Moderate',
        'reference': '[4] Chou, K-C. (2001). PROTEINS, 43, 246-255.',
        'doi': '10.1002/prot.1035',
        'descriptor_sets': ['classic', 'extended', 'all'],
        'in_default_set': True,
    },
    'amphiphilic_pseudo_amino_acid_composition': {
        'description': 'Amphiphilic Pseudo Amino Acid Composition (APAAComp) — extends PAAComp with hydrophobicity and hydrophilicity distribution patterns. First 20 values are weighted AAComp; remaining 2·lamda values encode amphiphilic correlations.',
        'module': 'composition',
        'category': 'Composition',
        'abbreviation': 'APAAComp',
        "parameters": {
            'sequence': 'str — protein sequence',
            'lamda': 'int (default=30) — rank of correlation (output = 20 + 2·lamda features)',
            'weight': 'float (default=0.5) — weighting factor',
            'properties': 'list[dict] (default=[hydrophobicity_, hydrophilicity_]) — property value dicts',
        },
        'output_shape': '1 x (20 + 2·lamda)  (default 1 x 80)',
        'speed': 'Moderate',
        'reference': '[5] Chou, K-C. (2005). Bioinformatics, 21(1), 10-19.',
        'doi': '10.1093/bioinformatics/bth466',
        'descriptor_sets': ['classic', 'extended', 'all'],
        'in_default_set': True,
    },
    'geary_autocorrelation': {
        'description': 'Geary Autocorrelation (GAuto) — uses the squared difference of property values for residue pairs separated by lag, instead of vector products. Generates lag × n_properties features.',
        'module': 'autocorrelation',
        'category': 'Autocorrelation',
        'abbreviation': 'GAuto',
        "parameters": {
            'sequence': 'str — protein sequence',
            'lag': 'int (default=30) — max lag between residue pairs',
            'properties': 'list[str] (default=8 AAIndex codes) — AAIndex1 accession numbers',
            'normalize': 'bool (default=True) — standardise property values (z-score) before calculation',
        },
        'output_shape': '1 x (lag × n_properties)  (default 1 x 240)',
        'speed': 'Moderate',
        'reference': '[2] Ong et al. (2007). BMC Bioinformatics, 8, 300.',
        'doi': '10.1186/1471-2105-8-300',
        'alias_of': 'autocorrelation',
        'descriptor_sets': ['classic', 'extended', 'all'],
        'in_default_set': True,
    },
    'conjoint_triad': {
        'description': 'Conjoint triad: frequencies of the 343 triads of a 7-class alphabet, normalised (f - min) / max as in Shen et al. 2007 (normalize=False for raw counts); gap > 0 gives KSCTriad.',
        'module': 'conjoint_triad',
        'category': 'Conjoint Triad',
        'abbreviation': 'CTriad',
        "parameters": {
            'sequence': 'str — protein sequence',
            'gap': 'int (default=0) — >0 gives KSCTriad (columns suffixed _gap{g})',
            'normalize': 'bool (default=True) — (f - min) / max normalisation (Shen 2007); False for raw counts',
        },
        'output_shape': '1 x 343·(gap+1)  (default 1 x 343)',
        'speed': 'Fast',
        'reference': '[1] Shen et al. (2007). PNAS, 104(11), 4337-4341.',
        'doi': '10.1073/pnas.0607879104',
        'descriptor_sets': ['classic', 'extended', 'all'],
        'in_default_set': True,
    },
    'ctd': {
        'description': 'Composition, Transition and Distribution (CTD) of 3-class property alphabets for 7 properties (or the 13 of iFeature); parts= selects C, T and/or D.',
        'module': 'ctd',
        'category': 'CTD',
        'abbreviation': 'CTD',
        "parameters": {
            'sequence': 'str — protein sequence',
            'property': "str (default='hydrophobicity') — CTD property (used when all_ctd=False)",
            'all_ctd': 'bool (default=True) — if True, compute over all 7 properties',
            'property_set': "str (default='default') — 'default' (7 properties) or 'ifeature13' (13 properties, 1 x 273)",
            'parts': "str (default='CTD') — any combination of 'C', 'T' and 'D'",
        },
        'output_shape': '1 x 147 (7 properties; 1 x 273 for ifeature13)',
        'speed': 'Fast (single) / Moderate (all)',
        'reference': '[1] Dubchak et al. (1995). PNAS, 92, 8700-8704.',
        'doi': '10.1073/pnas.92.19.8700',
        'descriptor_sets': ['classic', 'extended', 'all'],
        'in_default_set': True,
    },
    'quasi_sequence_order': {
        'description': "Quasi-sequence-order: amino acid fractions and coupling numbers weighted by 1 + w sum(tau); the 20 + lag values sum to 1. distance_matrix='both' concatenates Schneider-Wrede and Grantham (as in all_descriptors).",
        'module': 'sequence_order',
        'category': 'Sequence Order',
        'abbreviation': 'QSO',
        "parameters": {
            'sequence': 'str — protein sequence',
            'lag': 'int (default=30) — maximum gap; output features = lag + 20',
            'weight': 'float (default=0.1) — weighting factor for SOCN terms',
            'distance_matrix': "str (default='schneider-wrede') — 'schneider-wrede', 'grantham' or 'both' (used by all_descriptors)",
        },
        'output_shape': "1 x (lag + 20)  (default 1 x 50; 1 x 100 with distance_matrix='both')",
        'speed': 'Moderate',
        'reference': '[1] Chou, K-C. (2000). Biochem. Biophys. Res. Commun., 278, 477-483.',
        'doi': '10.1006/bbrc.2000.3815',
        'descriptor_sets': ['classic', 'extended', 'all'],
        'in_default_set': True,
    },
    'natural_vector': {
        'description': 'Natural vector: count, mean position and normalised second central moment of each amino acid.',
        'module': 'composition',
        'category': 'Composition',
        'abbreviation': 'NV',
        "parameters": {
            'sequence': 'str — protein sequence',
            'normalize': 'bool (default=False) — divide n_k and mu_k by the sequence length',
        },
        'output_shape': '1 x 60',
        'speed': 'Fast',
        'reference': 'Yu et al. (2013). J. Theor. Biol., 318, 197-204.',
        'doi': '10.1016/j.jtbi.2012.11.005',
        'descriptor_sets': ['extended', 'all'],
        'in_default_set': False,
    },
    'atomic_composition': {
        'description': "Percentage of C, H, N, O and S atoms over the residues' free amino-acid formulas.",
        'module': 'composition',
        'category': 'Composition',
        'abbreviation': 'ATC',
        "parameters": {
            'sequence': 'str — protein sequence',
            'subtract_water': 'bool (default=False) — remove (L-1) H2O',
        },
        'output_shape': '1 x 5',
        'speed': 'Fast',
        'reference': 'Pande et al. (2023). J. Comput. Biol., 30(2), 204-222 (Pfeature).',
        'doi': '10.1089/cmb.2022.0241',
        'descriptor_sets': ['extended', 'all'],
        'in_default_set': False,
    },
    'cksaagp': {
        'description': 'Composition of k-spaced amino acid group pairs over the iFeature 5-group alphabet.',
        'module': 'composition',
        'category': 'Composition',
        'abbreviation': 'CKSAAGP',
        "parameters": {
            'sequence': 'str — protein sequence',
            'gap': 'int (default=5) — maximum spacing',
        },
        'output_shape': '1 x 25·(gap+1)  (default 1 x 150)',
        'speed': 'Fast',
        'reference': 'Chen et al. (2018). Bioinformatics, 34(14), 2499-2502 (iFeature).',
        'doi': '10.1093/bioinformatics/bty140',
        'alias_of': 'cksaap',
        'descriptor_sets': ['extended', 'all'],
        'in_default_set': False,
    },
    'charge_fractions': {
        'description': 'Fraction of charged residues (FCR), of positive and negative residues, net charge per residue (NCPR) and charge asymmetry (sigma).',
        'module': 'physicochemical',
        'category': 'Physicochemical',
        'abbreviation': 'FCR',
        "parameters": {
            'sequence': 'str — protein sequence',
        },
        'output_shape': '1 x 5',
        'speed': 'Fast',
        'reference': 'Das & Pappu (2013). PNAS, 110(33), 13392-13397.',
        'doi': '10.1073/pnas.1304749110',
        'descriptor_sets': ['extended', 'all'],
        'in_default_set': False,
    },
    'kappa': {
        'description': 'Das & Pappu kappa: segregation of oppositely charged residues (0 well mixed, 1 segregated), averaged over blob sizes; NaN without charged residues.',
        'module': 'physicochemical',
        'category': 'Physicochemical',
        'abbreviation': 'Kappa',
        "parameters": {
            'sequence': 'str — protein sequence',
            'blob_sizes': 'tuple (default=(5, 6)) — blob sizes averaged into Kappa',
        },
        'output_shape': '1 x 1',
        'speed': 'Moderate',
        'reference': 'Das & Pappu (2013). PNAS, 110(33), 13392-13397.',
        'doi': '10.1073/pnas.1304749110',
        'descriptor_sets': ['extended', 'all'],
        'in_default_set': False,
    },
    'sequence_charge_decoration': {
        'description': 'Sequence charge decoration: (1/L) sum_{m>n} q_m q_n sqrt(m-n).',
        'module': 'physicochemical',
        'category': 'Physicochemical',
        'abbreviation': 'SCD',
        "parameters": {
            'sequence': 'str — protein sequence',
        },
        'output_shape': '1 x 1',
        'speed': 'Fast',
        'reference': 'Sawle & Ghosh (2015). J. Chem. Phys., 143(8), 085101.',
        'doi': '10.1063/1.4929391',
        'descriptor_sets': ['extended', 'all'],
        'in_default_set': False,
    },
    'hydrophobicity_profile': {
        'description': 'Sliding-window profile of any single-value scale (hydropathy, flexibility, Parker, Janin, TOP-IDP, ...) summarised by mean, max, min and std; optional threshold fraction and segment (e.g. transmembrane) features.',
        'module': 'physicochemical',
        'category': 'Physicochemical',
        'abbreviation': 'HP',
        "parameters": {
            'sequence': 'str — protein sequence',
            'scale': "str (default='kyte_doolittle') — any single-value protpy scale or AAindex code",
            'window': 'int (default=9) — window length',
            'edge': 'float (default=1.0) — relative end weight',
            'threshold': 'float | None (default=None) — adds fraction of windows above threshold',
            'segments': 'bool (default=False) — add segment count, fraction and positive-inside bias (needs threshold)',
        },
        'output_shape': '1 x 4 (5 with threshold, 8 with segments)',
        'speed': 'Fast',
        'reference': 'Kyte & Doolittle (1982). J. Mol. Biol., 157(1), 105-132.',
        'doi': '10.1016/0022-2836(82)90515-0',
        'descriptor_sets': ['extended', 'all'],
        'in_default_set': False,
    },
    'thermostability_indices': {
        'description': 'IVYWREL index and charged-versus-polar (CvP) bias, composition-based thermostability indicators.',
        'module': 'physicochemical',
        'category': 'Physicochemical',
        'abbreviation': 'Thermo',
        "parameters": {
            'sequence': 'str — protein sequence',
        },
        'output_shape': '1 x 2',
        'speed': 'Fast',
        'reference': 'Zeldovich, Berezovsky & Shakhnovich (2007). PLoS Comput. Biol., 3(1), e5.',
        'doi': '10.1371/journal.pcbi.0030005',
        'descriptor_sets': ['extended', 'all'],
        'in_default_set': False,
    },
    'wootton_federhen_complexity': {
        'description': 'Wootton-Federhen local complexity: mean and minimum windowed entropy K2 and the fraction of residues in low-complexity windows (trigger step only, not masked regions).',
        'module': 'complexity',
        'category': 'Complexity',
        'abbreviation': 'WF',
        "parameters": {
            'sequence': 'str — protein sequence',
            'window': 'int (default=12) — window length',
            'k2_trigger': 'float (default=2.2) — K2 trigger (bits)',
        },
        'output_shape': '1 x 3',
        'speed': 'Fast',
        'reference': 'Wootton & Federhen (1993). Comput. Chem., 17(2), 149-163.',
        'doi': '10.1016/0097-8485(93)85006-X',
        'descriptor_sets': ['extended', 'all'],
        'in_default_set': False,
    },
    'lempel_ziv_complexity': {
        'description': 'Normalised Lempel-Ziv (LZ76) complexity c(n) log_20(n) / n.',
        'module': 'complexity',
        'category': 'Complexity',
        'abbreviation': 'LZ',
        "parameters": {
            'sequence': 'str — protein sequence',
        },
        'output_shape': '1 x 1',
        'speed': 'Fast',
        'reference': 'Lempel & Ziv (1976). IEEE Trans. Inf. Theory, 22(1), 75-81.',
        'doi': '10.1109/TIT.1976.1055501',
        'descriptor_sets': ['extended', 'all'],
        'in_default_set': False,
    },
    'spectral_descriptors': {
        'description': 'Binned FFT power spectrum (default EIIP) or DCT coefficients 1..n (default Kyte-Doolittle) of a property signal.',
        'module': 'spectral',
        'category': 'Spectral',
        'abbreviation': 'Spectral',
        "parameters": {
            'sequence': 'str — protein sequence',
            'transform': "str (default='fft') — 'fft' (binned power spectrum) or 'dct' (coefficients 1..n)",
            'scale': "str | dict | None (default=None) — single-value scale; None uses 'eiip' (FFT) or 'kyte_doolittle' (DCT)",
            'n': 'int (default=10) — number of FFT bins or DCT coefficients',
        },
        'output_shape': '1 x n  (default 1 x 10)',
        'speed': 'Fast',
        'reference': 'Cosic (1994). IEEE Trans. Biomed. Eng., 41(12), 1101-1114.',
        'doi': '10.1109/10.335859',
        'descriptor_sets': ['extended', 'all'],
        'in_default_set': False,
    },
    'aggrescan': {
        'description': 'AGGRESCAN aggregation hot-spot descriptors from the a3v propensity profile (hot spots: >=5 residues above -0.02, no Pro).',
        'module': 'structural',
        'category': 'Structural',
        'abbreviation': 'AGG',
        "parameters": {
            'sequence': 'str — protein sequence',
            'window': "int | str (default=5) — odd window length, or 'psico' for the server's length rule (5/7/9/11)",
            'full': 'bool (default=False) — return all 12 AGGRESCAN outputs instead of the 5 normalised ones',
        },
        'output_shape': '1 x 5 (1 x 12 with full=True)',
        'speed': 'Fast',
        'reference': 'Conchillo-Solé et al. (2007). BMC Bioinformatics, 8, 65.',
        'doi': '10.1186/1471-2105-8-65',
        'descriptor_sets': ['extended', 'all'],
        'in_default_set': False,
    },
    'charge_hydropathy': {
        'description': 'Uversky charge-hydropathy score HQ = 2.785<H> - abs(<R>) - 1.151 of the whole sequence, or its FoldIndex sliding-window profile (window=51).',
        'module': 'structural',
        'category': 'Structural',
        'abbreviation': 'CH',
        "parameters": {
            'sequence': 'str — protein sequence',
            'window': 'int | None (default=None) — None scores the whole sequence; e.g. 51 gives the FoldIndex profile',
        },
        'output_shape': '1 x 1 (window=None) or 1 x 4',
        'speed': 'Fast',
        'reference': 'Uversky, Gillespie & Fink (2000). Proteins, 41(3), 415-427.',
        'doi': '10.1002/1097-0134(20001115)41:3<415::AID-PROT130>3.0.CO;2-7',
        'descriptor_sets': ['extended', 'all'],
        'in_default_set': False,
    },
    'plaac': {
        'description': 'PLAAC prion-like domain scores: max window LLR, HMM core score, PrD length and presence.',
        'module': 'structural',
        'category': 'Structural',
        'abbreviation': 'PLAAC',
        "parameters": {
            'sequence': 'str — protein sequence',
            'core_length': 'int (default=60) — minimal core length',
            'fill_value': 'float (default=nan) — value of PLAAC_LLR / PLAAC_COREscore where undefined',
        },
        'output_shape': '1 x 4',
        'speed': 'Fast',
        'reference': 'Lancaster et al. (2014). Bioinformatics, 30(17), 2501-2502.',
        'doi': '10.1093/bioinformatics/btu310',
        'descriptor_sets': ['extended', 'all'],
        'in_default_set': False,
    },
    'cksaap': {
        'description': 'Composition of k-spaced amino acid pairs: CKSAAP(XY,g) = N_g(XY) / (L-g-1) x 100 for g = 0..gap; alphabet= groups the residues (CKSAAGP, distance pair).',
        'module': 'composition',
        'category': 'Composition',
        'abbreviation': 'CKSAAP',
        "parameters": {
            'sequence': 'str — protein sequence',
            'gap': 'int (default=5) — maximum spacing',
            'alphabet': 'int | str | list | None (default=None) — reduced alphabet (see kmer_composition)',
        },
        'output_shape': '1 x 400·(gap+1)  (default 1 x 2400), or G²·(gap+1) with an alphabet',
        'speed': 'Moderate',
        'reference': 'Chen, Kurgan & Ruan (2007). BMC Struct. Biol., 7, 25.',
        'doi': '10.1186/1472-6807-7-25',
        'descriptor_sets': ['all'],
        'in_default_set': False,
    },
    'asdc': {
        'description': 'Adaptive skip dipeptide composition: composition of all ordered pairs (s_i, s_j), i < j.',
        'module': 'composition',
        'category': 'Composition',
        'abbreviation': 'ASDC',
        "parameters": {
            'sequence': 'str — protein sequence',
        },
        'output_shape': '1 x 400',
        'speed': 'Fast',
        'reference': 'Wei et al. (2018). Bioinformatics, 34(23), 4007-4016 (ACPred-FL).',
        'doi': '10.1093/bioinformatics/bty451',
        'descriptor_sets': ['all'],
        'in_default_set': False,
    },
    'distance_pair': {
        'description': 'PseAAC of distance-pairs over a reduced cluster profile (cp13/cp14/cp19/cp20).',
        'module': 'composition',
        'category': 'Composition',
        'abbreviation': 'DP',
        "parameters": {
            'sequence': 'str — protein sequence',
            'rng': 'int (default=3) — maximum pair distance',
            'alphabet': "str (default='cp14') — cluster profile",
        },
        'output_shape': '1 x (G + G²·rng)  (default 1 x 602)',
        'speed': 'Fast',
        'reference': 'Liu et al. (2014). PLoS ONE, 9(9), e106691 (iDNA-Prot|dis).',
        'doi': '10.1371/journal.pone.0106691',
        'alias_of': 'cksaap',
        'descriptor_sets': ['all'],
        'in_default_set': False,
    },
    'auto_cross_covariance': {
        'description': 'Auto/cross-covariance of an amino-acid scale: sum z_j(i) z_k(i+l) / (L-l).',
        'module': 'scales',
        'category': 'Scales',
        'abbreviation': 'ACC',
        "parameters": {
            'sequence': 'str — protein sequence',
            'scale': "str | list[str] (default='z_scales_5') — factor scale or AAindex code(s)",
            'lag': 'int (default=5) — maximum lag',
            'mode': "str (default='acc') — 'acc', 'ac' or 'cc'",
        },
        'output_shape': '1 x d²·lag  (default 1 x 125)',
        'speed': 'Fast',
        'reference': 'Wold et al. (1993). Anal. Chim. Acta, 277(2), 239-253.',
        'doi': '10.1016/0003-2670(93)80437-P',
        'descriptor_sets': ['all'],
        'in_default_set': False,
    },
    'moreaubroto_autocorrelation': {
        'description': 'MoreauBroto Autocorrelation (MBAuto) — topological descriptor based on the product of property values for residue pairs separated by lag positions. Generates lag × n_properties features.',
        'module': 'autocorrelation',
        'category': 'Autocorrelation',
        'abbreviation': 'MBAuto',
        "parameters": {
            'sequence': 'str — protein sequence',
            'lag': 'int (default=30) — max lag between residue pairs',
            'properties': 'list[str] (default=8 AAIndex codes) — AAIndex1 accession numbers',
            'normalize': 'bool (default=True) — standardise property values (z-score) before calculation',
        },
        'output_shape': '1 x (lag × n_properties)  (default 1 x 240)',
        'speed': 'Moderate',
        'reference': '[1] Hollas (2003). J. Math. Chem., 33(2), 91-101.',
        'doi': '10.1023/A:1023247831238',
        'alias_of': 'autocorrelation',
        'descriptor_sets': ['all'],
        'in_default_set': False,
    },
    'moran_autocorrelation': {
        'description': 'Moran Autocorrelation (MAuto) — uses property deviations from the mean rather than raw values to compute autocorrelation. Generates lag × n_properties features.',
        'module': 'autocorrelation',
        'category': 'Autocorrelation',
        'abbreviation': 'MAuto',
        "parameters": {
            'sequence': 'str — protein sequence',
            'lag': 'int (default=30) — max lag between residue pairs',
            'properties': 'list[str] (default=8 AAIndex codes) — AAIndex1 accession numbers',
            'normalize': 'bool (default=True) — standardise property values (z-score) before calculation',
        },
        'output_shape': '1 x (lag × n_properties)  (default 1 x 240)',
        'speed': 'Moderate',
        'reference': '[2] Ong et al. (2007). BMC Bioinformatics, 8, 300.',
        'doi': '10.1186/1471-2105-8-300',
        'alias_of': 'autocorrelation',
        'descriptor_sets': ['all'],
        'in_default_set': False,
    },
    'sequence_order_coupling_number': {
        'description': 'Sequence order coupling numbers tau_d = sum of squared Schneider-Wrede or Grantham distances of residues d apart; mostly a length measure unless normalize=True.',
        'module': 'sequence_order',
        'category': 'Sequence Order',
        'abbreviation': 'SOCN',
        "parameters": {
            'sequence': 'str — protein sequence',
            'lag': 'int (default=30) — maximum gap; output features = lag',
            'distance_matrix': "str (default='schneider-wrede') — 'schneider-wrede', 'grantham' or 'both'",
            'normalize': 'bool (default=False) — divide tau_d by L - d',
        },
        'output_shape': "1 x lag  (default 1 x 30; 1 x 60 with distance_matrix='both')",
        'speed': 'Moderate',
        'reference': '[1] Chou, K-C. (2000). Biochem. Biophys. Res. Commun., 278, 477-483.',
        'doi': '10.1006/bbrc.2000.3815',
        'descriptor_sets': ['all'],
        'in_default_set': False,
    },
    'tripeptide_composition': {
        'description': "Percentage of each of the 8000 tripeptides; alias of kmer_composition(k=3). Opt-in only (include=['tripeptide_composition']).",
        'module': 'composition',
        'category': 'Composition',
        'abbreviation': 'TPComp',
        "parameters": {
            'sequence': 'str — protein sequence',
        },
        'output_shape': '1 x 8000',
        'speed': 'Moderate',
        'reference': '[1] Gromiha, M. M. (2010). Protein Bioinformatics. Elsevier.',
        'doi': None,
        'alias_of': 'kmer_composition',
        'descriptor_sets': [],
        'in_default_set': False,
    },
    'kmer_composition': {
        'description': 'Percentage of each overlapping k-mer (k = 1-3) of the 20 amino acids, or of a reduced alphabet (size 2/3/4/6, ifeature5, cp13/14/19/20 or a custom grouping).',
        'module': 'composition',
        'category': 'Composition',
        'abbreviation': 'kmer',
        "parameters": {
            'sequence': 'str — protein sequence',
            'k': 'int (default=2) — k-mer length, 1-3',
            'alphabet': "int | str | list | None (default=None) — size 2/3/4/6, 'ifeature5', 'cp13'/'cp14'/'cp19'/'cp20' or a list of residue groups",
        },
        'output_shape': '1 x 20^k (default 1 x 400), or 1 x G^k with an alphabet',
        'speed': 'Fast (k≤3) / Slow (k=4)',
        'reference': '[1] Gromiha, M. M. (2010). Protein Bioinformatics. Elsevier.',
        'doi': None,
        'descriptor_sets': [],
        'in_default_set': False,
    },
    'dde': {
        'description': 'Dipeptide deviation from expected mean: (Dc - Tm) / sqrt(Tv), with Tm from codon counts.',
        'module': 'composition',
        'category': 'Composition',
        'abbreviation': 'DDE',
        "parameters": {
            'sequence': 'str — protein sequence (length >= 2)',
        },
        'output_shape': '1 x 400',
        'speed': 'Fast',
        'reference': 'Saravanan & Gautham (2015). OMICS, 19(10), 648-658.',
        'doi': '10.1089/omi.2015.0095',
        'alias_of': 'dipeptide_composition',
        'descriptor_sets': [],
        'in_default_set': False,
    },
    'autocorrelation': {
        'description': 'Moreau-Broto, Moran or Geary autocorrelation of standardised AAindex properties for lags 1..lag.',
        'module': 'autocorrelation',
        'category': 'Autocorrelation',
        'abbreviation': 'Auto',
        "parameters": {
            'sequence': 'str — protein sequence',
            'method': "str (default='geary') — 'moreaubroto', 'moran' or 'geary'",
            'lag': 'int (default=30) — max lag between residue pairs',
            'properties': 'list[str] (default=8 AAIndex codes) — AAIndex1 accession numbers',
            'normalize': 'bool (default=True) — standardise property values (z-score) before calculation',
        },
        'output_shape': '1 x (n_properties × lag)  (default 1 x 240)',
        'speed': 'Moderate',
        'reference': '[2] Ong, S. A. K. et al. (2007). BMC Bioinformatics, 8, 300.',
        'doi': '10.1186/1471-2105-8-300',
        'descriptor_sets': [],
        'in_default_set': False,
    },
    'ctd_': {
        'description': 'Combined CTD descriptor (CTD-C + CTD-T + CTD-D) — concatenation of Composition, Transition and Distribution for one or all 7 properties. Single property = 21 features; all properties = 147 features.',
        'module': 'ctd',
        'category': 'CTD',
        'abbreviation': 'CTD',
        "parameters": {
            'sequence': 'str — protein sequence',
            'property': "str (default='hydrophobicity') — CTD property (used when all_ctd=False)",
            'all_ctd': 'bool (default=True) — if True, compute over all 7 properties',
            'property_set': "str (default='default') — 'default' (7 properties) or 'ifeature13' (13 properties, 1 x 273)",
        },
        'output_shape': '1 x 21  (single property) / 1 x 147  (all 7 properties)',
        'speed': 'Fast (single) / Moderate (all)',
        'reference': '[1] Dubchak et al. (1995). PNAS, 92, 8700-8704.',
        'doi': '10.1073/pnas.92.19.8700',
        'alias_of': 'ctd',
        'descriptor_sets': [],
        'in_default_set': False,
    },
    'ctd_composition': {
        'description': 'CTD Composition (CTD-C) — fraction of residues belonging to each of three property-defined classes (class 1, 2 and 3). 3 features per property.',
        'module': 'ctd',
        'category': 'CTD',
        'abbreviation': 'CTD-C',
        "parameters": {
            'sequence': 'str — protein sequence',
            'property': "str (default='hydrophobicity') — one of: hydrophobicity, normalized_vdwv, polarity, charge, secondary_struct, solvent_accessibility, polarizability",
        },
        'output_shape': '1 x 3  (per property)',
        'speed': 'Fast',
        'reference': '[1] Dubchak et al. (1995). PNAS, 92, 8700-8704.',
        'doi': '10.1073/pnas.92.19.8700',
        'alias_of': 'ctd',
        'descriptor_sets': [],
        'in_default_set': False,
    },
    'ctd_transition': {
        'description': 'CTD Transition (CTD-T) — percentage frequency of transitions between different property classes for consecutive residues. 3 features per property.',
        'module': 'ctd',
        'category': 'CTD',
        'abbreviation': 'CTD-T',
        "parameters": {
            'sequence': 'str — protein sequence',
            'property': "str (default='hydrophobicity') — CTD physicochemical property",
        },
        'output_shape': '1 x 3  (per property)',
        'speed': 'Fast',
        'reference': '[1] Dubchak et al. (1995). PNAS, 92, 8700-8704.',
        'doi': '10.1073/pnas.92.19.8700',
        'alias_of': 'ctd',
        'descriptor_sets': [],
        'in_default_set': False,
    },
    'ctd_distribution': {
        'description': 'CTD Distribution (CTD-D) — chain positions where the 1st, 25%, 50%, 75% and 100% of residues belonging to a property class are located. 15 features per property.',
        'module': 'ctd',
        'category': 'CTD',
        'abbreviation': 'CTD-D',
        "parameters": {
            'sequence': 'str — protein sequence',
            'property': "str (default='hydrophobicity') — CTD physicochemical property",
        },
        'output_shape': '1 x 15  (per property)',
        'speed': 'Fast',
        'reference': '[1] Dubchak et al. (1995). PNAS, 92, 8700-8704.',
        'doi': '10.1073/pnas.92.19.8700',
        'alias_of': 'ctd',
        'descriptor_sets': [],
        'in_default_set': False,
    },
    'aa_scale_descriptors': {
        'description': 'Mean (or per-position) values of an amino acid factor scale: z-scales, VHSE, Atchley, T-scales, ST-scales, ProtFP, FASGAI, BLOSUM indices or MS-WHIM.',
        'module': 'scales',
        'category': 'Scales',
        'abbreviation': 'Scale',
        "parameters": {
            'sequence': 'str — protein sequence',
            'scale': "str — factor scale: 'z_scales_5', 'vhse', 'atchley_factors', 't_scales', 'st_scales', 'protfp', 'fasgai', 'blosum_indices' or 'ms_whim'",
            'domains': 'int (default=1) — number of near-equal segments',
            'mode': "str (default='mean') — 'mean' or 'position' (per-residue values, aligned sequences)",
        },
        'output_shape': "1 x d·domains (mode='mean') or 1 x L·d (mode='position')",
        'speed': 'Fast',
        'reference': 'van Westen et al. (2013). J. Cheminform., 5, 41.',
        'doi': '10.1186/1758-2946-5-41',
        'descriptor_sets': [],
        'in_default_set': False,
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
        "Spectral", "Structural". Case-insensitive.

    Returns
    =======
    :names: list[str]
        Sorted list of matching descriptor names.

    Examples
    ========
    >>> import protpy
    >>> protpy.list_descriptors("CTD")
    ['ctd', 'ctd_', 'ctd_composition', 'ctd_distribution', 'ctd_transition']
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
