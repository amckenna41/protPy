################################################################################
#############                   Composition                       ##############
################################################################################

from __future__ import annotations

import re
import itertools
from collections import Counter
import numpy as np
import pandas as pd
import math
from aaindex import aaindex1

"""
References
==========
[1] Gromiha, M. M. (2010). Protein Sequence Analysis. In M. M. Gromiha (Ed.), 
    Protein Bioinformatics (pp. 29–62). Elsevier.
[2] Hua, S. and Sun, Z. (2001) Support vector machine approach for protein
    subcellular localization prediction. Bioinformatics, 17, 721-728.
[3] Grassmann, J., Reczko, M., Suhai, S. and Edler, L. (1999) Protein fold
    class prediction: new methods of statistical classification. Proc Int Conf
    Intell Syst Mol Biol, 106-112.
[4] Kuo-Chen Chou. Prediction of Protein Cellular Attributes Using Pseudo-Amino Acid
    Composition. PROTEINS: Structure, Function, and Genetics, 2001, 43: 246-255.
[5] Kuo-Chen Chou, Using amphiphilic pseudo amino acid composition to predict enzyme 
    subfamily classes, Bioinformatics, Volume 21, Issue 1, January 2005, Pages 10–19, 
    https://doi.org/10.1093/bioinformatics/bth466
[6] http://www.csbio.sjtu.edu.cn/bioinf/PseAAC/ParaValue.htm.
[7] Reczko, M. and Bohr, H. (1994) The DEF data base of sequence based protein
    fold class predictions. Nucleic Acids Res, 22, 3616-3619.
[8] Kyte, J., & Doolittle, R. F. (1982). A simple method for displaying the
    hydropathic character of a protein. Journal of Molecular Biology, 157(1),
    105-132. https://doi.org/10.1016/0022-2836(82)90515-0
[9] Lobry, J. R., & Gautier, C. (1994). Hydrophobicity, expressivity and
    aromaticity are the major trends of amino-acid usage in 999 Escherichia coli
    chromosome-encoded genes. Nucleic Acids Research, 22(15), 3174-3180.
[10] Guruprasad, K., Reddy, B. V. B., & Pandit, M. W. (1990). Correlation
    between stability of a protein and its dipeptide composition: a novel
    approach for predicting in vivo stability of a protein from its primary
    sequence. Protein Engineering, 4(2), 155-161.
[11] Bjellqvist, B., et al. (1993). The focusing positions of polypeptides in
    immobilized pH gradients can be predicted from their amino acid sequences.
    Electrophoresis, 14(1), 1023-1031.
[12] Gasteiger, E., et al. (2005). Protein Identification and Analysis Tools on
    the ExPASy Server. In J. M. Walker (Ed.), The Proteomics Protocols Handbook
    (pp. 571-607). Humana Press.
[13] Cameselle-Teijeiro, J. (1979). The Henderson-Hasselbalch equation. 
    Biochemical Education, 7(3), 69-70.
[14] Chou, P. Y., & Fasman, G. D. (1974). Conformational parameters for amino
    acids in helical, beta-sheet, and random coil regions calculated from
    proteins. Biochemistry, 13(2), 211-222.
[15] Prosite database of protein families and domains:
    https://prosite.expasy.org/
[16] Ikai, A. J. (1980). Thermostability and aliphatic index of globular
    proteins. Journal of Biochemistry, 88(6), 1895-1898.
[17] Gasteiger, E., et al. (2005). Protein Identification and Analysis Tools
    on the ExPASy Server. In J. M. Walker (Ed.), The Proteomics Protocols
    Handbook (pp. 571-607). Humana Press. (extinction coefficient method)
[18] Boman, H. G. (2003). Antibacterial peptides: basic facts and emerging
    concepts. Journal of Internal Medicine, 254(3), 197-215.
[19] Eisenberg, D., Weiss, R. M., & Terwilliger, T. C. (1982). The helical
    hydrophobic moment: a measure of the amphiphilicity of a helix. Nature,
    299, 371-374. https://doi.org/10.1038/299371a0
[20] Chen, K., Kurgan, L. A., & Ruan, J. (2007). Prediction of flexible/rigid
    regions from protein sequences using k-spaced amino acid pairs. BMC
    Structural Biology, 7, 25. https://doi.org/10.1186/1472-6807-7-25
[21] Chen, Z. et al. (2018). iFeature: a Python package and web server for
    features extraction and selection from protein and peptide sequences.
    Bioinformatics, 34(14), 2499-2502. https://doi.org/10.1093/bioinformatics/bty140
[22] Saravanan, V., & Gautham, N. (2015). Harnessing computational biology for
    exact linear B-cell epitope prediction: a novel amino acid composition-based
    feature descriptor. OMICS, 19(10), 648-658. https://doi.org/10.1089/omi.2015.0095
[23] Wei, L. et al. (2018). ACPred-FL: a sequence-based predictor using effective
    feature representation to improve the prediction of anti-cancer peptides.
    Bioinformatics, 34(23), 4007-4016. https://doi.org/10.1093/bioinformatics/bty451
[24] Liu, B. et al. (2014). iDNA-Prot|dis: identifying DNA-binding proteins by
    incorporating amino acid distance-pairs and reduced alphabet profile into
    the general pseudo amino acid composition. PLoS ONE, 9(9), e106691.
    https://doi.org/10.1371/journal.pone.0106691
[25] Cannata, N. et al. (2002). Simplifying amino acid alphabets by means of a
    branch and bound algorithm and substitution matrices. Bioinformatics, 18(8),
    1102-1108. https://doi.org/10.1093/bioinformatics/18.8.1102
[26] Yu, C. et al. (2013). Protein space: a natural method for realizing the
    nature of protein universe. Journal of Theoretical Biology, 318, 197-204.
    https://doi.org/10.1016/j.jtbi.2012.11.005
[27] Pande, A. et al. (2023). Computing wide range of protein/peptide features
    from their sequence and structure. Journal of Computational Biology, 30(2),
    204-222. https://doi.org/10.1089/cmb.2022.0241
[28] Kumar, M., Verma, R., & Raghava, G. P. S. (2006). Prediction of
    mitochondrial proteins using support vector machine and hidden Markov
    model. Journal of Biological Chemistry, 281(9), 5357-5363.
    https://doi.org/10.1074/jbc.M511061200
[29] Lata, S., Sharma, B. K., & Raghava, G. P. S. (2007). Analysis and
    prediction of antibacterial peptides. BMC Bioinformatics, 8, 263.
    https://doi.org/10.1186/1471-2105-8-263
"""
#list of amino acids
from ._constants import amino_acids, _validate_sequence, _validate_int, _validate_weight
from ._helpers import (_load_json_table, _read_json, _gapped_pair_counts, _kmer_counts, _to_df, _property_matrix,
                       _minmax, _scale_range)

#physicochemical values for all amino acids for hydrophobicity, hydrophilicty and residue mass
#values taken from http://www.csbio.sjtu.edu.cn/bioinf/PseAAC/ParaValue.htm
#property values differ from their equivalent in the aaindex therefore hard-coding them below
hydrophobicity_ = {
        "A": 0.62, "C": 0.29, "D": -0.90, "E": -0.74, "F": 1.19, "G": 0.48, 
        "H": -0.40, "I": 1.38, "K": -1.50, "L": 1.06, "M": 0.64, "N": -0.78,
        "P": 0.12, "Q": -0.85, "R": -2.53, "S": -0.18, "T": -0.05, "V": 1.08, 
        "W": 0.81, "Y": 0.26
        }

hydrophilicity_ = {
        "A": -0.5, "C": -1.0, "D": 3.0, "E": 3.0, "F": -2.5,  "G": 0.0, 
        "H": -0.5, "I": -1.8, "K": 3.0, "L": -1.8, "M": -1.3, "N": 0.2, 
        "P": 0.0, "Q": 0.2, "R": 3.0, "S": 0.3, "T": -0.4, "V": -1.5, 
        "W": -3.4, "Y": -2.3
        }

residue_mass_ = {
        "A": 15.0, "C": 47.0, "D": 59.0, "E": 73.0, "F": 91.0, "G": 1.000, 
        "H": 82.0, "I": 57.0, "K": 73.0, "L": 57.0,  "M": 75.0, "N": 58.0, 
        "P": 42.0, "Q": 72.0, "R": 101.0, "S": 31.0, "T": 45.0, "V": 43.0, 
        "W": 130.0, "Y": 107.0
        }

################################# Amino Acid Composition ##################################

def amino_acid_composition(sequence: str, region: str | None = None, n: int = 25) -> pd.DataFrame:
    """
    Calculate Amino Acid Composition (AAComp) of protein sequence. AAComp
    describes the fraction of each amino acid type within a protein sequence,
    and is calculated as:

    AA_Comp(s) = AA(t)/N(s)

    where AA_Comp(s) is the AAComp of protein sequence s, AA(t) is the number
    of amino acid types t (where t = 1,2,..,20) and N(s) is the length of the
    sequences [1].

    Parameters
    ==========
    :sequence: str
        protein sequence.
    :region: str | None (default=None)
        restrict the calculation to part of the sequence: None (whole sequence),
        "nterm" (first ``n`` residues, columns prefixed ``NT_``), "cterm" (last
        ``n`` residues, ``CT_``) or "split" (N-terminal, middle and C-terminal
        blocks concatenated, prefixed ``NT_``, ``MID_``, ``CT_``; split
        composition [28, 29]).
    :n: int (default=25)
        terminal block length used by ``region``.

    Returns
    =======
    :amino_acid_composition_df: pd.DataFrame
        pandas dataframe of AAComp for protein sequence. Dataframe will
        be of the shape 1 x 20, where 20 is the number of features 
        calculated from the descriptor (for the 20 amino acids) and 1 is 
        the input sequence.
    """
    #validate input protein sequence
    sequence = _validate_sequence(sequence)
    if region is not None:
        return _regional(amino_acid_composition, sequence, region, n)

    aa_composition = {}

    #iterate through each amino acid, calculating each composition
    for aa in amino_acids:
        aa_composition[aa] = round(float(sequence.count(aa)) / len(sequence) * 100, 3)

    #transform values and columns to DataFrame
    amino_acid_composition_df = pd.DataFrame([list(aa_composition.values())], 
        columns=list(aa_composition.keys()))

    return amino_acid_composition_df

################################## Dipeptide Composition ##################################

def dipeptide_composition(sequence: str, region: str | None = None, n: int = 25) -> pd.DataFrame:
    """
    Calculate Dipeptide Composition (DPComp) for protein sequence.
    Dipeptide composition is the fraction of each dipeptide type within a
    protein sequence. With dipeptides being of length 2 and there being 20
    canonical amino acids this creates 20^2 different combinations, thus a
    400-Dimensional vector will be produced such that:

    DPComp(s,t) = AA(s,t) / N -1

    where DPComp(s,t) is the dipeptide composition of the protein sequence
    for amino acid type s and t (where s and t = 1,2,..,20), AA(s,t) is the number
    of dipeptides represented by amino acid type s and t and N is the total number
    of dipeptides [1].

    Parameters
    ==========
    :sequence: str
        protein sequence.
    :region: str | None (default=None)
        restrict the calculation to part of the sequence: None (whole sequence),
        "nterm" (first ``n`` residues, columns prefixed ``NT_``), "cterm" (last
        ``n`` residues, ``CT_``) or "split" (N-terminal, middle and C-terminal
        blocks concatenated, prefixed ``NT_``, ``MID_``, ``CT_``; split
        composition [28, 29]).
    :n: int (default=25)
        terminal block length used by ``region``.

    Returns
    =======
    :dipeptide_composition_df: pd.DataFrame
        pandas dataframe of dipeptide composition for protein sequence. Dataframe will
        be of the shape 1 x 400, where 400 is the number of features calculated 
        from the descriptor (20^2 for the 20 canonical amino acids) and 1 is the 
        input sequence.
    """
    #validate input protein sequence
    sequence = _validate_sequence(sequence)
    if region is not None:
        return _regional(dipeptide_composition, sequence, region, n)

    dipep_composition = {}

    #count overlapping dipeptides with a sliding window (str.count misses overlaps, e.g. "AA" in "AAA")
    total_dipeptides = len(sequence) - 1
    counts = Counter(sequence[i:i+2] for i in range(total_dipeptides))

    #iterate through each amino acid, calculating each dipeptide combination's composition
    for i in amino_acids:
        for j in amino_acids:
            dipep = i + j
            dipep_composition[dipep] = round(
                counts[dipep] / total_dipeptides * 100, 2) if total_dipeptides > 0 else 0.0

    #transform values and columns to DataFrame
    dipeptide_composition_df = pd.DataFrame([list(dipep_composition.values())], 
        columns=list(dipep_composition.keys()))

    return dipeptide_composition_df

def _regional(func, sequence: str, region: str, n: int) -> pd.DataFrame:
    """Apply ``func`` to the N-terminal, C-terminal or split (N/middle/C) regions of ``sequence``."""
    n = _validate_int("n", n, 1)
    parts = {"nterm": [("NT_", sequence[:n])], "cterm": [("CT_", sequence[-n:])],
             "split": [("NT_", sequence[:n]), ("MID_", sequence[n:-n]), ("CT_", sequence[-n:])]}
    if region not in parts:
        raise ValueError(f"Invalid region '{region}'. Valid regions are: None, {list(parts)}.")
    #the middle block of a split must be non-empty
    if region == "split" and len(sequence) <= 2 * n:
        raise ValueError(f"region='split' needs a sequence longer than 2n ({2 * n}), got length {len(sequence)}.")
    frames = []
    for prefix, part in parts[region]:
        df = func(part)
        frames.append(df.add_prefix(prefix))
    return pd.concat(frames, axis=1)

################################## Tripeptide Composition #################################

def tripeptide_composition(sequence: str) -> pd.DataFrame:
    """
    Calculate Tripeptide Composition (TPComp) of protein sequence.
    Tripeptide composition is the fraction of each tripeptide type within a
    protein sequence. With tripeptides being of length 3 and there being 20
    canonical amino acids this creates 20^3 different combinations, thus a
    8000-Dimensional vector will be produced such that:

    TPComp(s,t,u) = AA(s,t,u) / N - 1

    where TPComp(s,t,u) is the tripeptide composition of the protein sequence
    for amino acid type s, t and u (where s, t and u = 1,2,..,20), AA(s,t,u) is
    the number of tripeptides represented by amino acid type s and t, and N is
    the total number of tripeptides [1].

    Parameters
    ==========
    :sequence: str
        protein sequence in str form.

    Returns
    =======
    :tripeptide_composition_df: pd.DataFrame
        pandas DataFrame of tripeptide composition for protein sequence. Dataframe will
        be of the shape 1 x 8000, where 8000 is the number of features calculated
        from the descriptor (20^3 for the 20 canonical amino acids) and 1 is the
        input sequence. Values are percentages (0–100), consistent with
        dipeptide_composition.
    """
    #validate input protein sequence
    sequence = _validate_sequence(sequence)

    tripeptide_composition = {}
    tripeptides = []    

    #get list of tripeptides
    for i in amino_acids:
        for j in amino_acids:
            for k in amino_acids:
                tripeptides.append(i + j + k)

    #get fractional frequency of each tripeptide (percentage), consistent with dipeptide_composition
    total_tripeptides = len(sequence) - 2
    counts = Counter(sequence[i:i+3] for i in range(max(total_tripeptides, 0)))
    for i in tripeptides:
        tripeptide_composition[i] = round(counts[i] / total_tripeptides * 100, 3) if total_tripeptides > 0 else 0.0

    #transform values and columns to DataFrame
    tripeptide_composition_df = pd.DataFrame([list(tripeptide_composition.values())], 
        columns=list(tripeptide_composition.keys()))

    return tripeptide_composition_df

################################## GRAVY ##################################################

#Kyte-Doolittle hydropathy scale for all 20 standard amino acids [8]
kyte_doolittle_ = {
    "A":  1.8, "C":  2.5, "D": -3.5, "E": -3.5, "F":  2.8, "G": -0.4,
    "H": -3.2, "I":  4.5, "K": -3.9, "L":  3.8, "M":  1.9, "N": -3.5,
    "P": -1.6, "Q": -3.5, "R": -4.5, "S": -0.8, "T": -0.7, "V":  4.2,
    "W": -0.9, "Y": -1.3
}

def gravy(sequence: str, scale: str = "kyte_doolittle", normalized: bool = False) -> pd.DataFrame:
    """
    Calculate the Grand Average of Hydropathy (GRAVY) for a protein sequence.
    GRAVY is the sum of the Kyte-Doolittle hydropathy values of all residues
    divided by the sequence length [8]:

    GRAVY = (sum of hydropathy values) / N

    A positive GRAVY value indicates an overall hydrophobic protein; a negative
    value indicates an overall hydrophilic protein.

    Another single-value scale can be averaged with ``scale`` ("hopp_woods",
    "eisenberg", "wimley_white", "engelman", "domcut" or any AAindex1 accession
    code); the column is then named ``GRAVY_{scale}``.

    Parameters
    ==========
    :sequence: str
        protein sequence.
    :scale: str (default="kyte_doolittle")
        hydrophobicity scale to average.
    :normalized: bool (default=False)
        rescale to [0, 1] with the scale's minimum and maximum residue values
        (GRAVY is a mean, so it lies between them); columns get a ``_norm`` suffix.

    Returns
    =======
    :gravy_df: pd.DataFrame
        pandas DataFrame containing the GRAVY score for the input sequence.
        Dataframe will be of the shape 1 x 1.
    """
    #validate input protein sequence
    sequence = _validate_sequence(sequence)

    if scale != "kyte_doolittle":
        values = _property_matrix(sequence, scale)
        if values.shape[1] != 1:
            raise ValueError(f"Scale '{scale}' is multi-dimensional; use a single-value scale.")
        df = _to_df([values.mean()], [f"GRAVY_{scale}"])
        return _minmax(df, *_scale_range(scale)) if normalized else df

    #sum Kyte-Doolittle hydropathy values and divide by sequence length
    gravy_score = round(sum(kyte_doolittle_[aa] for aa in sequence) / len(sequence), 3)

    #transform to DataFrame
    gravy_df = pd.DataFrame([[gravy_score]], columns=["GRAVY"])

    if normalized:
        return _minmax(gravy_df, min(kyte_doolittle_.values()), max(kyte_doolittle_.values()))
    return gravy_df

################################## Aromaticity ############################################

def aromaticity(sequence: str) -> pd.DataFrame:
    """
    Calculate Aromaticity of a protein sequence. Aromaticity is the fraction of
    aromatic amino acids (Phe, Trp, Tyr) in the sequence [9]:

    Aromaticity = (F + W + Y) / N

    where F, W and Y are the counts of Phe, Trp and Tyr residues, and N is the
    total sequence length.

    Parameters
    ==========
    :sequence: str
        protein sequence.

    Returns
    =======
    :aromaticity_df: pd.DataFrame
        pandas DataFrame of shape 1 x 1 containing the aromaticity score.
    """
    #validate input protein sequence
    sequence = _validate_sequence(sequence)

    #count aromatic residues (F, W, Y) and divide by sequence length
    aromatic_count = sum(sequence.count(aa) for aa in ("F", "W", "Y"))
    aromaticity_score = round(aromatic_count / len(sequence), 3)

    #transform to DataFrame
    aromaticity_df = pd.DataFrame([[aromaticity_score]], columns=["Aromaticity"])

    return aromaticity_df

############################### Instability Index #########################################

#DIWV (dipeptide instability weight values) from Guruprasad et al. (1990) [10]
#missing dipeptide pairs default to 1.0 (neutral contribution)
diwv_ = {
    "WW": 1.0,  "WC": 1.0,  "WM": 24.68, "WH": 24.68, "WY": 1.0,  "WF": 1.0,  "WQ": 1.0,
    "WN": 13.34,"WK": 1.0,  "WD": 1.0,   "WT": -14.03,"WS": 1.0,  "WR": 1.0,  "WV": -7.49,
    "WP": 1.0,  "WA": -14.03,"WI": 1.0,  "WL": 13.34, "WG": -9.37,"WE": 1.0,
    "CW": 24.68,"CC": -6.54, "CM": 33.6,  "CH": 33.6,  "CY": 1.0,  "CF": 1.0,  "CQ": -6.54,
    "CN": 1.0,  "CK": 1.0,  "CD": 20.26, "CT": 33.6,  "CS": 1.0,  "CR": 1.0,  "CV": -6.54,
    "CP": 20.26,"CA": 1.0,  "CI": 1.0,   "CL": 20.26, "CG": 1.0,  "CE": 1.0,
    "MW": 1.0,  "MC": 1.0,  "MM": -1.88, "MH": 58.28, "MY": 24.68,"MF": 1.0,  "MQ": -6.54,
    "MN": 1.0,  "MK": 1.0,  "MD": 1.0,   "MT": -1.88, "MS": 44.94,"MR": -6.54,"MV": 1.0,
    "MP": 44.94,"MA": 13.34,"MI": 1.0,   "ML": 1.0,   "MG": 1.0,  "ME": 1.0,
    "HW": -1.88,"HC": 1.0,  "HM": 1.0,   "HH": -1.88, "HY": 44.94,"HF": -9.37,"HQ": 1.0,
    "HN": 24.68,"HK": 24.68,"HD": 1.0,   "HT": -6.54, "HS": -6.54,"HR": 1.0,  "HV": 1.0,
    "HP": -1.88,"HA": 1.0,  "HI": 44.94, "HL": 1.0,   "HG": -9.37,"HE": 1.0,
    "YW": -9.37,"YC": 1.0,  "YM": 44.94, "YH": 13.34, "YY": 13.34,"YF": 1.0,  "YQ": 1.0,
    "YN": 1.0,  "YK": -9.37,"YD": 24.68, "YT": -7.49, "YS": 1.0,  "YR": -15.91,"YV": 1.0,
    "YP": 13.34,"YA": 24.68,"YI": 1.0,   "YL": 1.0,   "YG": -7.49,"YE": -6.54,
    "FW": 1.0,  "FC": 1.0,  "FM": 1.0,   "FH": 1.0,   "FY": 33.6, "FF": 1.0,  "FQ": 1.0,
    "FN": 1.0,  "FK": -14.03,"FD": 13.34,"FT": 1.0,   "FS": 1.0,  "FR": 1.0,  "FV": 1.0,
    "FP": 20.26,"FA": 1.0,  "FI": 1.0,   "FL": 1.0,   "FG": 1.0,  "FE": 1.0,
    "QW": 1.0,  "QC": -6.54,"QM": 1.0,   "QH": 1.0,   "QY": -6.54,"QF": -6.54,"QQ": 20.26,
    "QN": 1.0,  "QK": 1.0,  "QD": 20.26, "QT": 1.0,   "QS": 44.94,"QR": 1.0,  "QV": -6.54,
    "QP": 20.26,"QA": 1.0,  "QI": 1.0,   "QL": 1.0,   "QG": 1.0,  "QE": 20.26,
    "NW": -9.37,"NC": -1.88,"NM": 1.0,   "NH": 1.0,   "NY": 1.0,  "NF": -14.03,"NQ": -6.54,
    "NN": 1.0,  "NK": 24.68,"ND": 1.0,   "NT": -7.49, "NS": 1.0,  "NR": 1.0,  "NV": 1.0,
    "NP": -1.88,"NA": 1.0,  "NI": 44.94, "NL": 1.0,   "NG": -14.03,"NE": 1.0,
    "KW": 1.0,  "KC": 1.0,  "KM": 33.6,  "KH": -1.88, "KY": 1.0,  "KF": 1.0,  "KQ": 24.68,
    "KN": 1.0,  "KK": 1.0,  "KD": 1.0,   "KT": 1.0,   "KS": 1.0,  "KR": 33.6, "KV": -7.49,
    "KP": -6.54,"KA": 1.0,  "KI": -7.49, "KL": -7.49, "KG": -7.49,"KE": 1.0,
    "DW": 1.0,  "DC": 1.0,  "DM": 1.0,   "DH": 1.0,   "DY": 1.0,  "DF": 1.0,  "DQ": 1.0,
    "DN": 42.67,"DK": -7.49,"DD": 1.0,   "DT": 1.0,   "DS": 1.0,  "DR": -6.54,"DV": 1.0,
    "DP": 1.0,  "DA": 1.0,  "DI": 1.0,   "DL": 1.0,   "DG": 1.0,  "DE": 1.0,
    "TW": -14.03,"TC": 1.0, "TM": 1.0,   "TH": 1.0,   "TY": 1.0,  "TF": 13.34,"TQ": -6.54,
    "TN": -14.03,"TK": 1.0, "TD": 1.0,   "TT": 1.0,   "TS": 1.0,  "TR": 1.0,  "TV": 1.0,
    "TP": 1.0,  "TA": 13.34,"TI": 1.0,   "TL": 1.0,   "TG": -7.49,"TE": 20.26,
    "SW": 1.0,  "SC": 33.6, "SM": 1.0,   "SH": 1.0,   "SY": 1.0,  "SF": 1.0,  "SQ": 20.26,
    "SN": 1.0,  "SK": 1.0,  "SD": 1.0,   "ST": 1.0,   "SS": 20.26,"SR": 20.26,"SV": 1.0,
    "SP": 44.94,"SA": 20.26,"SI": 1.0,   "SL": 1.0,   "SG": 1.0,  "SE": 20.26,
    "RW": 58.28,"RC": 1.0,  "RM": 1.0,   "RH": 20.26, "RY": -6.54,"RF": 1.0,  "RQ": 20.26,
    "RN": 13.34,"RK": 1.0,  "RD": -6.54, "RT": 1.0,   "RS": 44.94,"RR": 58.28,"RV": 1.0,
    "RP": 20.26,"RA": 1.0,  "RI": 1.0,   "RL": 1.0,   "RG": -7.49,"RE": 1.0,
    "VW": 1.0,  "VC": 1.0,  "VM": 1.0,   "VH": 1.0,   "VY": -6.54,"VF": 1.0,  "VQ": 1.0,
    "VN": 1.0,  "VK": -7.49,"VD": -14.03,"VT": -7.49, "VS": 1.0,  "VR": 1.0,  "VV": 1.0,
    "VP": 20.26,"VA": 1.0,  "VI": 1.0,   "VL": 1.0,   "VG": -7.49,"VE": 1.0,
    "PW": -6.54,"PC": -6.54,"PM": -6.54, "PH": 1.0,   "PY": 1.0,  "PF": 20.26,"PQ": 20.26,
    "PN": 1.0,  "PK": 1.0,  "PD": -6.54, "PT": 1.0,   "PS": 20.26,"PR": -6.54,"PV": 20.26,
    "PP": 20.26,"PA": 20.26,"PI": 1.0,   "PL": 1.0,   "PG": 1.0,  "PE": 18.38,
    "AW": -14.03,"AC": 44.94,"AM": 1.0,  "AH": -7.49, "AY": 1.0,  "AF": 1.0,  "AQ": 1.0,
    "AN": 1.0,  "AK": 1.0,  "AD": -7.49, "AT": 1.0,   "AS": 1.0,  "AR": 1.0,  "AV": 1.0,
    "AP": 20.26,"AA": 1.0,  "AI": 1.0,   "AL": 1.0,   "AG": 1.0,  "AE": 1.0,
    "IW": 1.0,  "IC": 1.0,  "IM": 1.0,   "IH": 13.34, "IY": 1.0,  "IF": 1.0,  "IQ": 1.0,
    "IN": 1.0,  "IK": -7.49,"ID": 1.0,   "IT": 1.0,   "IS": 1.0,  "IR": 1.0,  "IV": -7.49,
    "IP": 1.0,  "IA": 1.0,  "II": 1.0,   "IL": 1.0,   "IG": 1.0,  "IE": 44.94,
    "LW": 24.68,"LC": 1.0,  "LM": 1.0,   "LH": 1.0,   "LY": 1.0,  "LF": 1.0,  "LQ": 33.6,
    "LN": 1.0,  "LK": -7.49,"LD": 1.0,   "LT": 1.0,   "LS": 1.0,  "LR": 20.26,"LV": 1.0,
    "LP": 20.26,"LA": 1.0,  "LI": 1.0,   "LL": 1.0,   "LG": 1.0,  "LE": 1.0,
    "GW": 13.34,"GC": 1.0,  "GM": 1.0,   "GH": 1.0,   "GY": -7.49,"GF": 1.0,  "GQ": 1.0,
    "GN": -7.49,"GK": -7.49,"GD": 1.0,   "GT": -7.49, "GS": 1.0,  "GR": 1.0,  "GV": 1.0,
    "GP": 1.0,  "GA": 1.0,  "GI": 1.0,   "GL": 1.0,   "GG": 13.34,"GE": -6.54,
    "EW": -14.03,"EC": 44.94,"EM": 1.0,  "EH": -6.54, "EY": 1.0,  "EF": 1.0,  "EQ": 20.26,
    "EN": 1.0,  "EK": 1.0,  "ED": 1.0,   "ET": -6.54, "ES": 20.26,"ER": 1.0,  "EV": 1.0,
    "EP": 20.26,"EA": 1.0,  "EI": 1.0,   "EL": 1.0,   "EG": 1.0,  "EE": 20.26
}

def instability_index(sequence: str, normalized: bool = False) -> pd.DataFrame:
    """
    Calculate the Instability Index (II) of a protein sequence. The instability
    index is a measure of in vivo stability, computed as a weighted sum of
    dipeptide instability weight values (DIWV) [10]:

    II = (10 / N) * sum(DIWV(x_i, x_{i+1})) for i in 1..N-1

    Proteins with II < 40 are predicted to be stable; II >= 40 indicates
    instability.

    Parameters
    ==========
    :sequence: str
        protein sequence.
    :normalized: bool (default=False)
        rescale to [0, 1] with the bounds 10 x min(DIWV, 0) and 10 x max(DIWV, 0)
        (II = 10 x sum of L - 1 dipeptide weights / L); column ``InstabilityIndex_norm``.

    Returns
    =======
    :instability_index_df: pd.DataFrame
        pandas DataFrame of shape 1 x 1 containing the instability index.
    """
    #validate input protein sequence
    sequence = _validate_sequence(sequence)

    #sum DIWV values for all consecutive dipeptides; default to 1.0 for unlisted pairs
    dipeptide_sum = sum(diwv_.get(sequence[i:i+2], 1.0) for i in range(len(sequence) - 1))
    ii = round((dipeptide_sum * 10.0) / len(sequence), 3)

    #transform to DataFrame
    instability_index_df = pd.DataFrame([[ii]], columns=["InstabilityIndex"])

    if normalized:
        weights = list(diwv_.values()) + [1.0]
        return _minmax(instability_index_df, 10 * min(min(weights), 0), 10 * max(max(weights), 0))
    return instability_index_df

################################# Isoelectric Point ######################################

#pKa values for ionisable groups used in the isoelectric point calculation [11]
pka_values_ = {
    "D": 3.9, "E": 4.1, "H": 6.0, "C": 8.3, "Y": 10.1, "K": 10.5, "R": 12.5,
    "N_term": 8.0, "C_term": 3.1
}

def _pka_set(pka_scale: str) -> dict:
    """Return the pKa set named ``pka_scale`` from residue_tables.json."""
    sets = _load_json_table("residue_tables.json", "pka_sets")["values"]
    if pka_scale not in sets:
        raise ValueError(f"Invalid pka_scale '{pka_scale}'. Valid scales are: {list(sets)}.")
    return sets[pka_scale]

def isoelectric_point(sequence: str, pka_scale: str = "protpy", normalized: bool = False) -> pd.DataFrame:
    """
    Calculate the theoretical Isoelectric Point (pI) of a protein sequence.
    pI is the pH at which the net charge of the protein is zero. It is computed
    by bisecting the pH range 0-14 until positive and negative charges balance [11].

    ``pka_scale`` selects the pKa values: "protpy" (default, unchanged), "emboss"
    (EMBOSS Epk.dat), "bjellqvist" (Bjellqvist 1993/1994, including the
    residue-specific N-terminal pKa values for A, M, S, P, T, V, E and C-terminal
    values for D and E) or "ipc_protein" (Kozlowski 2016). The same charge model
    (K, R, H, N-terminus positive; D, E, C, Y, C-terminus negative) is used for
    every scale.

    Parameters
    ==========
    :sequence: str
        protein sequence.
    :pka_scale: str (default="protpy")
        pKa value set: "protpy", "emboss", "bjellqvist" or "ipc_protein".
    :normalized: bool (default=False)
        rescale the pI from the searched pH range 0-14 to [0, 1] (pI / 14);
        column ``IsoelectricPoint_norm``.

    Returns
    =======
    :isoelectric_point_df: pd.DataFrame
        pandas DataFrame of shape 1 x 1 containing the pI value.
    """
    #validate input protein sequence
    sequence = _validate_sequence(sequence)

    if pka_scale == "protpy":
        pka = pka_values_
    else:
        pka = dict(_pka_set(pka_scale))
        #residue-specific terminal pKa values (Bjellqvist)
        pka["N_term"] = pka.get("N_term_by_residue", {}).get(sequence[0], pka["N_term"])
        pka["C_term"] = pka.get("C_term_by_residue", {}).get(sequence[-1], pka["C_term"])

    #count ionisable residues
    counts = {aa: sequence.count(aa) for aa in "DEHCYKR"}

    def _net_charge(ph: float) -> float:
        pos = (counts["K"] / (1 + 10 ** (ph - pka["K"])) +
               counts["R"] / (1 + 10 ** (ph - pka["R"])) +
               counts["H"] / (1 + 10 ** (ph - pka["H"])) +
               1.0 / (1 + 10 ** (ph - pka["N_term"])))
        neg = (counts["D"] / (1 + 10 ** (pka["D"] - ph)) +
               counts["E"] / (1 + 10 ** (pka["E"] - ph)) +
               counts["C"] / (1 + 10 ** (pka["C"] - ph)) +
               counts["Y"] / (1 + 10 ** (pka["Y"] - ph)) +
               1.0 / (1 + 10 ** (pka["C_term"] - ph)))
        return pos - neg

    #net charge decreases monotonically with pH, so bisect over the full 0-14 range
    low, high = 0.0, 14.0
    for _ in range(50):
        mid = (low + high) / 2
        if _net_charge(mid) > 0:
            low = mid
        else:
            high = mid
    ph = (low + high) / 2

    #transform to DataFrame
    isoelectric_point_df = pd.DataFrame([[round(ph, 3)]], columns=["IsoelectricPoint"])

    if normalized:
        return _minmax(isoelectric_point_df, 0.0, 14.0)
    return isoelectric_point_df

############################### Molecular Weight ##########################################

#average residue molecular weights (Da) for each of the 20 standard amino acids [12]
residue_molecular_weight_ = {
    "A": 89.094,  "C": 121.154, "D": 133.104, "E": 147.130, "F": 165.192,
    "G": 75.067,  "H": 155.156, "I": 131.174, "K": 146.189, "L": 131.174,
    "M": 149.208, "N": 132.119, "P": 115.132, "Q": 146.146, "R": 174.203,
    "S": 105.093, "T": 119.119, "V": 117.148, "W": 204.228, "Y": 181.191
}

def molecular_weight(sequence: str, normalized: bool = False) -> pd.DataFrame:
    """
    Calculate the Molecular Weight (MW) of a protein sequence. MW is computed
    as the sum of average residue masses minus one water molecule per peptide
    bond (18.015 Da) [12]:

    MW = sum(residue_mass_i) - (N - 1) * 18.015

    Parameters
    ==========
    :sequence: str
        protein sequence.
    :normalized: bool (default=False)
        rescale to [0, 1] between the lightest (all-Gly) and heaviest (all-Trp)
        sequences of the same length, which reduces to
        (mean residue mass - min mass) / (max mass - min mass); column
        ``MolecularWeight_norm``.

    Returns
    =======
    :molecular_weight_df: pd.DataFrame
        pandas DataFrame of shape 1 x 1 containing the molecular weight (Da).
    """
    #validate input protein sequence
    sequence = _validate_sequence(sequence)

    #sum residue masses and subtract one water per peptide bond
    mw = sum(residue_molecular_weight_[aa] for aa in sequence) - (len(sequence) - 1) * 18.015
    molecular_weight_df = pd.DataFrame([[round(mw, 3)]], columns=["MolecularWeight"])

    if normalized:
        L = len(sequence)
        lightest = L * min(residue_molecular_weight_.values()) - (L - 1) * 18.015
        heaviest = L * max(residue_molecular_weight_.values()) - (L - 1) * 18.015
        return _minmax(molecular_weight_df, lightest, heaviest)
    return molecular_weight_df

############################## Charge Distribution ########################################

def charge_distribution(sequence: str, ph: float = 7.4, pka_scale: str = "protpy",
    normalized: bool = False) -> pd.DataFrame:
    """
    Calculate Charge Distribution of a protein sequence at a given pH. Positive
    charge is contributed by Lys (K), Arg (R) and His (H); negative charge by
    Asp (D) and Glu (E). Henderson-Hasselbalch equations are used [13]:

    positive = sum(count(aa) / (1 + 10^(pH - pKa))) for aa in {K, R, H}
    negative = sum(count(aa) / (1 + 10^(pKa - pH))) for aa in {D, E}
    net      = positive - negative

    Parameters
    ==========
    :sequence: str
        protein sequence.
    :ph: float (default=7.4)
        pH at which to evaluate charges.
    :pka_scale: str (default="protpy")
        pKa value set for K, R, H, D and E: "protpy", "emboss", "bjellqvist" or
        "ipc_protein" (see ``isoelectric_point``). Termini are not included.
    :normalized: bool (default=False)
        rescale to [0, 1] at the given pH: PositiveCharge / (L x f+),
        NegativeCharge / (L x f-) and (NetCharge + L x f-) / (L x (f+ + f-)), where
        f+ is the charge of the most basic residue (highest pKa) and f- that of the
        most acidic (lowest pKa); columns get a ``_norm`` suffix.

    Returns
    =======
    :charge_distribution_df: pd.DataFrame
        pandas DataFrame of shape 1 x 3 with columns
        ["PositiveCharge", "NegativeCharge", "NetCharge"].
    """
    #validate input protein sequence
    sequence = _validate_sequence(sequence)

    #pKa values for charged residues
    pka_charged_ = {"K": 10.5, "R": 12.5, "H": 6.0, "D": 3.9, "E": 4.1}
    if pka_scale != "protpy":
        pka_charged_ = {aa: _pka_set(pka_scale)[aa] for aa in pka_charged_}

    #compute positive charge (K, R, H) using Henderson-Hasselbalch
    pos = sum(sequence.count(aa) / (1 + 10 ** (ph - pka_charged_[aa])) for aa in ("K", "R", "H"))
    #compute negative charge (D, E)
    neg = sum(sequence.count(aa) / (1 + 10 ** (pka_charged_[aa] - ph)) for aa in ("D", "E"))

    charge_distribution_df = pd.DataFrame(
        [[round(pos, 3), round(neg, 3), round(pos - neg, 3)]],
        columns=["PositiveCharge", "NegativeCharge", "NetCharge"]
    )

    if normalized:
        L = len(sequence)
        f_pos = max(1 / (1 + 10 ** (ph - pka_charged_[aa])) for aa in ("K", "R", "H"))
        f_neg = max(1 / (1 + 10 ** (pka_charged_[aa] - ph)) for aa in ("D", "E"))
        return _minmax(charge_distribution_df, [0, 0, -L * f_neg], [L * f_pos, L * f_neg, L * f_pos])
    return charge_distribution_df

####################### Hydrophobic/Polar/Charged Composition ############################

#residue class membership for hydrophobic, polar and charged groupings
hydrophobic_residues_ = frozenset("ACFILMVWY")
polar_residues_ = frozenset("GNQST")
charged_residues_ = frozenset("DERHK")

def hydrophobic_polar_charged_composition(sequence: str) -> pd.DataFrame:
    """
    Calculate Hydrophobic, Polar and Charged Composition of a protein sequence.
    Residues are grouped into three physicochemical classes [1]:

    - Hydrophobic (nonpolar): A, C, F, I, L, M, V, W, Y
    - Polar (uncharged):       G, N, Q, S, T
    - Charged:                 D, E, R, H, K

    Each value is expressed as a percentage of the total sequence length.

    Parameters
    ==========
    :sequence: str
        protein sequence.

    Returns
    =======
    :hpc_df: pd.DataFrame
        pandas DataFrame of shape 1 x 3 with columns
        ["Hydrophobic", "Polar", "Charged"].
    """
    #validate input protein sequence
    sequence = _validate_sequence(sequence)

    n = len(sequence)
    #compute fraction of each class as percentage
    h = round(sum(1 for aa in sequence if aa in hydrophobic_residues_) / n * 100, 3)
    p = round(sum(1 for aa in sequence if aa in polar_residues_) / n * 100, 3)
    c = round(sum(1 for aa in sequence if aa in charged_residues_) / n * 100, 3)

    hpc_df = pd.DataFrame([[h, p, c]], columns=["Hydrophobic", "Polar", "Charged"])

    return hpc_df

########################## Secondary Structure Propensity ################################

#Chou-Fasman propensity values for alpha-helix, beta-sheet and coil [14]
helix_propensity_ = {
    "A": 1.45, "C": 0.77, "D": 0.98, "E": 1.53, "F": 1.12, "G": 0.53, "H": 1.24,
    "I": 1.00, "K": 1.07, "L": 1.34, "M": 1.20, "N": 0.73, "P": 0.59, "Q": 1.17,
    "R": 0.79, "S": 0.79, "T": 0.82, "V": 1.14, "W": 1.14, "Y": 0.61
}
sheet_propensity_ = {
    "A": 0.97, "C": 1.30, "D": 0.80, "E": 0.26, "F": 1.28, "G": 0.81, "H": 0.71,
    "I": 1.60, "K": 0.74, "L": 1.22, "M": 1.67, "N": 0.65, "P": 0.62, "Q": 1.23,
    "R": 0.90, "S": 0.72, "T": 1.20, "V": 1.65, "W": 1.19, "Y": 1.29
}
coil_propensity_ = {
    "A": 0.72, "C": 1.17, "D": 1.41, "E": 0.80, "F": 0.58, "G": 1.64, "H": 1.00,
    "I": 0.51, "K": 1.01, "L": 0.57, "M": 0.67, "N": 1.68, "P": 1.91, "Q": 0.98,
    "R": 1.24, "S": 1.33, "T": 1.03, "V": 0.50, "W": 0.99, "Y": 1.29
}

def secondary_structure_propensity(sequence: str, normalized: bool = False) -> pd.DataFrame:
    """
    Calculate Secondary Structure Propensity of a protein sequence. Each residue
    is assigned Chou-Fasman propensity values for alpha-helix (H), beta-sheet (E)
    and coil/loop (C) conformations [14]. The mean propensity for each secondary
    structure class across the whole sequence is returned.

    Parameters
    ==========
    :sequence: str
        protein sequence.
    :normalized: bool (default=False)
        rescale each column to [0, 1] with the minimum and maximum residue
        propensity of its table; columns get a ``_norm`` suffix.

    Returns
    =======
    :ssp_df: pd.DataFrame
        pandas DataFrame of shape 1 x 3 with columns ["Helix", "Sheet", "Coil"].
    """
    #validate input protein sequence
    sequence = _validate_sequence(sequence)

    n = len(sequence)
    #average propensity per secondary structure class
    h = round(sum(helix_propensity_[aa] for aa in sequence) / n, 3)
    e = round(sum(sheet_propensity_[aa] for aa in sequence) / n, 3)
    c = round(sum(coil_propensity_[aa] for aa in sequence) / n, 3)

    ssp_df = pd.DataFrame([[h, e, c]], columns=["Helix", "Sheet", "Coil"])

    if normalized:
        tables = (helix_propensity_, sheet_propensity_, coil_propensity_)
        return _minmax(ssp_df, [min(t.values()) for t in tables], [max(t.values()) for t in tables])
    return ssp_df

################################# k-mer Composition ######################################

def kmer_composition(sequence: str, k: int = 2) -> pd.DataFrame:
    """
    Calculate k-mer Composition of a protein sequence. A k-mer is a subsequence
    of length k. For each of the 20^k possible k-mers the frequency relative to
    the total number of k-mers in the sequence is calculated [1].

    Parameters
    ==========
    :sequence: str
        protein sequence.
    :k: int (default=2)
        length of k-mer subsequences (k >= 1).

    Returns
    =======
    :kmer_composition_df: pd.DataFrame
        pandas DataFrame of shape 1 x 20^k containing fractional k-mer frequencies.
    """
    #validate input protein sequence
    sequence = _validate_sequence(sequence)

    #validate k; raise error for impractically large k
    k = _validate_int("k", k, 1)
    if k > 4:
        raise ValueError(
            f"k={k} would generate {20**k:,} columns (20^k). Maximum supported k is 4."
        )

    #generate all possible k-mers from the 20 canonical amino acids
    all_kmers = ["".join(p) for p in itertools.product(amino_acids, repeat=k)]
    total = len(sequence) - k + 1

    #count every (overlapping) k-mer in one pass, divide by total k-mer count
    counts = Counter(sequence[i:i+k] for i in range(max(total, 0)))
    kmer_comp = {}
    for kmer in all_kmers:
        kmer_comp[kmer] = round(counts[kmer] / total * 100, 3) if total > 0 else 0.0

    kmer_composition_df = pd.DataFrame([list(kmer_comp.values())], columns=list(kmer_comp.keys()))

    return kmer_composition_df

########################## Reduced Alphabet Composition ##################################

#predefined reduced alphabets mapping 20 amino acids to fewer groups
#key = alphabet size, value = list of frozensets (each frozenset = one group)
reduced_alphabets_ = {
    2: [frozenset("ACFILMVWY"), frozenset("DEGNHKPQRST")],
    3: [frozenset("ACFILMVWY"), frozenset("GNQST"), frozenset("DERHKP")],
    4: [frozenset("ACFILMVW"), frozenset("GNQSTY"), frozenset("DERHK"), frozenset("P")],
    6: [
        frozenset("ACST"),   # small polar/aliphatic
        frozenset("FILMVWY"),# hydrophobic/aromatic
        frozenset("DE"),      # acidic
        frozenset("KRH"),     # basic
        frozenset("NQ"),      # amide
        frozenset("GP"),      # special
    ]
}

def reduced_alphabet_composition(sequence: str, alphabet_size: int = 6, alphabet: str | None = None,
    k: int = 1) -> pd.DataFrame:
    """
    Calculate Reduced Alphabet Composition of a protein sequence. The 20 standard
    amino acids are clustered into a smaller alphabet of physicochemically similar
    groups. The percentage of residues belonging to each group is returned [1].

    Supported alphabet sizes: 2, 3, 4, 6 (default=6).
    Invalid sizes are reset to 6.

    A named alphabet can be used instead with ``alphabet``. ``alphabet="ifeature5"``
    (aliphatic GAVLMI, aromatic FYW, positive KRH, negative DE, uncharged STCPNQ)
    gives grouped amino acid composition (GAAC, k=1, 1 x 5), grouped dipeptide
    composition (GDPC, k=2, 1 x 25) and grouped tripeptide composition (GTPC,
    k=3, 1 x 125) [21]. The cluster profiles "cp13", "cp14", "cp19" and "cp20" [24]
    are also accepted. With k > 1 overlapping k-mers of groups are counted and
    divided by (L - k + 1):

    GkPC(g1..gk) = N(g1..gk) / (L - k + 1) x 100

    Parameters
    ==========
    :sequence: str
        protein sequence.
    :alphabet_size: int (default=6)
        number of amino acid groups (2, 3, 4 or 6); ignored when ``alphabet`` is set.
    :alphabet: str | None (default=None)
        named alphabet: "ifeature5", "cp13", "cp14", "cp19" or "cp20".
    :k: int (default=1)
        group k-mer length, 1-3.

    Returns
    =======
    :reduced_alphabet_df: pd.DataFrame
        pandas DataFrame of shape 1 x G^k, where G is the number of groups. Values
        are percentages.
    """
    #validate input protein sequence
    sequence = _validate_sequence(sequence)
    k = _validate_int("k", k, 1)
    if k > 3:
        raise ValueError(f"k must be 1, 2 or 3, got {k}.")

    if alphabet is None:
        #default to 6 if an unsupported size is requested
        if alphabet_size not in reduced_alphabets_:
            alphabet_size = 6
        groups = [(str(idx + 1), group) for idx, group in enumerate(reduced_alphabets_[alphabet_size])]
        prefix = "ReducedAlphabet"
    else:
        groups = _named_alphabet(alphabet)
        prefix = {1: "GAAC", 2: "GDPC", 3: "GTPC"}[k] if alphabet == "ifeature5" else f"RA_{alphabet}"

    n = len(sequence)
    if k == 1 and alphabet is None:
        values = []
        col_names = []
        #compute percentage composition for each group
        for label, group in groups:
            values.append(round(sum(1 for aa in sequence if aa in group) / n * 100, 3))
            col_names.append(f"{prefix}_{label}")
        return pd.DataFrame([values], columns=col_names)

    mapping = {aa: label for label, group in groups for aa in group}
    counts = _kmer_counts(sequence, k, mapping)
    total = n - k + 1
    labels = [label for label, _ in groups]
    combos = list(itertools.product(labels, repeat=k))
    values = [counts[c] / total * 100 if total > 0 else 0.0 for c in combos]
    return _to_df(values, [f"{prefix}_" + ".".join(c) for c in combos])


def _named_alphabet(alphabet: str) -> list[tuple[str, str]]:
    """Return [(label, residues), ...] for a named partition alphabet in reduced_alphabets.json."""
    valid = ["ifeature5", "cp13", "cp14", "cp19", "cp20"]
    if alphabet not in valid:
        raise ValueError(f"Invalid alphabet '{alphabet}'. Valid alphabets are: {valid}.")
    values = _load_json_table("reduced_alphabets.json", alphabet)["values"]
    if isinstance(values, dict):
        return list(values.items())
    return [(group, group) for group in values]

############################## Motif Composition ##########################################

#default set of biologically relevant motifs with regex patterns [15]
default_motifs_ = {
    "NxS/T_glycosylation": r"N[^P][ST]",
    "RGD_integrin":        r"RGD",
    "KDEL_retention":      r"KDEL",
    "CxxC_zinc_finger":    r"C.{2}C",
    "CAAX_prenylation":    r"C[ACFILMVWY]{2}[ACEILMQV]$",
    "cAMP_PKA":            r"RR.S",
    "dileucine_sorting":   r"[DE]...L[LI]",
    "PEST_region":         r"P.{1,10}[ESD].{1,10}T"
}

def motif_composition(sequence: str, motifs: dict[str, str] | None = None) -> pd.DataFrame:
    """
    Calculate Motif Composition of a protein sequence. For each motif in the
    provided dictionary, the number of regex matches within the sequence is
    counted [15]. If no motifs are supplied, a set of biologically relevant
    default motifs is used.

    Parameters
    ==========
    :sequence: str
        protein sequence.
    :motifs: dict or None (default=None)
        dictionary mapping motif names (str) to regex patterns (str). If None, the built-in default_motifs dict is used.

    Returns
    =======
    :motif_composition_df: pd.DataFrame
        pandas DataFrame of shape 1 x M, where M is the number of motifs.
        Each value is the integer count of matches for that motif.
    """
    #validate input protein sequence
    sequence = _validate_sequence(sequence)

    #use default motifs if none provided
    if motifs is None:
        motifs = default_motifs_

    motif_counts = {}
    #count overlapping matches for each motif pattern
    for name, pattern in motifs.items():
        motif_counts[name] = len(re.findall(r'(?=' + pattern + r')', sequence))

    motif_composition_df = pd.DataFrame([list(motif_counts.values())], columns=list(motif_counts.keys()))

    return motif_composition_df

######################### Amino Acid Pair Composition #####################################

def amino_acid_pair_composition(sequence: str) -> pd.DataFrame:
    """
    Calculate Amino Acid Pair Composition (PairComp) of a protein sequence.
    For each pair of consecutive residues (i, i+1), the fractional frequency
    of all 400 possible dipeptide types (20 x 20) is computed [1]:

    PairComp(s,t) = count(s,t) / (N - 1) * 100

    Unlike the standard dipeptide composition, results are annotated with the
    physicochemical class of each residue in the pair (Hydrophobic/Polar/Charged).

    Parameters
    ==========
    :sequence: str
        protein sequence.

    Returns
    =======
    :pair_composition_df: pd.DataFrame
        pandas DataFrame of shape 1 x 400, where each column label has the form
        'XY_Class1-Class2' (e.g. 'AL_Hydrophobic-Hydrophobic').
    """
    #validate input protein sequence
    sequence = _validate_sequence(sequence)

    #map each amino acid to its physicochemical class
    def _aa_class(aa):
        if aa in hydrophobic_residues_:
            return "Hydrophobic"
        elif aa in polar_residues_:
            return "Polar"
        return "Charged"

    total = len(sequence) - 1
    counts = Counter(sequence[k:k+2] for k in range(total))
    pair_comp = {}
    #build the 400 annotated columns and compute frequencies
    for i in amino_acids:
        for j in amino_acids:
            label = f"{i}{j}_{_aa_class(i)}-{_aa_class(j)}"
            pair_comp[label] = round(counts[i + j] / total * 100, 3) if total > 0 else 0.0

    pair_composition_df = pd.DataFrame([list(pair_comp.values())], columns=list(pair_comp.keys()))

    return pair_composition_df

################################## Aliphatic Index ########################################

def aliphatic_index(sequence: str, normalized: bool = False) -> pd.DataFrame:
    """
    Calculate the Aliphatic Index (AI) of a protein sequence. The aliphatic
    index is a measure of thermostability, defined as the relative volume
    occupied by aliphatic side chains (Ala, Val, Ile, Leu) [16]:

    AI = Xala + 2.9 * Xval + 3.9 * (Xile + Xleu)

    where X values are mole percentages of each residue. Higher values
    indicate greater thermostability.

    Parameters
    ==========
    :sequence: str
        protein sequence.
    :normalized: bool (default=False)
        rescale to [0, 1] by the maximum of 390 (an all-Ile/Leu sequence);
        column ``AliphaticIndex_norm``.

    Returns
    =======
    :aliphatic_index_df: pd.DataFrame
        pandas DataFrame of shape 1 x 1 containing the aliphatic index score.
    """
    #validate input protein sequence
    sequence = _validate_sequence(sequence)

    n = len(sequence)
    #compute mole percentages for aliphatic residues
    xa = sequence.count("A") / n * 100
    xv = sequence.count("V") / n * 100
    xi = sequence.count("I") / n * 100
    xl = sequence.count("L") / n * 100

    #apply Ikai (1980) aliphatic index formula
    ai = round(xa + 2.9 * xv + 3.9 * (xi + xl), 3)

    aliphatic_index_df = pd.DataFrame([[ai]], columns=["AliphaticIndex"])

    if normalized:
        return _minmax(aliphatic_index_df, 0.0, 390.0)
    return aliphatic_index_df

############################### Extinction Coefficient ####################################

def extinction_coefficient(sequence: str, absorbance: bool = False, normalized: bool = False) -> pd.DataFrame:
    """
    Calculate the molar Extinction Coefficient of a protein sequence at 280 nm.
    The coefficient is derived from the number of Trp (W), Tyr (Y) and
    Cys (C) residues using the method of Gasteiger et al. (2005) [17]:

    EC = (W * 5500) + (Y * 1490) + (SS * 125)

    where SS is the number of disulfide bonds (Cys_count // 2) for the
    oxidised form, and 0 for the reduced form. Both are returned.

    With ``absorbance=True`` the A280 of a 1 g/L (0.1%) solution is added,
    Abs 0.1% = EC / MW, for both forms (Pace et al. 1995, Gill & von Hippel 1989).

    Parameters
    ==========
    :sequence: str
        protein sequence.
    :absorbance: bool (default=False)
        also return ``Abs01_Reduced`` and ``Abs01_Oxidized``.
    :normalized: bool (default=False)
        rescale to [0, 1]: coefficients by the all-Trp maximum 5500 x L, and
        absorbances by the per-residue maximum 5500 / (204.228 - 18.015) (Trp's
        coefficient over its residue mass); columns get a ``_norm`` suffix.

    Returns
    =======
    :extinction_coefficient_df: pd.DataFrame
        pandas DataFrame of shape 1 x 2 (1 x 4 with ``absorbance``) with columns
        ``ExtCoeff_Reduced`` and ``ExtCoeff_Oxidized`` (M⁻¹cm⁻¹).
    """
    #validate input protein sequence
    sequence = _validate_sequence(sequence)

    w = sequence.count("W")
    y = sequence.count("Y")
    c = sequence.count("C")

    #reduced: no disulfide bonds
    ec_reduced = w * 5500 + y * 1490
    #oxidized: each pair of Cys contributes one disulfide bond (125 M⁻¹cm⁻¹)
    ec_oxidized = ec_reduced + (c // 2) * 125

    extinction_coefficient_df = pd.DataFrame(
        [[ec_reduced, ec_oxidized]],
        columns=["ExtCoeff_Reduced", "ExtCoeff_Oxidized"]
    )

    if absorbance:
        mw = molecular_weight(sequence)["MolecularWeight"].iloc[0]
        extinction_coefficient_df["Abs01_Reduced"] = round(ec_reduced / mw, 3)
        extinction_coefficient_df["Abs01_Oxidized"] = round(ec_oxidized / mw, 3)

    if normalized:
        hi = [5500 * len(sequence)] * 2 + [5500 / (residue_molecular_weight_["W"] - 18.015)] * (2 if absorbance else 0)
        return _minmax(extinction_coefficient_df, 0.0, hi)
    return extinction_coefficient_df

################################ Boman Index ##############################################

#side-chain solubility values (Radzicka & Wolfenden) used by Boman (2003) [18]
boman_scale_ = {
    "A":  1.81, "C":  1.28, "D": -8.72, "E": -6.81, "F":  2.98, "G":  0.94,
    "H": -4.66, "I":  4.92, "K": -5.55, "L":  4.92, "M":  2.35, "N": -6.64,
    "P":  0.00, "Q": -5.54, "R": -14.92, "S": -3.40, "T": -2.57, "V":  4.04,
    "W":  2.33, "Y": -0.14
}

def boman_index(sequence: str, normalized: bool = False) -> pd.DataFrame:
    """
    Calculate the Boman Index (potential protein interaction index) for
    a protein sequence. The index is the sum of the solubility values of
    all residues divided by the sequence length, negated [18]:

    BI = -sum(solubility(aa_i)) / N

    Higher values indicate a greater tendency to bind other proteins or
    cell membranes (e.g. antimicrobial peptides). Values >2.48 suggest
    high interaction potential.

    Parameters
    ==========
    :sequence: str
        protein sequence.
    :normalized: bool (default=False)
        rescale to [0, 1] between -max and -min of the solubility scale
        (BI is a negated mean); column ``BomanIndex_norm``.

    Returns
    =======
    :boman_index_df: pd.DataFrame
        pandas DataFrame of shape 1 x 1 containing the Boman Index score.
    """
    #validate input protein sequence
    sequence = _validate_sequence(sequence)

    #negated mean of the solubility values
    bi = round(-sum(boman_scale_[aa] for aa in sequence) / len(sequence), 3)

    boman_index_df = pd.DataFrame([[bi]], columns=["BomanIndex"])

    if normalized:
        return _minmax(boman_index_df, -max(boman_scale_.values()), -min(boman_scale_.values()))
    return boman_index_df

############################### Aggregation Propensity ####################################

def aggregation_propensity(sequence: str, window: int = 5, hydrophobicity_threshold: float = 2.0, charge_threshold: int = 0) -> pd.DataFrame:
    """
    Calculate the Aggregation Propensity of a protein sequence. Aggregation-
    prone regions (APRs) are identified using a sliding window: a window is
    classified as an APR when its mean Kyte-Doolittle hydrophobicity exceeds
    ``hydrophobicity_threshold`` and the number of charged residues (D, E, K,
    R) is at most ``charge_threshold`` [8].

    Parameters
    ==========
    :sequence: str
        protein sequence.
    :window: int (default=5)
        sliding window length in residues.
    :hydrophobicity_threshold: float (default=2.0)
        minimum mean Kyte-Doolittle hydrophobicity for an APR window.
    :charge_threshold: int (default=0)
        maximum number of charged residues (D, E, K, R) allowed in an APR.

    Returns
    =======
    :aggregation_propensity_df: pd.DataFrame
        pandas DataFrame of shape 1 x 2 with columns:
        ``AggregProneRegions`` — number of non-overlapping-scored APR windows;
        ``AggregProneFraction`` — percentage of residues covered by any APR.
    """
    #validate input protein sequence
    sequence = _validate_sequence(sequence)

    window = _validate_int("window", window, 1)

    n = len(sequence)
    charged = frozenset("DEKR")
    #track which residue positions are covered by at least one APR
    covered = [False] * n
    apr_count = 0

    #score each window
    for i in range(n - window + 1):
        seg = sequence[i:i + window]
        avg_hydro = sum(kyte_doolittle_[aa] for aa in seg) / window
        n_charged = sum(1 for aa in seg if aa in charged)
        if avg_hydro >= hydrophobicity_threshold and n_charged <= charge_threshold:
            apr_count += 1
            for j in range(i, i + window):
                covered[j] = True

    apr_fraction = round(sum(covered) / n * 100, 3)

    aggregation_propensity_df = pd.DataFrame(
        [[apr_count, apr_fraction]],
        columns=["AggregProneRegions", "AggregProneFraction"]
    )

    return aggregation_propensity_df

################################## Shannon Entropy ########################################

def shannon_entropy(sequence: str, normalized: bool = False) -> pd.DataFrame:
    """
    Calculate the Shannon Entropy of a protein sequence. Shannon entropy
    quantifies the diversity of the amino acid composition using the
    information-theoretic formula:

        H = -sum_i ( p_i * log2(p_i) )

    where p_i is the fractional frequency of amino acid i in the sequence.
    A value near zero indicates a low-complexity or repetitive sequence; the
    theoretical maximum is log2(20) ≈ 4.322 bits for a perfectly uniform
    distribution across all 20 canonical amino acids. Shannon entropy is
    widely used as a sequence quality filter and complexity measure in
    machine-learning pipelines [1].

    Parameters
    ==========
    :sequence: str
        protein sequence.
    :normalized: bool (default=False)
        rescale to [0, 1] by the maximum log2(20) bits; column
        ``ShannonEntropy_norm``.

    Returns
    =======
    :shannon_entropy_df: pd.DataFrame
        pandas DataFrame of shape 1 x 1 with column ``ShannonEntropy``.
        Value is rounded to 3 decimal places.
    """
    #validate input protein sequence
    sequence = _validate_sequence(sequence)

    n = len(sequence)
    #count occurrences of each amino acid present in the sequence
    counts = {}
    for aa in sequence:
        counts[aa] = counts.get(aa, 0) + 1

    #compute entropy using Shannon formula: H = -sum(p_i * log2(p_i))
    entropy = 0.0
    for count in counts.values():
        p = count / n
        entropy -= p * math.log2(p)

    shannon_entropy_df = pd.DataFrame([[round(entropy, 3)]], columns=["ShannonEntropy"])

    if normalized:
        return _minmax(shannon_entropy_df, 0.0, math.log2(20))
    return shannon_entropy_df

############################### Hydrophobic Moment ########################################

def hydrophobic_moment(sequence: str, window: int = 11, angle: float | str = 100,
    scale: str | None = None, normalized: bool = False) -> pd.DataFrame:
    """
    Calculate the mean and maximum Hydrophobic Moment of a protein sequence.
    The hydrophobic moment measures the amphiphilicity of a helical segment —
    the directional asymmetry of hydrophobicity projected around a helical
    axis [19].

    For each window of length ``window``, the moment is:

        mu = sqrt( (sum H_i * sin(i * angle))^2 + (sum H_i * cos(i * angle))^2 ) / window

    where H_i is the hydrophobicity of residue i and angle is in degrees.
    Default angle of 100° corresponds to an alpha-helix; 160° is used for
    beta-strands (Eisenberg et al. 1984).

    By default H is protpy's normalised hydrophobicity table (A = 0.62), which is
    the normalised consensus scale of Eisenberg, Schwarz, Komarony & Wall (1984,
    J Mol Biol 179:125-142) as listed by ExPASy ProtScale. ``scale`` selects
    another single-value scale, e.g. "eisenberg" (AAindex EISD840101, Eisenberg
    1984 Annu Rev Biochem, A = 0.25), "kyte_doolittle" or any AAindex1 code.

    Parameters
    ==========
    :sequence: str
        protein sequence.
    :window: int (default=11)
        sliding window length in residues.
    :angle: float | str (default=100)
        residue rotation angle in degrees, or "helix" (100) / "sheet" (160).
    :scale: str | None (default=None)
        hydrophobicity scale; None uses protpy's default table.
    :normalized: bool (default=False)
        rescale to [0, 1] by the largest absolute value of the hydrophobicity
        scale (a window moment cannot exceed it); columns get a ``_norm`` suffix.

    Returns
    =======
    :hydrophobic_moment_df: pd.DataFrame
        pandas DataFrame of shape 1 x 2 with columns
        ``HydrophobicMoment_Mean`` and ``HydrophobicMoment_Max``.
    """
    #validate input protein sequence
    sequence = _validate_sequence(sequence)

    window = _validate_int("window", window, 1)
    if isinstance(angle, str):
        if angle not in ("helix", "sheet"):
            raise ValueError(f"angle must be a number, 'helix' or 'sheet', got '{angle}'.")
        angle = {"helix": 100, "sheet": 160}[angle]
    if scale is None:
        h_scale = hydrophobicity_
    else:
        values = _property_matrix("".join(amino_acids), scale)
        if values.shape[1] != 1:
            raise ValueError(f"Scale '{scale}' is multi-dimensional; use a single-value scale.")
        h_scale = dict(zip(amino_acids, values[:, 0]))

    n = len(sequence)
    #use full sequence as single window when shorter than default window
    effective_window = min(window, n)
    angle_rad = angle * math.pi / 180.0
    moments = []

    #compute hydrophobic moment for each window position using Eisenberg scale
    for i in range(n - effective_window + 1):
        seg = sequence[i:i + effective_window]
        hsin = sum(h_scale[aa] * math.sin(j * angle_rad) for j, aa in enumerate(seg))
        hcos = sum(h_scale[aa] * math.cos(j * angle_rad) for j, aa in enumerate(seg))
        mu = math.sqrt(hsin ** 2 + hcos ** 2) / effective_window
        moments.append(mu)

    mean_moment = round(sum(moments) / len(moments), 3)
    max_moment = round(max(moments), 3)

    hydrophobic_moment_df = pd.DataFrame(
        [[mean_moment, max_moment]],
        columns=["HydrophobicMoment_Mean", "HydrophobicMoment_Max"]
    )

    if normalized:
        return _minmax(hydrophobic_moment_df, 0.0, max(abs(v) for v in h_scale.values()))
    return hydrophobic_moment_df

############################## Pseudo Amino Acid Composition ##############################

def pseudo_amino_acid_composition(sequence: str, lamda: int = 30, weight: float = 0.05, properties: list[str] | None = None,
    normalized: bool = False) -> pd.DataFrame:
    """
    Pseudo amino acid composition (PAAComp) combines the vanilla amino acid composition 
    descriptor with additional local features, such as correlation between residues of 
    a certain distance, as amino acid composition doesn't take into account sequence 
    order info. The pseudo components of the descriptor are a series rank-different 
    correlation factors [4, 5]. 
    
    The first 20 components are a weighted sum of the amino acid composition and 30 are 
    physicochemical square correlations as dictated by the lamda and properties parameters. 
    This generates an output of [1, (20 + lamda)] = 1 x 50 when using the default lamda 
    of 30. By default, the physicochemical properties used are hydrophobicity and 
    hydrophillicity, with a lamda of 30 and weight of 0.05.

    Parameters
    ==========
    :lamda: int
        rank of correlation. Number of calculable descriptors depends on lamda value.
    :weight: float
        weighting factor.
    :properties: list 
        list of dicts of physicochemical/structural property values for amino acids.
    :normalized: bool (default=False)
        rescale to [0, 1] by dividing by 100 (components are non-negative
        percentages summing to 100); columns get a ``_norm`` suffix.

    Returns
    =======
    :pseudo_amino_acid_composition_df: pd.Dataframe
        output dataframe of calculated pseudo amino acid composition descriptors 
        for input sequence. Dataframe will be of the shape 1 x N, where N is the 
        number of features calculated from the descriptor (20 + lamda) and 1 is 
        the input sequence. By default, the shape will be 1 x 50 (using default lamda=30).
    """
    #validate input protein sequence
    sequence = _validate_sequence(sequence)

    #validate lamda & weight; lamda falls back to its default of 30 when > sequence length
    lamda = _validate_int("lamda", lamda)
    if lamda > len(sequence):
        lamda = 30
    weight = _validate_weight(weight)

    #normalize None to empty list before further processing
    if properties is None:
        properties = []

    #cast properties to list 
    if not (isinstance(properties, list)):
        properties = [properties]
    
    property_values = []

    #by default use properties (hydrphocity, hydrophilicty, residue mass),
    #otherwise get the individual property values from aaindex1 using accession numbers
    if not properties:
        property_values = [hydrophobicity_, hydrophilicity_, residue_mass_]
    else:
        for prop in properties:
            if (prop not in aaindex1.record_codes()):
                raise ValueError(f"Accession number not found in aaindex1: {prop}.")
            #keep only the 20 canonical amino acids (drops the '-' gap entry)
            property_values.append(
                {aa: v for aa, v in aaindex1[prop].values.items() if aa in amino_acids})

    psuedo_aac = {}

    #calculate amino acid composition
    aa_comp = amino_acid_composition(sequence)

    #compute correlation factors once; reuse for both the first-20 and last-lamda features
    rightpart = [sequence_order_correlation_factor(sequence, i + 1, property_values)
                 for i in range(lamda)]

    #shared normalization denominator: 1 + w * Σθ_k  (Chou 2001)
    temp = 1 + weight * sum(rightpart)

    #compute first 20 pseudo amino acid descriptor components (weighted AAComp)
    for index, i in enumerate(amino_acids):
        psuedo_aac["PAAC_" + str(index + 1)] = round(aa_comp[i].values[0] / temp, 3)

    #compute last lamda pseudo amino acid descriptor components (sequence-order correlations)
    for index in range(20, 20 + lamda):
        psuedo_aac["PAAC_" + str(index + 1)] = round(
            weight * rightpart[index - 20] / temp * 100, 3)

    #transform values and columns to DataFrame
    psuedo_amino_acid_composition_df = pd.DataFrame([list(psuedo_aac.values())], 
        columns=list(psuedo_aac.keys()))

    if normalized:
        return _minmax(psuedo_amino_acid_composition_df, 0.0, 100.0)
    return psuedo_amino_acid_composition_df 

def sequence_order_correlation_factor(sequence: str, k: int = 1, properties: list[dict[str, float]] | None = None) -> float:
    """
    Calculating sequence order correlation factor with gap equal to k based on
    the given input properities for a protein sequence.

    Parameters
    ==========
    :k: int
        gap between amino acids in the sequence.
    :properties: list
        list of dicts of physicochemical/structural property values for amino acids.

    Returns
    =======
    :result: float
        correlation factor value for the sequence using the correlation function
        for the particular sequence.
    """
    # default to empty list if not provided (avoid mutable default argument)
    if properties is None:
        properties = []

    # If there are no valid residue pairs for the requested gap, return 0.0.
    if len(sequence) <= k:
        return 0.0

    #normalize each property once rather than for every residue pair
    normalized_props = [normalize_property(prop) for prop in properties]
    correlation_factor = []

    #iterate through sequence calculating the correlation (as in correlation_function)
    #between its amino acids with gap=k
    for i in range(len(sequence) - k):
        aa1 = sequence[i]
        aa2 = sequence[i + k]
        theta = sum(math.pow(norm[aa1] - norm[aa2], 2) for norm in normalized_props)
        correlation_factor.append(round(theta / len(normalized_props), 3))

    #round array of correlation factors 
    correlation_factor_ = round(sum(correlation_factor) / (len(sequence) - k), 3)

    return correlation_factor_    

def correlation_function(Ri: str = "S", Rj: str = "D", properties: list[dict[str, float]] | None = None) -> float:
    """
    Calculate the correlation between 2 input amino acids based on the 
    physicochemical/structural properties from the protein sequence.

    Parameters
    ==========
    :Ri: str
        1st amino acid.
    :Rj: str
        2nd amino acid.
    :properties: list 
        list of dicts of physicochemical/structural property values for amino acids.

    Returns
    =======
    :correlation: float
        correlation value for the two input amino acids based on input properties.
    """
    # default to empty list if not provided (avoid mutable default argument)
    if properties is None:
        properties = []

    theta = 0.0

    #iterate over each property, calculating the correlation from the normalzied values
    for i in range(len(properties)):
        temp = normalize_property(properties[i])
        theta = theta + math.pow(temp[Ri] - temp[Rj], 2)    

    #calculate correlation for all properties
    correlation = round(theta / len(properties), 3)
    
    return correlation

def normalize_property(properties: dict[str, float]) -> dict[str, float]:
    """
    Normalize physicochemical/structural property values using
    their mean and standard deviation.

    Parameters
    ==========
    :properties: dict
        dictionary of amino acid values for a physicochemical property.

    Returns
    =======
    :normalized_vals: dict
        dict of normalized property values.
    """
    vals = list(properties.values())
    if not vals:
        return {}

    mean = sum(vals) / len(vals)
    std = _std(vals, mean, ddof=0)
    if std == 0:
        return {aa: 0.0 for aa, _ in properties.items()}
    return {aa: (v - mean) / std for aa, v in properties.items()}

def _std(prop: list[float], mean: float, ddof: int = 1) -> float:
    """
    Calculate standard deviation from list of physicochemical values using
    the mean of their values.

    Parameters
    ==========
    :prop: dict
        dictionary of property values for each amino acid.
    :mean: float
        mean of physicochemical property values.
    :ddof: int (default=1)
        delta degrees of freedom.
        
    Returns
    =======
    :std_: float
        calculated standard deviation. 
    """
    temp = [math.pow(i - mean, 2) for i in prop]
    std_ = math.sqrt(sum(temp) / (len(prop) - ddof))
    return std_

######################## Amphiphilic Pseudo Amino Acid Composition ########################

def amphiphilic_pseudo_amino_acid_composition(sequence: str, lamda: int = 30, weight: float = 0.5,
    properties: list[dict[str, float]] | None = None) -> pd.DataFrame:
    """
    Amphiphillic pseudo amino acid composition (APAAComp) has the same form as the amino
    acid composition, but contains much more information that is related to the 
    sequence order of a protein and the distribution of the hydrophobic and 
    hydrophilic amino acids along its chain [5].
    
    The first 20 numbers in the descriptor are the components of the conventional amino 
    acid composition; the next 2*λ numbers are a set of correlation factors that reflect
    different hydrophobicity and hydrophilicity distribution patterns along a protein chain.

    Parameters
    ==========
    :sequence: str
        protein sequence.
    :lamda: int
        rank of correlation. Number of calculable descriptors depends on lambda value.
    :weight: float
        weighting factor.
    :properties: list[dict] | None (default=None → [``hydrophobicity_``, ``hydrophilicity_``])
        list of property dicts mapping each amino acid to a float value.
        When None the default hydrophobicity and hydrophilicity scales are used.

    Returns
    =======
    :amp_pseudo_amino_acid_composition_df: pd.Dataframe
        output dataframe of calculated amphiphilic pseudo amino acid composition
        descriptors for input sequence. Dataframe will be of the shape 1 x N, where
        N is 20 + n_props * lamda (n_props = number of properties supplied).
        With the default two properties and lamda=30 the shape is 1 x 80.
    """
    if properties is None:
        properties = [hydrophobicity_, hydrophilicity_]

    #validate input protein sequence
    sequence = _validate_sequence(sequence)

    #validate lamda & weight; lamda falls back to its default of 30 when > sequence length
    lamda = _validate_int("lamda", lamda)
    if lamda > len(sequence):
        lamda = 30
    weight = _validate_weight(weight)

    #cast properties to list
    if not isinstance(properties, list):
        properties = [properties]

    n_props = len(properties)
    amp_pseudo_amino_acid_composition = {}

    #compute all per-lag correlation factors in a single pass:
    #rightpart = [p0_lag1, p1_lag1, p0_lag2, p1_lag2, ..., p{n-1}_lag{lamda}]
    rightpart: list[float] = []
    for i in range(lamda):
        rightpart.extend(
            amphiphilic_sequence_order_correlation_factor(sequence, k=i + 1, properties=properties))

    #calculate amino acid composition
    aa_composition = amino_acid_composition(sequence)

    #shared normalization denominator: 1 + w * Σ(all correlation factors)
    temp = 1 + weight * sum(rightpart)

    #compute first 20 pseudo amino acid descriptor components (weighted AAComp)
    for index, i in enumerate(amino_acids):
        amp_pseudo_amino_acid_composition["APAAC_" + str(index + 1)] = round(
            aa_composition[i].values[0] / temp, 3)

    #compute remaining n_props * lamda amphiphilic composition values
    for index in range(20, 20 + n_props * lamda):
        amp_pseudo_amino_acid_composition["APAAC_" + str(index + 1)] = round(
            weight * rightpart[index - 20] / temp * 100, 3)

    #transform values and columns to DataFrame
    amp_pseudo_amino_acid_composition_df = pd.DataFrame([list(amp_pseudo_amino_acid_composition.values())], 
        columns=list(amp_pseudo_amino_acid_composition.keys()))

    return amp_pseudo_amino_acid_composition_df

def amphiphilic_sequence_order_correlation_factor(sequence: str, k: int = 1,
    properties: list[dict[str, float]] | None = None) -> list[float]:
    """
    Calculate Amphiphilic sequence order correlation factor for sequence
    with gap=k.

    Parameters
    ==========
    :sequence: str
        protein sequence.
    :k: int (default=1)
        gap between amino acids in the sequence.
    :properties: list[dict] | None (default=None → [``hydrophobicity_``, ``hydrophilicity_``])
        list of property dicts. Pre-normalized once per call.

    Returns
    =======
    :correlation_factor: list[float]
        one correlation factor per property.
    """
    if properties is None:
        properties = [hydrophobicity_, hydrophilicity_]

    #no residue pairs for the requested gap: every correlation factor is 0
    denom = len(sequence) - k
    if denom <= 0:
        return [0.0] * len(properties)

    #pre-normalize all properties once rather than inside every pair evaluation
    normalized_props = [normalize_property(prop) for prop in properties]
    n_props = len(properties)
    accumulators: list[list[float]] = [[] for _ in range(n_props)]

    #iterate over sequence pairs with the given gap k
    for i in range(len(sequence) - k):
        aa1, aa2 = sequence[i], sequence[i + k]
        for p_idx, norm in enumerate(normalized_props):
            accumulators[p_idx].append(round(norm[aa1] * norm[aa2], 3))

    return [round(sum(acc) / denom, 3) for acc in accumulators]

def amphiphilic_correllation_function(Ri: str = "S", Rj: str = "D",
    normalized_props: list[dict[str, float]] | None = None) -> list[float]:
    """
    Calculate correlation value based on physicochemical property values
    for input amino acids Ri and Rj.

    Parameters
    ==========
    :Ri: str
        1st amino acid.
    :Rj: str
        2nd amino acid.
    :normalized_props: list[dict] | None
        pre-normalized property dicts. Defaults to normalized hydrophobicity
        and hydrophilicity when None.

    Returns
    =======
    :thetas: list[float]
        one correlation value per property.
    """
    if normalized_props is None:
        normalized_props = [normalize_property(hydrophobicity_),
                            normalize_property(hydrophilicity_)]
    return [round(prop[Ri] * prop[Rj], 3) for prop in normalized_props]
####################### Composition of k-Spaced Amino Acid Pairs ######################

def cksaap(sequence: str, gap: int = 5) -> pd.DataFrame:
    """
    Calculate the Composition of k-Spaced Amino Acid Pairs (CKSAAP) of a protein
    sequence [20]. For each spacing g = 0..gap, every residue pair (s_i, s_{i+g+1})
    is counted and divided by the number of such pairs:

    CKSAAP(XY, g) = N_g(XY) / (L - g - 1) x 100

    g = 0 is the dipeptide composition. Values are percentages; a spacing with no
    pairs (L - g - 1 <= 0) gives zeros.

    Parameters
    ==========
    :sequence: str
        protein sequence.
    :gap: int (default=5)
        maximum number of residues between the two residues of a pair.

    Returns
    =======
    :cksaap_df: pd.DataFrame
        pandas DataFrame of shape 1 x 400·(gap+1), columns ``CKSAAP_{XY}_gap{g}``.
    """
    #validate input protein sequence
    sequence = _validate_sequence(sequence)
    gap = _validate_int("gap", gap)

    pairs = list(itertools.product(amino_acids, repeat=2))
    values, columns = [], []
    for g in range(gap + 1):
        counts = _gapped_pair_counts(sequence, g)
        total = len(sequence) - g - 1
        values += [counts[p] / total * 100 if total > 0 else 0.0 for p in pairs]
        columns += [f"CKSAAP_{x}{y}_gap{g}" for x, y in pairs]
    return _to_df(values, columns)

################# Composition of k-Spaced Amino Acid Group Pairs #######################

def cksaagp(sequence: str, gap: int = 5) -> pd.DataFrame:
    """
    Calculate the Composition of k-Spaced Amino Acid Group Pairs (CKSAAGP) [21].
    As CKSAAP, but residues are first mapped to the five iFeature groups:
    aliphatic (GAVLMI), aromatic (FYW), positive (KRH), negative (DE) and
    uncharged (STCPNQ).

    CKSAAGP(G1G2, g) = N_g(G1G2) / (L - g - 1) x 100

    Values are percentages and the 25 values for each spacing sum to 100
    (0 when the spacing has no pairs).

    Parameters
    ==========
    :sequence: str
        protein sequence.
    :gap: int (default=5)
        maximum number of residues between the two residues of a pair.

    Returns
    =======
    :cksaagp_df: pd.DataFrame
        pandas DataFrame of shape 1 x 25·(gap+1), columns ``CKSAAGP_{G1}.{G2}_gap{g}``.
    """
    #validate input protein sequence
    sequence = _validate_sequence(sequence)
    gap = _validate_int("gap", gap)

    groups = _named_alphabet("ifeature5")
    mapping = {aa: label for label, residues in groups for aa in residues}
    pairs = list(itertools.product([label for label, _ in groups], repeat=2))
    values, columns = [], []
    for g in range(gap + 1):
        counts = _gapped_pair_counts(sequence, g, mapping)
        total = len(sequence) - g - 1
        values += [counts[p] / total * 100 if total > 0 else 0.0 for p in pairs]
        columns += [f"CKSAAGP_{a}.{b}_gap{g}" for a, b in pairs]
    return _to_df(values, columns)

############################ Dipeptide Deviation from Expected Mean ######################

def dde(sequence: str) -> pd.DataFrame:
    """
    Calculate the Dipeptide Deviation from Expected Mean (DDE) [22]. The observed
    dipeptide composition is standardised against the composition expected from
    the number of codons coding each residue:

    Dc(XY) = N(XY) / (L - 1)
    Tm(XY) = (C_X / 61) x (C_Y / 61)
    Tv(XY) = Tm(XY) x (1 - Tm(XY)) / (L - 1)
    DDE(XY) = (Dc(XY) - Tm(XY)) / sqrt(Tv(XY))

    where C_X is the number of sense codons for residue X in the standard genetic
    code (61 in total). DDE is a standardised statistic, not a percentage.

    Parameters
    ==========
    :sequence: str
        protein sequence, at least 2 residues.

    Returns
    =======
    :dde_df: pd.DataFrame
        pandas DataFrame of shape 1 x 400, columns ``DDE_{XY}``.
    """
    #validate input protein sequence
    sequence = _validate_sequence(sequence)
    if len(sequence) < 2:
        raise ValueError("DDE needs a sequence of at least 2 residues.")

    codons = _load_json_table("residue_tables.json", "codon_counts")["values"]
    counts = _gapped_pair_counts(sequence, 0)
    total = len(sequence) - 1
    values, columns = [], []
    for x, y in itertools.product(amino_acids, repeat=2):
        tm = (codons[x] / 61) * (codons[y] / 61)
        tv = tm * (1 - tm) / total
        values.append((counts[(x, y)] / total - tm) / math.sqrt(tv))
        columns.append(f"DDE_{x}{y}")
    return _to_df(values, columns)

################################ Adaptive Skip Dipeptide Composition #####################

def asdc(sequence: str) -> pd.DataFrame:
    """
    Calculate the Adaptive Skip Dipeptide Composition (ASDC) [23]: the composition
    of every ordered residue pair (s_i, s_j) with i < j, over all gaps:

    ASDC(XY) = sum_g N_g(XY) / sum_g sum_XY N_g(XY) x 100

    The denominator is the total number of ordered pairs, L(L-1)/2. Values are
    percentages and sum to 100 (zeros for a single residue).

    Parameters
    ==========
    :sequence: str
        protein sequence.

    Returns
    =======
    :asdc_df: pd.DataFrame
        pandas DataFrame of shape 1 x 400, columns ``ASDC_{XY}``.
    """
    #validate input protein sequence
    sequence = _validate_sequence(sequence)

    index = {aa: i for i, aa in enumerate(amino_acids)}
    #integer matmul: exact counts, and avoids spurious BLAS float warnings
    onehot = np.zeros((len(sequence), 20), dtype=np.int64)
    onehot[np.arange(len(sequence)), [index[aa] for aa in sequence]] = 1
    #rows of `before` count each residue type strictly before position j
    before = np.cumsum(onehot, axis=0) - onehot
    pair_counts = before.T @ onehot
    total = len(sequence) * (len(sequence) - 1) / 2
    values = (pair_counts / total * 100).ravel() if total else np.zeros(400)
    return _to_df(list(values), [f"ASDC_{x}{y}" for x, y in itertools.product(amino_acids, repeat=2)])

#################################### Distance Pair #######################################

def distance_pair(sequence: str, rng: int = 3, alphabet: str = "cp14") -> pd.DataFrame:
    """
    Calculate the PseAAC of distance-pairs over a reduced amino acid alphabet
    (iDNA-Prot|dis) [24, 25]. Residues are mapped to the clusters of a cluster
    profile cp(n) and the vector is

    [f^0(1..G), f^1(1..G^2), ..., f^d(1..G^2)]

    where f^0 is the cluster composition (count / L) and f^d, d = 1..rng, is the
    composition of cluster pairs (s_i, s_{i+d}) (count / (L - d)); d = 1 is the
    nearest-neighbour pair. This follows Eqs. 5-12 of [24], giving G + G^2·rng
    values (e.g. 14 + 14x14x3 = 602 for cp14, rng=3). Values are percentages.

    Alphabets: cp13, cp14, cp19, cp20 (no reduction) from Eq. 10 of [24]. The
    printed cp(14) lists E twice; the stray E in the first cluster is dropped,
    the only reading with 14 clusters.

    Parameters
    ==========
    :sequence: str
        protein sequence.
    :rng: int (default=3)
        maximum pair distance d.
    :alphabet: str (default="cp14")
        cluster profile: "cp13", "cp14", "cp19" or "cp20".

    Returns
    =======
    :distance_pair_df: pd.DataFrame
        pandas DataFrame of shape 1 x (G + G^2·rng), columns ``DP_{c}`` and
        ``DP_{c1}.{c2}_d{d}`` where c are cluster labels (their residues).
    """
    #validate input protein sequence
    sequence = _validate_sequence(sequence)
    rng = _validate_int("rng", rng)
    if alphabet not in ("cp13", "cp14", "cp19", "cp20"):
        raise ValueError(f"Invalid alphabet '{alphabet}'. Valid alphabets are: ['cp13', 'cp14', 'cp19', 'cp20'].")

    labels = [label for label, _ in _named_alphabet(alphabet)]
    mapping = {aa: label for label in labels for aa in label}
    comp = Counter(mapping[aa] for aa in sequence)
    values = [comp[c] / len(sequence) * 100 for c in labels]
    columns = [f"DP_{c}" for c in labels]
    pairs = list(itertools.product(labels, repeat=2))
    for d in range(1, rng + 1):
        counts = _gapped_pair_counts(sequence, d - 1, mapping)
        total = len(sequence) - d
        values += [counts[p] / total * 100 if total > 0 else 0.0 for p in pairs]
        columns += [f"DP_{a}.{b}_d{d}" for a, b in pairs]
    return _to_df(values, columns)

#################################### Natural Vector ######################################

def natural_vector(sequence: str) -> pd.DataFrame:
    """
    Calculate the Natural Vector of a protein sequence [26]. For each amino acid k:

    n_k  = number of occurrences of k
    mu_k = mean (1-based) position of k
    D2_k = sum_i (pos_i - mu_k)^2 / (n_k x L)

    mu_k and D2_k are 0 when k does not occur. Values are not percentages.

    Parameters
    ==========
    :sequence: str
        protein sequence.

    Returns
    =======
    :natural_vector_df: pd.DataFrame
        pandas DataFrame of shape 1 x 60, columns ``NV_n_{X}``, ``NV_mu_{X}`` and
        ``NV_D2_{X}``.
    """
    #validate input protein sequence
    sequence = _validate_sequence(sequence)

    L = len(sequence)
    positions = {aa: [] for aa in amino_acids}
    for i, aa in enumerate(sequence, start=1):
        positions[aa].append(i)
    n = [len(positions[aa]) for aa in amino_acids]
    mu = [sum(positions[aa]) / len(positions[aa]) if positions[aa] else 0.0 for aa in amino_acids]
    d2 = [sum((p - m) ** 2 for p in positions[aa]) / (len(positions[aa]) * L) if positions[aa] else 0.0
          for aa, m in zip(amino_acids, mu)]
    columns = ([f"NV_n_{aa}" for aa in amino_acids] + [f"NV_mu_{aa}" for aa in amino_acids]
               + [f"NV_D2_{aa}" for aa in amino_acids])
    return _to_df(n + mu + d2, columns)

#################################### Atomic Composition ##################################

def atomic_composition(sequence: str, subtract_water: bool = False) -> pd.DataFrame:
    """
    Calculate the Atomic Composition (ATC) of a protein sequence [27]: the
    percentage of C, H, N, O and S atoms, summed over the molecular formulas of
    the residues' free amino acids.

    ATC(e) = N_e / (N_C + N_H + N_N + N_O + N_S) x 100

    By default free amino-acid formulas are summed without removing water,
    matching Pfeature [27]. With ``subtract_water=True``, (L - 1) H2O are removed
    for the peptide bonds, giving the formula of the linear peptide (e.g. "GG"
    -> C4H8N2O3).

    Parameters
    ==========
    :sequence: str
        protein sequence.
    :subtract_water: bool (default=False)
        remove one H2O per peptide bond.

    Returns
    =======
    :atomic_composition_df: pd.DataFrame
        pandas DataFrame of shape 1 x 5, columns ``ATC_C``, ``ATC_H``, ``ATC_N``,
        ``ATC_O`` and ``ATC_S``. Values are percentages.
    """
    #validate input protein sequence
    sequence = _validate_sequence(sequence)

    formulas = _load_json_table("residue_tables.json", "atomic_formulas")["values"]
    atoms = {e: sum(formulas[aa][e] for aa in sequence) for e in "CHNOS"}
    if subtract_water:
        atoms["H"] -= 2 * (len(sequence) - 1)
        atoms["O"] -= len(sequence) - 1
    total = sum(atoms.values())
    return _to_df([atoms[e] / total * 100 for e in "CHNOS"], [f"ATC_{e}" for e in "CHNOS"])
