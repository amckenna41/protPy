################################################################################
#############                   Composition                       ##############
################################################################################

from __future__ import annotations

import re
import itertools
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
"""
#list of amino acids
from ._constants import amino_acids, _validate_sequence

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

def amino_acid_composition(sequence: str) -> pd.DataFrame:
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

    amino_acid_composition = {}

    #iterate through each amino acid, calculating each composition
    for aa in amino_acids:
        amino_acid_composition[aa] = round(float(sequence.count(aa)) / len(sequence) * 100, 3)

    #transform values and columns to DataFrame
    amino_acid_composition_df = pd.DataFrame([list(amino_acid_composition.values())], 
        columns=list(amino_acid_composition.keys()))

    return amino_acid_composition_df

################################## Dipeptide Composition ##################################

def dipeptide_composition(sequence: str) -> pd.DataFrame:
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

    dipeptide_composition = {}

    #iterate through each amino acid, calculating each dipeptide combination's composition
    for i in amino_acids:
        for j in amino_acids:
            dipep = i + j
            dipeptide_composition[dipep] = round(
                float(sequence.count(dipep)) / (len(sequence)-1) * 100, 2)

    #transform values and columns to DataFrame
    dipeptide_composition_df = pd.DataFrame([list(dipeptide_composition.values())], 
        columns=list(dipeptide_composition.keys()))

    return dipeptide_composition_df

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
        input sequence.
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

    #get frequency of each tripeptide in the sequence
    for i in tripeptides:
        tripeptide_composition[i] = len(re.findall(i, sequence))

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

def gravy(sequence: str) -> pd.DataFrame:
    """
    Calculate the Grand Average of Hydropathy (GRAVY) for a protein sequence.
    GRAVY is the sum of the Kyte-Doolittle hydropathy values of all residues
    divided by the sequence length [8]:

    GRAVY = (sum of hydropathy values) / N

    A positive GRAVY value indicates an overall hydrophobic protein; a negative
    value indicates an overall hydrophilic protein.

    Parameters
    ==========
    :sequence: str
        protein sequence.

    Returns
    =======
    :gravy_df: pd.DataFrame
        pandas DataFrame containing the GRAVY score for the input sequence.
        Dataframe will be of the shape 1 x 1.
    """
    #validate input protein sequence
    sequence = _validate_sequence(sequence)

    #sum Kyte-Doolittle hydropathy values and divide by sequence length
    gravy_score = round(sum(kyte_doolittle_[aa] for aa in sequence) / len(sequence), 3)

    #transform to DataFrame
    gravy_df = pd.DataFrame([[gravy_score]], columns=["GRAVY"])

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

def instability_index(sequence: str) -> pd.DataFrame:
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

    return instability_index_df

################################# Isoelectric Point ######################################

#pKa values for ionisable groups used in the isoelectric point calculation [11]
pka_values_ = {
    "D": 3.9, "E": 4.1, "H": 6.0, "C": 8.3, "Y": 10.1, "K": 10.5, "R": 12.5,
    "N_term": 8.0, "C_term": 3.1
}

def isoelectric_point(sequence: str) -> pd.DataFrame:
    """
    Calculate the theoretical Isoelectric Point (pI) of a protein sequence.
    pI is the pH at which the net charge of the protein is zero. It is computed
    by iteratively adjusting pH until positive and negative charges balance [11].

    Parameters
    ==========
    :sequence: str
        protein sequence.

    Returns
    =======
    :isoelectric_point_df: pd.DataFrame
        pandas DataFrame of shape 1 x 1 containing the pI value.
    """
    #validate input protein sequence
    sequence = _validate_sequence(sequence)

    #count ionisable residues
    counts = {aa: sequence.count(aa) for aa in "DEHCYKR"}

    #iterative charge-balance: adjust pH by +/- 0.001 until net charge < 1e-4
    ph = 7.0
    for _ in range(5000):
        pos = (counts["K"] / (1 + 10 ** (ph - pka_values_["K"])) +
               counts["R"] / (1 + 10 ** (ph - pka_values_["R"])) +
               counts["H"] / (1 + 10 ** (ph - pka_values_["H"])) +
               1.0 / (1 + 10 ** (ph - pka_values_["N_term"])))
        neg = (counts["D"] / (1 + 10 ** (pka_values_["D"] - ph)) +
               counts["E"] / (1 + 10 ** (pka_values_["E"] - ph)) +
               counts["C"] / (1 + 10 ** (pka_values_["C"] - ph)) +
               counts["Y"] / (1 + 10 ** (pka_values_["Y"] - ph)) +
               1.0 / (1 + 10 ** (pka_values_["C_term"] - ph)))
        net_charge = pos - neg
        if abs(net_charge) < 1e-4:
            break
        ph += 0.001 if net_charge > 0 else -0.001

    #transform to DataFrame
    isoelectric_point_df = pd.DataFrame([[round(ph, 3)]], columns=["IsoelectricPoint"])

    return isoelectric_point_df

############################### Molecular Weight ##########################################

#average residue molecular weights (Da) for each of the 20 standard amino acids [12]
residue_molecular_weight_ = {
    "A": 89.094,  "C": 121.154, "D": 133.104, "E": 147.130, "F": 165.192,
    "G": 75.032,  "H": 155.156, "I": 131.174, "K": 146.189, "L": 131.174,
    "M": 149.208, "N": 132.119, "P": 115.132, "Q": 146.146, "R": 174.203,
    "S": 105.093, "T": 119.119, "V": 117.148, "W": 204.228, "Y": 181.191
}

def molecular_weight(sequence: str) -> pd.DataFrame:
    """
    Calculate the Molecular Weight (MW) of a protein sequence. MW is computed
    as the sum of average residue masses minus one water molecule per peptide
    bond (18.015 Da) [12]:

    MW = sum(residue_mass_i) - (N - 1) * 18.015

    Parameters
    ==========
    :sequence: str
        protein sequence.

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

    return molecular_weight_df

############################## Charge Distribution ########################################

def charge_distribution(sequence: str, ph: float = 7.4) -> pd.DataFrame:
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

    #compute positive charge (K, R, H) using Henderson-Hasselbalch
    pos = sum(sequence.count(aa) / (1 + 10 ** (ph - pka_charged_[aa])) for aa in ("K", "R", "H"))
    #compute negative charge (D, E)
    neg = sum(sequence.count(aa) / (1 + 10 ** (pka_charged_[aa] - ph)) for aa in ("D", "E"))

    charge_distribution_df = pd.DataFrame(
        [[round(pos, 3), round(neg, 3), round(pos - neg, 3)]],
        columns=["PositiveCharge", "NegativeCharge", "NetCharge"]
    )

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

def secondary_structure_propensity(sequence: str) -> pd.DataFrame:
    """
    Calculate Secondary Structure Propensity of a protein sequence. Each residue
    is assigned Chou-Fasman propensity values for alpha-helix (H), beta-sheet (E)
    and coil/loop (C) conformations [14]. The mean propensity for each secondary
    structure class across the whole sequence is returned.

    Parameters
    ==========
    :sequence: str
        protein sequence.

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

    #validate k; default to 2 if invalid; raise error for impractically large k
    if not isinstance(k, int) or k < 1:
        k = 2
    if k > 4:
        raise ValueError(
            f"k={k} would generate {20**k:,} columns (20^k). Maximum supported k is 4."
        )

    #generate all possible k-mers from the 20 canonical amino acids
    all_kmers = ["".join(p) for p in itertools.product(amino_acids, repeat=k)]
    total = len(sequence) - k + 1

    kmer_comp = {}
    #count each k-mer occurrence, divide by total k-mer count
    for kmer in all_kmers:
        cnt = sum(1 for i in range(total) if sequence[i:i+k] == kmer)
        kmer_comp[kmer] = round(cnt / total * 100, 3) if total > 0 else 0.0

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

def reduced_alphabet_composition(sequence: str, alphabet_size: int = 6) -> pd.DataFrame:
    """
    Calculate Reduced Alphabet Composition of a protein sequence. The 20 standard
    amino acids are clustered into a smaller alphabet of physicochemically similar
    groups. The fraction of residues belonging to each group is returned [1].

    Supported alphabet sizes: 2, 3, 4, 6 (default=6).
    Invalid sizes are reset to 6.

    Parameters
    ==========
    :sequence: str
        protein sequence.
    :alphabet_size: int (default=6)
        number of amino acid groups (2, 3, 4 or 6).

    Returns
    =======
    :reduced_alphabet_df: pd.DataFrame
        pandas DataFrame of shape 1 x alphabet_size.
    """
    #validate input protein sequence
    sequence = _validate_sequence(sequence)

    #default to 6 if an unsupported size is requested
    if alphabet_size not in reduced_alphabets_:
        alphabet_size = 6

    groups = reduced_alphabets_[alphabet_size]
    n = len(sequence)
    values = []
    col_names = []
    #compute fractional composition for each group
    for idx, group in enumerate(groups):
        freq = round(sum(1 for aa in sequence if aa in group) / n * 100, 3)
        values.append(freq)
        col_names.append(f"ReducedAlphabet_{idx + 1}")

    reduced_alphabet_df = pd.DataFrame([values], columns=col_names)

    return reduced_alphabet_df

############################## Motif Composition ##########################################

#default set of biologically relevant motifs with regex patterns [15]
default_motifs_ = {
    "NxS/T_glycosylation": r"N[^P][ST]",
    "RGD_integrin":        r"RGD",
    "KDEL_retention":      r"KDEL",
    "CxxC_zinc_finger":    r"C.{2}C",
    "CAAX_prenylation":    r"C[ACFILMVWY]{2}[ACEILMQV]$",
    "cAMP_PKA":            r"RR.S",
    "dileucine_sorting":   r"[DE]xxxL[LI]$",
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
    pair_comp = {}
    #build the 400 annotated columns and compute frequencies
    for i in amino_acids:
        for j in amino_acids:
            label = f"{i}{j}_{_aa_class(i)}-{_aa_class(j)}"
            count = sum(1 for k in range(total) if sequence[k] == i and sequence[k+1] == j)
            pair_comp[label] = round(count / total * 100, 3) if total > 0 else 0.0

    pair_composition_df = pd.DataFrame([list(pair_comp.values())], columns=list(pair_comp.keys()))

    return pair_composition_df

################################## Aliphatic Index ########################################

def aliphatic_index(sequence: str) -> pd.DataFrame:
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

    return aliphatic_index_df

############################### Extinction Coefficient ####################################

def extinction_coefficient(sequence: str) -> pd.DataFrame:
    """
    Calculate the molar Extinction Coefficient of a protein sequence at 280 nm.
    The coefficient is derived from the number of Trp (W), Tyr (Y) and
    Cys (C) residues using the method of Gasteiger et al. (2005) [17]:

    EC = (W * 5500) + (Y * 1490) + (SS * 125)

    where SS is the number of disulfide bonds (Cys_count // 2) for the
    oxidised form, and 0 for the reduced form. Both are returned.

    Parameters
    ==========
    :sequence: str
        protein sequence.

    Returns
    =======
    :extinction_coefficient_df: pd.DataFrame
        pandas DataFrame of shape 1 x 2 with columns
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

    return extinction_coefficient_df

################################ Boman Index ##############################################

#solubility scale from Boman (2003) [18]
boman_scale_ = {
    "A":  2.1, "C": -0.7, "D": -3.5, "E": -3.5, "F":  2.8, "G":  2.1,
    "H": -3.2, "I":  4.5, "K": -3.9, "L":  4.5, "M":  1.9, "N": -4.8,
    "P": -1.6, "Q": -3.5, "R": -0.8, "S": -0.8, "T": -0.7, "V":  4.2,
    "W": -0.9, "Y": -1.3
}

def boman_index(sequence: str) -> pd.DataFrame:
    """
    Calculate the Boman Index (potential protein interaction index) for
    a protein sequence. The index is the sum of the solubility values of
    all residues divided by the sequence length [18]:

    BI = sum(solubility(aa_i)) / N

    Higher values indicate a greater tendency to bind other proteins or
    cell membranes (e.g. antimicrobial peptides). Values >2.48 suggest
    high interaction potential.

    Parameters
    ==========
    :sequence: str
        protein sequence.

    Returns
    =======
    :boman_index_df: pd.DataFrame
        pandas DataFrame of shape 1 x 1 containing the Boman Index score.
    """
    #validate input protein sequence
    sequence = _validate_sequence(sequence)

    #sum solubility values and divide by sequence length
    bi = round(sum(boman_scale_[aa] for aa in sequence) / len(sequence), 3)

    boman_index_df = pd.DataFrame([[bi]], columns=["BomanIndex"])

    return boman_index_df

############################### Aggregation Propensity ####################################

def aggregation_propensity(sequence: str, window: int = 5, hydrophobicity_threshold: float = 2.0, charge_threshold: int = 1) -> pd.DataFrame:
    """
    Calculate the Aggregation Propensity of a protein sequence. Aggregation-
    prone regions (APRs) are identified using a sliding window: a window is
    classified as an APR when its mean Kyte-Doolittle hydrophobicity exceeds
    ``hydrophobicity_threshold`` and the number of charged residues (D, E, K,
    R) is below ``charge_threshold`` [8].

    Parameters
    ==========
    :sequence: str
        protein sequence.
    :window: int (default=5)
        sliding window length in residues.
    :hydrophobicity_threshold: float (default=2.0)
        minimum mean Kyte-Doolittle hydrophobicity for an APR window.
    :charge_threshold: int (default=1)
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
        if avg_hydro >= hydrophobicity_threshold and n_charged < charge_threshold:
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

def shannon_entropy(sequence: str) -> pd.DataFrame:
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

    return shannon_entropy_df

############################### Hydrophobic Moment ########################################

def hydrophobic_moment(sequence: str, window: int = 11, angle: float = 100) -> pd.DataFrame:
    """
    Calculate the mean and maximum Hydrophobic Moment of a protein sequence.
    The hydrophobic moment measures the amphiphilicity of a helical segment —
    the directional asymmetry of hydrophobicity projected around a helical
    axis [19].

    For each window of length ``window``, the moment is:

        mu = sqrt( (sum H_i * sin(i * angle))^2 + (sum H_i * cos(i * angle))^2 ) / window

    where H_i is the Eisenberg hydrophobicity of residue i and angle is in
    degrees. Default angle of 100° corresponds to an alpha-helix; 160° is
    commonly used for beta-strands.

    Parameters
    ==========
    :sequence: str
        protein sequence.
    :window: int (default=11)
        sliding window length in residues.
    :angle: float (default=100)
        residue rotation angle in degrees (100° = alpha-helix).

    Returns
    =======
    :hydrophobic_moment_df: pd.DataFrame
        pandas DataFrame of shape 1 x 2 with columns
        ``HydrophobicMoment_Mean`` and ``HydrophobicMoment_Max``.
    """
    #validate input protein sequence
    sequence = _validate_sequence(sequence)

    n = len(sequence)
    #use full sequence as single window when shorter than default window
    effective_window = min(window, n)
    angle_rad = angle * math.pi / 180.0
    moments = []

    #compute hydrophobic moment for each window position using Eisenberg scale
    for i in range(n - effective_window + 1):
        seg = sequence[i:i + effective_window]
        hsin = sum(hydrophobicity_[aa] * math.sin(j * angle_rad) for j, aa in enumerate(seg))
        hcos = sum(hydrophobicity_[aa] * math.cos(j * angle_rad) for j, aa in enumerate(seg))
        mu = math.sqrt(hsin ** 2 + hcos ** 2) / effective_window
        moments.append(mu)

    mean_moment = round(sum(moments) / len(moments), 3)
    max_moment = round(max(moments), 3)

    hydrophobic_moment_df = pd.DataFrame(
        [[mean_moment, max_moment]],
        columns=["HydrophobicMoment_Mean", "HydrophobicMoment_Max"]
    )

    return hydrophobic_moment_df

############################## Pseudo Amino Acid Composition ##############################

def pseudo_amino_acid_composition(sequence: str, lamda: int = 30, weight: float = 0.05, properties: list[str] | None = None) -> pd.DataFrame:
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

    #set lamda to its default value of 30 if <0, or > sequence len or not an int
    if ((lamda < 0) or (lamda > len(sequence)) or not isinstance(lamda, int)):
        lamda = 30  

    #validate weight, set default weight of 0.05 if invalid value input
    if ((weight < 0) or not (isinstance(weight, (float, int)))):
        weight = 0.05

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
            property_values.append(aaindex1[prop].values)

    right_part = 0.0
    psuedo_aac = {}
    rightpart = []

    #calculate correlation factor for protein sequence using propeties
    for i in range(lamda):
        right_part = right_part + sequence_order_correlation_factor(sequence, i+1, property_values)

    #calculate amino acid composition
    aa_comp = amino_acid_composition(sequence)
    
    #compute first 20 pseudo amino acid descriptor components based on property_values,
    #append descriptor values to dict
    temp = 1 + weight + right_part
    for index, i in enumerate(amino_acids):
        psuedo_aac["PAAC_" + str(index + 1)] = round(aa_comp[i].values[0] / temp, 3)

    #calculate correlation factor for protein sequence using propeties
    for i in range(lamda):
        rightpart.append(sequence_order_correlation_factor(sequence, i + 1, property_values))

    #compute last 20+lambda pseudo amino acid descriptor components based on property_values,
    #append descriptor values to dict
    temp = 1 + weight * sum(rightpart)
    for index in range(20, 20 + lamda):
        psuedo_aac["PAAC_" + str(index + 1)] = round(
            weight * rightpart[index - 20] / temp * 100, 3)

    #transform values and columns to DataFrame
    psuedo_amino_acid_composition_df = pd.DataFrame([list(psuedo_aac.values())], 
        columns=list(psuedo_aac.keys()))

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
    correlation_factor = []

    #iterate through sequence using correlation function to calculate sequence order
    #between its amino acids with gap=k 
    for i in range(len(sequence) - k):
        aa1 = sequence[i]
        aa2 = sequence[i + k]
        correlation_factor.append(correlation_function(aa1, aa2, properties))
    
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
    normalized_vals = {}

    #iterate over all property values, calculate their mean and standard dev
    for i, j in properties.items():
        mean = sum(properties.values()) / len(properties.values())
        normalized_vals[i] = (j - mean) / _std(properties.values(), mean, ddof=0)

    return normalized_vals

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
    properties: list[dict[str, float]] = [hydrophobicity_, hydrophilicity_]) -> pd.DataFrame:
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
    :properties: list (default=[hydrophobicity, hydrophilicity])
        list of dicts of physicochemical/structural property values for amino acids.

    Returns
    =======    
    :amp_pseudo_amino_acid_composition_df: pd.Dataframe
        output dataframe of calculated amphiphilic pseudo amino acid composition 
        descriptors for input sequence. Dataframe will be of the shape 1 x N, where 
        N is the number of features calculated from the descriptor (20 + 2*lambda) 
        and 1 is the input sequence. By default, the shape will be 1 x 80 
        (using default lamda=30).
    """
    #validate input protein sequence
    sequence = _validate_sequence(sequence)

    #set lamda to its default value of 30 if <0, or > sequence len or not an int
    if ((lamda < 0) or (lamda > len(sequence)) or not isinstance(lamda, int)):
        lamda = 30  

    #validate weight, set default weight of 0.5 if invalid value input
    if ((weight < 0) or not (isinstance(weight, (float, int)))):
        weight = 0.5

    #cast properties to list 
    if not (isinstance(properties, list)):
        properties = [properties]

    amp_pseudo_amino_acid_composition = {}

    #calculate correlation factor in protein sequence
    rightpart = 0.0
    for i in range(lamda):
        rightpart = rightpart + sum(
            amphiphilic_sequence_order_correlation_factor(sequence, k=i + 1))
    
    #calculate amino acid composition
    aa_composition = amino_acid_composition(sequence)

    #compute first 20 pseudo amino acid descriptor components based on properties
    temp = 1 + weight * rightpart
    for index, i in enumerate(amino_acids):
        amp_pseudo_amino_acid_composition["APAAC_" + str(index + 1)] = round(aa_composition[i].values[0] / temp, 3)

    #calculate correlation factor in protein sequence
    rightpart = []
    for i in range(lamda):
        temp = amphiphilic_sequence_order_correlation_factor(sequence, k=i + 1)
        rightpart.append(temp[0])
        rightpart.append(temp[1])

    #compute remaining (20 + 2 * lamda) amphiphillic composition values
    temp = 1 + weight * sum(rightpart)
    for index in range(20, 20 + 2 * lamda):
        amp_pseudo_amino_acid_composition["APAAC_" + str(index + 1)] = round(
            weight * rightpart[index - 20] / temp * 100, 3)

    #transform values and columns to DataFrame
    amp_pseudo_amino_acid_composition_df = pd.DataFrame([list(amp_pseudo_amino_acid_composition.values())], 
        columns=list(amp_pseudo_amino_acid_composition.keys()))

    return amp_pseudo_amino_acid_composition_df

def amphiphilic_sequence_order_correlation_factor(sequence: str, k: int = 1) -> list[float]:
    """ 
    Calculate Amphipillic sequence order correlation factor for sequence
    with gap=k.

    Parameters
    ==========
    :sequence: str
        protein sequence.
    :k: int (default=1)
        gap between amino acids in the sequence.

    Returns
    =======
    :correlation_factor: list
        list of correlation factors for both hydrophobicity and hydrophillicty.
    """
    correlation_factor_hydrophobicity = []
    correlation_factor_hydrophilicity = []

    #iterate over sequence, calculating correlation factors using both properties
    for i in range(len(sequence) - k):
        aa1 = sequence[i]
        aa2 = sequence[i + k]
        corelation_function = amphiphilic_correllation_function(aa1, aa2)
        correlation_factor_hydrophobicity.append(corelation_function[0])
        correlation_factor_hydrophilicity.append(corelation_function[1])

    #append hydrophobicity and hydrophilicity values to array
    correlation_factor = []
    correlation_factor.append(round(sum(correlation_factor_hydrophobicity) / (len(sequence) - k), 3))
    correlation_factor.append(round(sum(correlation_factor_hydrophilicity) / (len(sequence) - k), 3))

    return correlation_factor

def amphiphilic_correllation_function(Ri: str = "S", Rj: str = "D") -> tuple[float, float]:
    """
    Calculate correlation value based on hydrophobicity and hydrophillicity
    property values for input amino acids Ri and Rj in sequence.

    Parameters
    ==========
    :Ri: str
        1st amino acid.
    :Rj: str
        2nd amino acid.

    Returns
    =======
    :theta1, theta2: float
        correlation values for input property values.
    """
    #normalize properties
    _hydrophobicity = normalize_property(hydrophobicity_)
    _hydrophilicity = normalize_property(hydrophilicity_)

    #calculate correlation using hydrophobicity & hydrophilicty properties
    theta1 = round(_hydrophobicity[Ri] * _hydrophobicity[Rj], 3)
    theta2 = round(_hydrophilicity[Ri] * _hydrophilicity[Rj], 3)

    return theta1, theta2