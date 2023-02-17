################################################################################
###############              Quasi Sequence Order                ###############
################################################################################

import os
import json
import pandas as pd
import math
from json import JSONDecodeError
import sys
from . import composition

"""
References
----------

[1] Kuo-Chen Chou. Prediction of Protein Subcellar Locations by Incorporating
    Quasi-Sequence-Order Effect. Biochemical and Biophysical Research Communications,
    2000, 278, 477-483.
[2] Kuo-Chen Chou and Yu-Dong Cai. Prediction of Protein Sucellular Locations by
    GO-FunD-PseAA Predictor. Biochemical and Biophysical Research Communications,
    2004, 320, 1236-1239.
[3] Gisbert Schneider and Paul Wrede. The Rational Design of Amino Acid Sequences
    by Artifical Neural Networks and Simulated Molecular Evolution: Do Novo Design
    of an Idealized Leader Cleavge Site. Biophys Journal, 1994, 66, 335-344.
[4] Grantham, R. (1974-09-06). "Amino acid difference formula to help explain protein
    evolution". Science. 185 (4154): 862–864. Bibcode:1974Sci...185..862G.
    doi:10.1126/science.185.4154.862. ISSN 0036-8075. PMID 4843792. S2CID 35388307.
"""

#list of amino acids
amino_acids = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", 
    "Q", "R", "S", "T", "V", "W", "Y"]

def sequence_order_coupling_number_(sequence, d=1, distance_matrix="schneider-wrede-physiochemical-distance-matrix"):
    """
    Calculate Sequence Order Coupling Number (SOCN) features for input protein sequence.
    SOCN computes the dissimilarity between amino acid pairs. The distance between 
    amino acid pairs is determined by d which varies between 1 to lag. For each d, it 
    computes the sum of the dissimilarities of all amino acid pairs. The output will be 
    a single float value representing the SOCN.
    
    This function should not be confused with sequence_order_coupling_number() below 
    which calculates the multiple SOCN descriptor values for the input sequence 
    according to the value of lag; with the default of lag=30, 30 SOCN's will be 
    calculated.

    Parameters
    ----------
    :sequence : str
        protein sequence.
    :d : int (default=1)
        gap between two amino acids.
    :distance_matrix : str (default="schneider-wrede-physicochemical-distance-matrix")
        path to physiochemical distance matrix for calculating quasi sequence order.

    Returns
    -------
    :socn : float
        calculated sequence order coupling number value from sequence.
    """
    #check input sequence is a string, if not raise type error
    if not isinstance(sequence, str):
        raise TypeError('Input sequence must be a string, got input of type {}.'.format(type(sequence)))

    #uppercase protein sequence 
    sequence = sequence.upper()

    #if invalid amino acids in sequence, raise value error
    for aa in sequence:
        if (aa not in amino_acids):
            raise ValueError("Invalid amino acid in protein sequence: .".format(aa))

    #append extension if missing from input
    if (os.path.splitext(distance_matrix)[1] == ''):
        distance_matrix = distance_matrix + '.json'
    
    #get filepath to distance matrix json
    if not (os.path.isfile(distance_matrix)):
        if not (os.path.isfile(os.path.join(os.path.dirname(os.path.abspath(sys.modules['protpy'].__file__)), "data", distance_matrix))):
            raise OSError('Distance Matrix json ({}) not found.'.format(
                os.path.join(os.path.dirname(os.path.abspath(sys.modules['protpy'].__file__)), "data", distance_matrix)))
        else:
            distance_matrix = os.path.join(os.path.dirname(os.path.abspath(sys.modules['protpy'].__file__)), "data", distance_matrix)
    else:
        distance_matrix = distance_matrix
    
    #open distance matrix json if present
    try:
        with open(distance_matrix, "r") as f:
            distance_matrix = json.load(f)
    except:
        raise JSONDecodeError('Error getting config JSON file: {}.'.format(distance_matrix))

    tau = 0.0

    #iterate over length of sequence minus gap, incrementing the SOCN at each iteration using
    # the values from the distance matrix with the current and next amino acid, round to 2 dp
    for aa in range(len(sequence) - d):
        current_aa = sequence[aa]
        next_aa = sequence[aa + d]
        tau = tau + math.pow(distance_matrix[current_aa + next_aa], 2)

    return round(tau, 3)

def sequence_order_coupling_number(sequence, lag=30,
    distance_matrix="schneider-wrede-physiochemical-distance-matrix.json"):
    """
    Calculate Sequence Order Coupling Number (SOCN) features for input protein sequence.
    SOCN computes the dissimilarity between amino acid pairs. The distance between 
    amino acid pairs is determined by d which varies between 1 to lag. For each d, it 
    computes the sum of the dissimilarities of all amino acid pairs. The number of 
    output features can be calculated as N * 2, where N = lag, by default this value 
    is 30 which generates an output of 1 x 60.

    Parameters
    ----------
    :sequence : str
        protein sequence.
    :lag : int (default=30)
        maximum gap betwwen 2 amino acids; the length of the protein should be larger
        than lag. Default set to 30.
    :distance_matrix : str (default="schneider-wrede-physicochemical-distance-matrix")
        path to physiochemical distance matrix for calculating quasi sequence order.

    Returns
    -------
    :sequence_order_df : pd.Dataframe
        Dataframe of SOCN descriptor values for all protein sequences. Output
        will be of the shape N x 1, where N is the number of features calculated 
        from the descriptor (calculated as M * 2 where M = lag).
    """
    #check input sequence is a string, if not raise type error
    if not isinstance(sequence, str):
        raise TypeError('Input sequence must be a string, got input of type {}.'.format(type(sequence)))

    #uppercase protein sequence 
    sequence = sequence.upper()

    #if invalid amino acids in sequence, raise value error
    for aa in sequence:
        if (aa not in amino_acids):
            raise ValueError("Invalid amino acid in protein sequence: .".format(aa))

    #raise value error if int cant be parsed from input lag
    try:
        lag = int(lag)
    except:
        raise ValueError("Invalid lag value input, integer cannot be parsed: {}.".format(lag))

    #validate lag, set default lag if invalid value input
    if (not (isinstance(lag, int) or int(lag)>=len(sequence) or (int(lag)<0))):
        lag = 30

    #set dataframe column name prefixes based on input distance matrix
    col_prefix = ""
    if (os.path.splitext(distance_matrix)[0] == "schneider-wrede-physiochemical-distance-matrix"):
        col_prefix = "SOCN_SW_"
    elif (os.path.splitext(distance_matrix)[0] == "grantham-physiochemical-distance-matrix"):
        col_prefix = "SOCN_Grant_"

    #set default lag if invalid value input
    if (lag>=len(sequence) or (lag<0) or not (isinstance(lag, int))):
        lag = 30

    seq_order = {}

    #iterate over sequence with lag, calculating SOCN using input distance matrix values
    for i in range(lag):
        tau = sequence_order_coupling_number_(sequence, i+1, distance_matrix)

        #append SOCN column and value to dict
        seq_order[col_prefix + str(i + 1)] = round(tau, 3)

    #transform SOCN data into pandas dataframe
    sequence_order_df = pd.DataFrame([list(seq_order.values())], columns=list(seq_order.keys()))

    return sequence_order_df

def sequence_order_coupling_number_all(sequence, lag=30):
    """
    Calculate Sequence Order Coupling Number (SOCN) descriptor values of input protein 
    sequence using both matrices (schneider-wrede & grantham).

    Parameters
    ----------
    :sequence : str
        protein sequence.
    :lag : int (default=30)
        maximum gap betwwen 2 amino acids; the length of the protein should be larger
        than lag. Default set to 30.

    Returns
    -------
    :socn_all : pd.Dataframe
        Concatenated dataframe of SOCN descriptor values of input protein sequence using
        both distance matrices.
    """
    #check input sequence is a string, if not raise type error
    if not isinstance(sequence, str):
        raise TypeError('Input sequence must be a string, got input of type {}.'.format(type(sequence)))

    #uppercase protein sequence 
    sequence = sequence.upper()

    #raise value error if int cant be parsed from input lag
    try:
        lag = int(lag)
    except:
        raise ValueError("Invalid lag value input, integer cannot be parsed: {}.".format(lag))

    #validate lag, set default lag if invalid value input
    if (lag>=len(sequence) or (lag<0) or not (isinstance(lag, int))):
        lag = 30

    #calculate SOCN using schneider-wrede distance matrix
    socn_schneider = sequence_order_coupling_number(sequence, lag, "schneider-wrede-physiochemical-distance-matrix.json")
    
    #calculate SOCN using grantham distance matrix
    socn_grantham = sequence_order_coupling_number(sequence, lag, "grantham-physiochemical-distance-matrix.json")

    #merge 2 dataframes 
    socn_all = pd.concat([socn_schneider, socn_grantham], axis=1)

    return socn_all

def quasi_sequence_order(sequence, lag=30, weight=0.1,
    distance_matrix="schneider-wrede-physiochemical-distance-matrix.json"):
    """
    Calculate Quasi Sequence Order (QSO) features for the protein sequences.
    The quasi-sequence-order descriptors were proposed by K.C. Chou, et.al. [1].
    They are derived from the distance matrix between the 20 amino acids. By default,
    the Scheider-Wrede physicochemical distance matrix was used. Also utilised in
    the descriptor calculation is the Grantham chemical distance matrix. Both of
    these matrices are used by Grantham et. al. in the calculation
    of the descriptor [4]. 100 values are calculated per sequence, thus generating
    an output of 1 x 100 per sequence.

    Parameters
    ----------
    :sequence : str
        protein sequence.
    :lag : int (default=30)
        maximum gap betwwen 2 amino acids; the length of the protein should be larger
        than lag. Default set to 30.
    :weight: float (default = 0.1)
        weighting factor
    :distance_matrix : str (default="schneider-wrede-physicochemical-distance-matrix")
        path to physiochemical distance matrix for calculating quasi sequence order.

    Returns
    -------
    :quasi_sequence_order_df : pd.Dataframe
        dataframe of quasi-sequence-order descriptor values for the
        protein sequences, with output shape 1 x 100 where 100 is the 
        number of calculated features per sequence.
    """
    #check input sequence is a string, if not raise type error
    if not isinstance(sequence, str):
        raise TypeError('Input sequence must be a string, got input of type {}.'.format(type(sequence)))

    #uppercase protein sequence 
    sequence = sequence.upper()

    #if invalid amino acids in sequence, raise value error
    for aa in sequence:
        if (aa not in amino_acids):
            raise ValueError("Invalid amino acid in protein sequence: .".format(aa))

    #raise value error if int cant be parsed from input lag
    try:
        lag = int(lag)
    except:
        raise ValueError("Invalid lag value input, integer cannot be parsed: {}.".format(lag))

    #validate lag, set default lag if invalid value input
    if (lag>=len(sequence) or (lag<0) or not (isinstance(lag, int))):
        lag = 30

    #validate weight, set default weight if invalid value input
    if ((weight<0) or not (isinstance(lag, float)) or not (isinstance(lag, int))):
        weight = 0.1

    #set dataframe column name prefixes based on input distance matrix
    col_prefix = ""
    if (os.path.splitext(distance_matrix)[0] == "schneider-wrede-physiochemical-distance-matrix"):
        col_prefix = "QSO_SW"
    elif (os.path.splitext(distance_matrix)[0] == "grantham-physiochemical-distance-matrix"):
        col_prefix = "QSO_Grant"
    
    quasi_sequence_order = {}
    right_part = 0.0

    #calculate first 20 quasi sequence order using sequence order coupling number for
    #proteins using lag and specificed physiochemical distance matrix
    for aa in range(lag):
        right_part = right_part + sequence_order_coupling_number_(sequence, aa+1, distance_matrix)

    #get amino acid composition
    aa_comp = composition.amino_acid_composition(sequence)

    #iterate over amino acids, calculating descriptor value 
    temp = 1 + weight * right_part
    for index, aa in enumerate(amino_acids):
        # quasi_sequence_order[col_prefix + "_1_" + str(index + 1)] = round(aa_comp[aa].values[0] / temp, 6)
        quasi_sequence_order[col_prefix + str(index + 1)] = round(aa_comp[aa].values[0] / temp, 6)
    
    #calculate SOCN up until maxlag
    _right_part = []
    for i in range(lag):
        _right_part.append(
            sequence_order_coupling_number_(sequence, i+1, distance_matrix))

    #calculate the last maxlag descriptor values
    temp = 1 + weight * sum(_right_part)
    for index in range(20, 20 + lag):
        # quasi_sequence_order[col_prefix + "_2" + str(index + 1)] = round(
        #     weight * right_part[index - 20] / temp, 6)
        quasi_sequence_order[col_prefix + str(index + 1)] = round(
            weight * _right_part[index - 20] / temp, 6)
    
    #transform descriptor data into pandas dataframe
    quasi_sequence_order_df = pd.DataFrame([list(quasi_sequence_order.values())], 
        columns=list(quasi_sequence_order.keys()))
    
    return quasi_sequence_order_df

def quasi_sequence_order_all(sequence, lag=30, weight=0.1):
    """
    Calculate Quasi Sequence Order features for input protein sequence using 
    both physiochemical distance matrices. Concatenate into one output dataframe. 
 
    Parameters
    ----------
    :sequence : str
        protein sequence.
    :lag : int (default=30)
        maximum gap between 2 amino acids; the length of the protein should be larger
        than lag. Default set to 30.
    :weight : float (default=0.1)
        weighting factor.

    Returns
    -------
    :quasi_sequence_order_all_df : pd.Dataframe
        dataframe of quasi-sequence-order descriptor values for the
        protein sequences, with output shape 1 x 100 where 100 is the 
        number of calculated features per sequence.
    """
    #uppercase protein sequence 
    sequence = sequence.upper()

    #raise value error if int cant be parsed from input lag
    try:
        lag = int(lag)
    except:
        raise ValueError("Invalid lag value input, integer cannot be parsed: {}.".format(lag))

    #validate lag, set default lag if invalid value input
    if (lag>=len(sequence) or (lag<0) or not (isinstance(lag, int))):
        lag = 30

    #validate weight, set default weight if invalid value input
    if ((weight<0) or not (isinstance(weight, float)) or not (isinstance(weight, int))):
        weight = 0.1

    #calculate quasi seq order using schneider distance matrix
    quasi_schneider = quasi_sequence_order(sequence, lag, weight=weight, 
        distance_matrix="schneider-wrede-physiochemical-distance-matrix.json")
    
    #calculate quasi seq order using schneider distance matrix
    quasi_grantham = quasi_sequence_order(sequence, lag, weight=weight, 
        distance_matrix="grantham-physiochemical-distance-matrix.json")

    #concat quasi sequence order dataframes
    quasi_sequence_order_all_df = pd.concat([quasi_schneider, quasi_grantham], axis=1)

    return quasi_sequence_order_all_df