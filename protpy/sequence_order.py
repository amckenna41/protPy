################################################################################
###############              Quasi Sequence Order                ###############
################################################################################

from __future__ import annotations

import os
import json
import functools
from difflib import get_close_matches
import pandas as pd
import math
import sys
from . import composition

"""
References
==========
[1] Kuo-Chen Chou. Prediction of Protein Subcellar Locations by Incorporating
    Quasi-Sequence-Order Effect. Biochemical and Biophysical Research Communications,
    2000, 278, 477-483.
[2] Kuo-Chen Chou and Yu-Dong Cai. Prediction of Protein Subcellular Locations by
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
from ._constants import amino_acids, _validate_sequence, _validate_int, _validate_weight
from ._helpers import _minmax

@functools.lru_cache(maxsize=2)
def _load_distance_matrix(name: str) -> dict:
    """Load and cache a physicochemical distance matrix JSON by name."""
    path = os.path.join(
        os.path.dirname(os.path.abspath(sys.modules['protpy'].__file__)),
        "data", f"{name}-physicochemical-distance-matrix.json"
    )
    try:
        with open(path) as f:
            return json.load(f)
    except (OSError, json.JSONDecodeError):
        raise ValueError(f"Error parsing distance matrix file: {path}.")

############################# Sequence Order Coupling Number ##############################

def sequence_order_coupling_number_(sequence: str, d: int = 1, distance_matrix: str = "schneider-wrede",
    normalized: bool = False) -> float:
    """
    Calculate Sequence Order Coupling Number (SOCN) features for input protein sequence.
    SOCN computes the dissimilarity between amino acid pairs. The distance between 
    amino acid pairs is determined by d which varies between 1 to lag. For each d, it 
    computes the sum of the dissimilarities of all amino acid pairs. The output will be 
    a single float value representing the SOCN [1].
    
    This function should not be confused with sequence_order_coupling_number() below 
    which calculates the multiple SOCN descriptor values for the input sequence 
    according to the value of lag; with the default of lag=30, 30 SOCN's will be 
    calculated.

    Parameters
    ==========
    :sequence: str
        protein sequence.
    :d: int (default=1)
        gap between two amino acids.
    :distance_matrix: str (default="schneider-wrede")
        physicochemical distance matrix to use; accepts "schneider-wrede" or "grantham".
    :normalized: bool (default=False)
        rescale to [0, 1] by the maximum possible value (N - d) x max(distance)^2
        (tau sums N - d squared distances); 0 when N <= d. Rounded to 6 dp.

    Returns
    =======
    :socn: float
        calculated sequence order coupling number value from sequence.
    """
    #validate input protein sequence
    sequence = _validate_sequence(sequence)

    #validate input distance matrix is valid, get its closest match using closeness function
    dist_matrix_matches = get_close_matches(distance_matrix, ["schneider-wrede", "grantham"], cutoff=0.6)
    if (dist_matrix_matches != []):
        distance_matrix = dist_matrix_matches[0]  #set desc to closest descriptor match found
    else:
        raise ValueError(f"Could not find a match for the input distance matrix {distance_matrix} in"
            f" list of available valid distance matrices:\n{['schneider-wrede', 'grantham']}.")
    
    #load (or retrieve cached) distance matrix
    distance_matrix_data = _load_distance_matrix(distance_matrix)

    tau = 0.0

    #iterate over length of sequence minus gap, incrementing the SOCN at each iteration using
    #the values from the distance matrix with the current and next amino acid
    for aa in range(len(sequence) - d):
        current_aa = sequence[aa]
        next_aa = sequence[aa + d]
        tau = tau + math.pow(distance_matrix_data[current_aa + next_aa], 2)

    if normalized:
        bound = (len(sequence) - d) * max(distance_matrix_data.values()) ** 2
        return round(tau / bound, 6) if bound > 0 else 0.0
    return round(tau, 3)

def sequence_order_coupling_number(sequence: str, lag: int = 30,
    distance_matrix: str = "schneider-wrede", normalized: bool = False) -> pd.DataFrame:
    """
    Calculate Sequence Order Coupling Number (SOCN) features for input protein sequence.
    SOCN computes the dissimilarity between amino acid pairs. The distance between 
    amino acid pairs is determined by d which varies between 1 to lag. For each d, it 
    computes the sum of the dissimilarities of all amino acid pairs. The number of 
    output features can be calculated as N, where N = lag, by default this value 
    is 30 which generates an output of 1 x 30.

    Parameters
    ==========
    :sequence: str
        protein sequence.
    :lag: int (default=30)
        maximum gap between 2 amino acids; the protein length should be larger than lag.
    :distance_matrix: str (default="schneider-wrede")
        physicochemical distance matrix to use; accepts "schneider-wrede" or "grantham".
    :normalized: bool (default=False)
        rescale every SOCN to [0, 1] as in ``sequence_order_coupling_number_``;
        columns get a ``_norm`` suffix.

    Returns
    =======
    :sequence_order_df: pd.Dataframe
        Dataframe of SOCN descriptor values for all protein sequences. Output
        will be of the shape 1 x N, where N is the number of features calculated 
        from the descriptor and N=lag.
    """
    #validate input protein sequence
    sequence = _validate_sequence(sequence)

    #raise value error if a non-negative int cant be parsed from input lag
    lag = _validate_int("lag", lag)

    #set dataframe column name prefixes based on input distance matrix, raise error if invalid name input
    col_prefix = ""
    if (distance_matrix == "schneider-wrede"):
        col_prefix = "SOCN_SW"
    elif (distance_matrix == "grantham"):
        col_prefix = "SOCN_Grant"
    else:
        raise ValueError(f"Invalid distance matrix name input {distance_matrix}, values can be schneider-wrede or grantham.")

    #validate lag, set default lag of 30 if invalid value input
    if (lag>=len(sequence)):
        lag = 30

    seq_order = {}

    #iterate over sequence with lag, calculating SOCN using input distance matrix values
    for i in range(lag):
        tau = sequence_order_coupling_number_(sequence, i+1, distance_matrix, normalized=normalized)

        #append SOCN column and value to dict
        seq_order[col_prefix + str(i + 1) + ("_norm" if normalized else "")] = tau if normalized else round(tau, 3)

    #transform SOCN data into pandas dataframe
    sequence_order_df = pd.DataFrame([list(seq_order.values())], columns=list(seq_order.keys()))

    return sequence_order_df

def sequence_order_coupling_number_all(sequence: str, lag: int = 30, normalized: bool = False) -> pd.DataFrame:
    """
    Calculate Sequence Order Coupling Number (SOCN) descriptor values of input protein 
    sequence using both matrices (schneider-wrede & grantham) [3]. The distance between 
    amino acid pairs is determined by d which varies between 1 to lag. For each d, it 
    computes the sum of the dissimilarities of all amino acid pairs. Each matrix generates
    an output of 1 x N, where N is the lag, so using the two concatenated matrices, the 
    output will be 1 x (N * 2). With a default lag of 30, this will generate an output 
    of 1 x 60.

    Parameters
    ==========
    :sequence: str
        protein sequence.
    :lag: int (default=30)
        maximum gap between 2 amino acids; the protein length should be larger than lag.
    :normalized: bool (default=False)
        rescale every SOCN to [0, 1] (see ``sequence_order_coupling_number_``);
        columns get a ``_norm`` suffix.

    Returns
    =======
    :socn_all: pd.Dataframe
        Concatenated dataframe of SOCN descriptor values of input protein sequence using
        both distance matrices. The number of output features can be calculated as N * 2, 
        where N = lag, by default this value is 30 which generates an output of 1 x 30 
        for each matrix, with the two matrices the output will be 1 x (30*2), using the 
        default lag.
    """
    #validate input protein sequence
    sequence = _validate_sequence(sequence)

    #raise value error if a non-negative int cant be parsed from input lag
    lag = _validate_int("lag", lag)

    #validate lag, set default lag of 30 if invalid value input
    if (lag>=len(sequence)):
        lag = 30

    #calculate SOCN using schneider-wrede distance matrix
    socn_schneider = sequence_order_coupling_number(sequence, lag, "schneider-wrede", normalized)
    
    #calculate SOCN using grantham distance matrix
    socn_grantham = sequence_order_coupling_number(sequence, lag, "grantham", normalized)

    #merge 2 dataframes 
    socn_all = pd.concat([socn_schneider, socn_grantham], axis=1)

    return socn_all

################################## Quasi Sequence Order ###################################

def quasi_sequence_order(sequence: str, lag: int = 30, weight: float = 0.1,
    distance_matrix: str = "schneider-wrede", normalized: bool = False) -> pd.DataFrame:
    """
    Calculate Quasi Sequence Order (QSO) features for the protein sequences.
    The quasi-sequence-order descriptors were proposed by K.C. Chou, et.al. [1].
    They are derived from the distance matrix between the 20 amino acids. By default,
    the Scheider-Wrede physicochemical distance matrix was used. Also utilised in
    the descriptor calculation is the Grantham chemical distance matrix. Both of
    these matrices are used by Grantham et. al. in the calculation of the descriptor 
    [3, 4]. N + 20 values are calculated per sequence, where N is the lag, with a 
    default lag of 30, the output will be 1 x 50. There is also a weighting factor
    that can be assigned to determine that weight per amino acid.

    Parameters
    ==========
    :sequence: str
        protein sequence.
    :lag: int (default=30)
        maximum gap between 2 amino acids; the protein length should be larger than lag.
    :weight: float (default = 0.1)
        weighting factor.
    :distance_matrix: str (default="schneider-wrede")
        physicochemical distance matrix to use; accepts "schneider-wrede" or "grantham".
    :normalized: bool (default=False)
        rescale to [0, 1]: the first 20 components (amino-acid percentages over
        1 + w x sum(tau), at most 100) are divided by 100 and the lag components
        (w x tau_d / (1 + w x sum(tau)), always below 1) are kept; columns get a
        ``_norm`` suffix.

    Returns
    =======
    :quasi_sequence_order_df: pd.Dataframe
        dataframe of quasi-sequence-order descriptor values for the protein sequences,
        with output shape 1 x (N + 20), where N is the lag. With a default
        lag of 30 the output will be 1 x 50 per sequence.
    """
    #validate input protein sequence
    sequence = _validate_sequence(sequence)

    #raise value error if a non-negative int can't be parsed from input lag
    lag = _validate_int("lag", lag)

    #validate lag, set default lag of 30 if invalid value input
    if (lag>=len(sequence)):
        lag = 30

    #validate weight
    weight = _validate_weight(weight)

    #set dataframe column name prefixes based on input distance matrix, raise error if invalid name input
    col_prefix = ""
    if (distance_matrix == "schneider-wrede"):
        col_prefix = "QSO_SW"
    elif (distance_matrix == "grantham"):
        col_prefix = "QSO_Grant"
    else:
        raise ValueError(f"Invalid distance matrix name input {distance_matrix}, values can be schneider-wrede or grantham.")
    
    quasi_sequence_order = {}

    #calculate SOCN for each gap up to lag once, using the specified distance matrix
    _right_part = [sequence_order_coupling_number_(sequence, i+1, distance_matrix) for i in range(lag)]

    #get amino acid composition
    aa_comp = composition.amino_acid_composition(sequence)

    #iterate over amino acids, calculating first 20 descriptor values
    temp = 1 + weight * sum(_right_part)
    for index, aa in enumerate(amino_acids):
        quasi_sequence_order[col_prefix + str(index + 1)] = round(aa_comp[aa].values[0] / temp, 6)

    #calculate the last lag descriptor values
    for index in range(20, 20 + lag):
        quasi_sequence_order[col_prefix + str(index + 1)] = round(
            weight * _right_part[index - 20] / temp, 6)
    
    #transform descriptor data into pandas dataframe
    quasi_sequence_order_df = pd.DataFrame([list(quasi_sequence_order.values())], 
        columns=list(quasi_sequence_order.keys()))
    
    if normalized:
        #first 20 components: amino-acid percentages / (1 + w sum(tau)) <= 100;
        #lag components: w tau_d / (1 + w sum(tau)) < 1
        return _minmax(quasi_sequence_order_df, 0.0, [100.0] * 20 + [1.0] * (quasi_sequence_order_df.shape[1] - 20))
    return quasi_sequence_order_df

def quasi_sequence_order_all(sequence: str, lag: int = 30, weight: float = 0.1, normalized: bool = False) -> pd.DataFrame:
    """
    Calculate Quasi Sequence Order features for input protein sequence using 
    both physicochemical distance matrices. Concatenate into one output dataframe. 
    The output will be in the shape 1 x ((N + 20)*2), where ((N + 20)*2) is the 
    quasi sequence order output from one matrix and N is the lag. The output  
    is multiplied by two to take into account the 2 matrices being concatenated. 
    There is also a weighting factor that can be assigned to determine that 
    weight per amino acid.
    
    Parameters
    ==========
    :sequence: str
        protein sequence.
    :lag: int (default=30)
        maximum gap between 2 amino acids; the protein length should be larger than lag.
    :weight: float (default=0.1)
        weighting factor.
    :normalized: bool (default=False)
        see ``quasi_sequence_order``; columns get a ``_norm`` suffix.

    Returns
    =======
    :quasi_sequence_order_all_df: pd.Dataframe
        dataframe of quasi-sequence-order descriptor values for the
        protein sequences, with output shape 1 x ((N + 20)*2) where ((N + 20)*2) is the
        quasi sequence order output from one matrix and N is the lag. The output
        is multiplied by two to take into account the 2 matrices being concatenated. 
    """
    #validate input protein sequence
    sequence = _validate_sequence(sequence)

    #raise value error if a non-negative int can't be parsed from input lag
    lag = _validate_int("lag", lag)

    #validate lag, set default lag of 30 if invalid value input
    if (lag>=len(sequence)):
        lag = 30

    #validate weight
    weight = _validate_weight(weight)

    #calculate quasi seq order using schneider distance matrix
    quasi_schneider = quasi_sequence_order(sequence, lag, weight=weight, 
        distance_matrix="schneider-wrede", normalized=normalized)
    
    #calculate quasi seq order using schneider distance matrix
    quasi_grantham = quasi_sequence_order(sequence, lag, weight=weight, 
        distance_matrix="grantham", normalized=normalized)

    #concat quasi sequence order dataframes
    quasi_sequence_order_all_df = pd.concat([quasi_schneider, quasi_grantham], axis=1)

    return quasi_sequence_order_all_df