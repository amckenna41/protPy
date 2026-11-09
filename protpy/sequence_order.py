################################################################################
###############              Quasi Sequence Order                ###############
################################################################################

from __future__ import annotations

import os
import json
import functools

import numpy as np
import pandas as pd

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
from ._helpers import _check_lag, _to_df, _warn_deprecated

#sequence_order_coupling_number_ (a single SOCN lag) is internal and not exported
__all__ = ["sequence_order_coupling_number", "quasi_sequence_order",
    #deprecated, removed in protpy 2.0
    "sequence_order_coupling_number_all", "quasi_sequence_order_all"]

#column prefixes of the two distance matrices
_matrices = {"schneider-wrede": "SW", "grantham": "Grant"}

@functools.lru_cache(maxsize=2)
def _load_distance_matrix(name: str) -> dict:
    """Load and cache a physicochemical distance matrix JSON by name."""
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data",
                        f"{name}-physicochemical-distance-matrix.json")
    try:
        with open(path) as f:
            return json.load(f)
    except (OSError, json.JSONDecodeError):
        raise ValueError(f"Error parsing distance matrix file: {path}.")

def _matrix_names(distance_matrix: str) -> list[str]:
    """
    Get the distance matrices to use: "schneider-wrede", "grantham" or "both".

    Parameters
    ==========
    :distance_matrix: str
        distance matrix name or "both".

    Returns
    =======
    :names: list[str]
        names of the distance matrices.
    """
    if distance_matrix == "both":
        return list(_matrices)
    if distance_matrix not in _matrices:
        raise ValueError(f"Invalid distance matrix name input {distance_matrix}, values can be "
                         "schneider-wrede, grantham or both.")
    return [distance_matrix]

def _socn(sequence: str, lag: int, distance_matrix: str) -> np.ndarray:
    """
    Calculate the sequence order coupling numbers tau_1..tau_lag of a validated
    protein sequence for one distance matrix [1]:

    tau_d = sum_{i=1}^{L-d} dist(s_i, s_i+d)^2

    Parameters
    ==========
    :sequence: str
        validated protein sequence.
    :lag: int
        maximum gap d.
    :distance_matrix: str
        "schneider-wrede" or "grantham".

    Returns
    =======
    :tau: np.ndarray
        coupling numbers for d = 1..lag.
    """
    data = _load_distance_matrix(distance_matrix)
    matrix = np.array([[data[a + b] for b in amino_acids] for a in amino_acids]) ** 2
    index = {aa: i for i, aa in enumerate(amino_acids)}
    idx = np.array([index[aa] for aa in sequence])
    return np.array([matrix[idx[:-d], idx[d:]].sum() for d in range(1, lag + 1)])

############################# Sequence Order Coupling Number ##############################

def sequence_order_coupling_number_(sequence: str, d: int = 1, distance_matrix: str = "schneider-wrede") -> float:
    """
    Calculate the Sequence Order Coupling Number (SOCN) of a protein sequence for
    a single gap d [1].

    Deprecated since v1.6.0, removed in 2.0: protpy no longer uses it and it is
    not exported from the protpy namespace; use
    ``sequence_order_coupling_number(sequence, lag=d).iloc[0, -1]``.

    Parameters
    ==========
    :sequence: str
        protein sequence.
    :d: int (default=1)
        gap between two amino acids.
    :distance_matrix: str (default="schneider-wrede")
        physicochemical distance matrix to use; accepts "schneider-wrede" or "grantham".

    Returns
    =======
    :socn: float
        sequence order coupling number for gap d (0 when there are no pairs).
    """
    _warn_deprecated("sequence_order_coupling_number_", "sequence_order_coupling_number")
    #validate input protein sequence
    sequence = _validate_sequence(sequence)
    _matrix_names(distance_matrix)
    d = _validate_int("d", d, 1)
    return float(_socn(sequence, d, distance_matrix)[-1]) if d < len(sequence) else 0.0

def sequence_order_coupling_number(sequence: str, lag: int = 30,
    distance_matrix: str = "schneider-wrede", normalize: bool = False) -> pd.DataFrame:
    """
    Calculate Sequence Order Coupling Number (SOCN) features for input protein sequence.
    SOCN computes the dissimilarity between amino acid pairs separated by d
    residues, for each d = 1..lag [1]:

    tau_d = sum_{i=1}^{L-d} dist(s_i, s_i+d)^2

    where dist is the Schneider-Wrede [3] or Grantham [4] physicochemical
    distance between two amino acids and L is the length of the sequence. As tau_d
    is a sum over L - d pairs, on proteins of different lengths it is mostly a
    measure of length; ``normalize=True`` divides it by L - d (the mean squared
    distance), which removes this dependence. ``distance_matrix="both"`` returns
    the Schneider-Wrede and Grantham values side by side (formerly
    sequence_order_coupling_number_all).

    Parameters
    ==========
    :sequence: str
        protein sequence.
    :lag: int (default=30)
        maximum gap between 2 amino acids, less than the length of the sequence.
    :distance_matrix: str (default="schneider-wrede")
        physicochemical distance matrix: "schneider-wrede", "grantham" or "both".
    :normalize: bool (default=False)
        divide tau_d by the number of pairs, L - d.

    Returns
    =======
    :sequence_order_df: pd.Dataframe
        Dataframe of SOCN descriptor values, of the shape 1 x lag (1 x 2*lag for
        "both"). Column names are SOCN_SW<d> and SOCN_Grant<d>.
    """
    #validate input protein sequence
    sequence = _validate_sequence(sequence)
    names = _matrix_names(distance_matrix)
    lag = _check_lag("lag", lag, len(sequence))

    values, columns = [], []
    for name in names:
        tau = _socn(sequence, lag, name)
        if normalize:
            tau = tau / (len(sequence) - np.arange(1, lag + 1))
        values += list(tau)
        columns += [f"SOCN_{_matrices[name]}{d}" for d in range(1, lag + 1)]
    return _to_df(values, columns)

def sequence_order_coupling_number_all(sequence: str, lag: int = 30) -> pd.DataFrame:
    """
    Calculate Sequence Order Coupling Number (SOCN) descriptor values of input protein
    sequence using both distance matrices (schneider-wrede & grantham) [3, 4], of
    the shape 1 x (lag * 2).

    Deprecated since v1.6.0, removed in 2.0: use
    ``sequence_order_coupling_number(sequence, lag, distance_matrix="both")``.

    Parameters
    ==========
    :sequence: str
        protein sequence.
    :lag: int (default=30)
        maximum gap between 2 amino acids, less than the length of the sequence.

    Returns
    =======
    :socn_all: pd.Dataframe
        dataframe of SOCN values for both distance matrices, by default 1 x 60.
    """
    _warn_deprecated("sequence_order_coupling_number_all",
                     'sequence_order_coupling_number(sequence, lag, distance_matrix="both")')
    return sequence_order_coupling_number(sequence, lag, distance_matrix="both")

################################## Quasi Sequence Order ###################################

def quasi_sequence_order(sequence: str, lag: int = 30, weight: float = 0.1,
    distance_matrix: str = "schneider-wrede") -> pd.DataFrame:
    """
    Calculate Quasi Sequence Order (QSO) features for the protein sequences [1, 2].
    The quasi-sequence-order descriptors combine the amino acid composition with
    the sequence order coupling numbers tau_d (see
    sequence_order_coupling_number) of a physicochemical distance matrix, the
    Schneider-Wrede [3] or Grantham [4] matrix:

    X_r = f_r / (1 + w * sum_{d=1}^{lag} tau_d)          for the 20 amino acids

    X_{20+d} = w * tau_d / (1 + w * sum_{d=1}^{lag} tau_d)   for d = 1..lag

    where f_r is the fraction of amino acid r and w the weighting factor. The
    20 + lag values sum to 1. (Before v1.6.0 f_r was a percentage, so the first 20
    values were 100 times too large.) ``distance_matrix="both"`` returns the
    values for both matrices side by side (formerly quasi_sequence_order_all);
    the default output is 1 x 50.

    Parameters
    ==========
    :sequence: str
        protein sequence.
    :lag: int (default=30)
        maximum gap between 2 amino acids, less than the length of the sequence.
    :weight: float (default = 0.1)
        weighting factor.
    :distance_matrix: str (default="schneider-wrede")
        physicochemical distance matrix: "schneider-wrede", "grantham" or "both".

    Returns
    =======
    :quasi_sequence_order_df: pd.Dataframe
        dataframe of quasi-sequence-order descriptor values, of the shape
        1 x (lag + 20) (1 x 2*(lag + 20) for "both"). Column names are
        QSO_SW<n> and QSO_Grant<n>.
    """
    #validate input protein sequence
    sequence = _validate_sequence(sequence)
    names = _matrix_names(distance_matrix)
    lag = _check_lag("lag", lag, len(sequence))
    weight = _validate_weight(weight)

    comp = np.array([sequence.count(aa) for aa in amino_acids]) / len(sequence)
    values, columns = [], []
    for name in names:
        tau = _socn(sequence, lag, name)
        denominator = 1 + weight * tau.sum()
        values += list(np.concatenate([comp, weight * tau]) / denominator)
        columns += [f"QSO_{_matrices[name]}{n}" for n in range(1, 21 + lag)]
    return _to_df(values, columns)

def quasi_sequence_order_all(sequence: str, lag: int = 30, weight: float = 0.1) -> pd.DataFrame:
    """
    Calculate Quasi Sequence Order features for input protein sequence using
    both physicochemical distance matrices, of the shape 1 x ((lag + 20) * 2).

    Deprecated since v1.6.0, removed in 2.0: use
    ``quasi_sequence_order(sequence, lag, weight, distance_matrix="both")``.

    Parameters
    ==========
    :sequence: str
        protein sequence.
    :lag: int (default=30)
        maximum gap between 2 amino acids, less than the length of the sequence.
    :weight: float (default=0.1)
        weighting factor.

    Returns
    =======
    :quasi_sequence_order_all_df: pd.Dataframe
        dataframe of QSO values for both distance matrices, by default 1 x 100.
    """
    _warn_deprecated("quasi_sequence_order_all", 'quasi_sequence_order(sequence, lag, weight, distance_matrix="both")')
    return quasi_sequence_order(sequence, lag, weight, distance_matrix="both")
