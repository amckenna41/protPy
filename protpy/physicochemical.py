################################################################################
#############                  Physicochemical                    ##############
################################################################################

from __future__ import annotations

"""
Charge and charge-patterning descriptors (FCR, NCPR, kappa, SCD, charge
density).

Charges follow the Das & Pappu convention: K, R = +1; D, E = -1; H = 0.

References
==========
[1] Mao, A. H., Crick, S. L., Vitalis, A., Chicoine, C. L., & Pappu, R. V.
    (2010). Net charge per residue modulates conformational ensembles of
    intrinsically disordered proteins. PNAS, 107(18), 8183-8188.
    https://doi.org/10.1073/pnas.0911107107
[2] Das, R. K., & Pappu, R. V. (2013). Conformations of intrinsically disordered
    proteins are influenced by linear sequence distributions of oppositely
    charged residues. PNAS, 110(33), 13392-13397.
    https://doi.org/10.1073/pnas.1304749110
[3] Sawle, L., & Ghosh, K. (2015). A theoretical method to compute sequence
    dependent configurational properties in charged polymers and proteins.
    J Chem Phys, 143(8), 085101. https://doi.org/10.1063/1.4929391
[4] Müller, A. T., Gabernet, G., Hiss, J. A., & Schneider, G. (2017). modlAMP:
    Python for antimicrobial peptides. Bioinformatics, 33(17), 2753-2755.
    https://doi.org/10.1093/bioinformatics/btx285
"""

import numpy as np
import pandas as pd

from ._constants import _validate_sequence, _validate_int
from ._helpers import _charge_vector, _to_df
from .composition import charge_distribution, molecular_weight

__all__ = [
    "charge_fractions", "ncpr", "kappa", "sequence_charge_decoration", "charge_density",
]

################################## Charge Fractions #######################################

def charge_fractions(sequence: str) -> pd.DataFrame:
    """
    Calculate the charge fractions of a protein sequence [1, 2]:

    FCR    = (n+ + n-) / L          fraction of charged residues
    FPlus  = n+ / L,  FMinus = n- / L
    Sigma  = NCPR^2 / FCR           overall charge asymmetry (0 if FCR = 0)

    where n+ counts K/R and n- counts D/E (H is neutral). Values are fractions.

    Parameters
    ==========
    :sequence: str
        protein sequence.

    Returns
    =======
    :charge_fractions_df: pd.DataFrame
        pandas DataFrame of shape 1 x 4, columns ``FCR``, ``FPlus``, ``FMinus``, ``Sigma``.
    """
    #validate input protein sequence
    sequence = _validate_sequence(sequence)

    q = _charge_vector(sequence)
    f_plus, f_minus = np.mean(q > 0), np.mean(q < 0)
    fcr, ncpr_ = f_plus + f_minus, f_plus - f_minus
    sigma = ncpr_ ** 2 / fcr if fcr else 0.0
    return _to_df([fcr, f_plus, f_minus, sigma], ["FCR", "FPlus", "FMinus", "Sigma"])

#################################### NCPR #################################################

def ncpr(sequence: str) -> pd.DataFrame:
    """
    Calculate the Net Charge Per Residue (NCPR) of a protein sequence [1]:

    NCPR = (n+ - n-) / L

    with n+ counting K/R and n- counting D/E.

    Parameters
    ==========
    :sequence: str
        protein sequence.

    Returns
    =======
    :ncpr_df: pd.DataFrame
        pandas DataFrame of shape 1 x 1, column ``NCPR``.
    """
    #validate input protein sequence
    sequence = _validate_sequence(sequence)
    return _to_df([_charge_vector(sequence).mean()], ["NCPR"])

#################################### Kappa ################################################

def _blob_sigma(pos: np.ndarray, neg: np.ndarray, g: int) -> np.ndarray:
    """Charge asymmetry (f+ - f-)^2 / (f+ + f-) of every overlapping blob of size g (rows = candidates)."""
    kernel = np.ones(g)
    if pos.ndim == 1:
        p = np.convolve(pos, kernel, mode="valid")
        m = np.convolve(neg, kernel, mode="valid")
    else:
        #blob counts from prefix sums, one row per candidate sequence
        cp = np.concatenate([np.zeros((pos.shape[0], 1)), np.cumsum(pos, axis=1)], axis=1)
        cm = np.concatenate([np.zeros((neg.shape[0], 1)), np.cumsum(neg, axis=1)], axis=1)
        p, m = cp[:, g:] - cp[:, :-g], cm[:, g:] - cm[:, :-g]
    total = p + m
    return np.divide((p - m) ** 2, total * g, out=np.zeros_like(total, dtype=float), where=total > 0)


def _delta(pos: np.ndarray, neg: np.ndarray, g: int, sigma: float) -> np.ndarray:
    return np.mean((_blob_sigma(pos, neg, g) - sigma) ** 2, axis=-1)


def _delta_max(n_pos: int, n_neg: int, length: int, g: int, sigma: float) -> float:
    """
    Maximal delta over the most segregated arrangements of the composition: a block of
    positives and a block of negatives, with the neutral residues split across the three
    gaps around them. Every split is tried (O(n0^2) candidates), in chunks.
    """
    n0 = length - n_pos - n_neg
    best = 0.0
    idx = np.arange(length)
    for a in range(n0 + 1):
        #a neutrals, positive block, b neutrals, negative block, remaining neutrals
        b = np.arange(n0 - a + 1)[:, None]
        pos = (idx >= a) & (idx < a + n_pos)
        neg = (idx >= a + n_pos + b) & (idx < a + n_pos + b + n_neg)
        pos = np.broadcast_to(pos, neg.shape).astype(float)
        best = max(best, float(_delta(pos, neg.astype(float), g, sigma).max()))
    return best


def kappa(sequence: str, blob_sizes: tuple[int, ...] = (5, 6)) -> pd.DataFrame:
    """
    Calculate the kappa charge-patterning parameter of Das & Pappu [2]. For blob
    size g, the sequence is split into the L - g + 1 overlapping blobs and

    sigma_b = (f+ - f-)^2 / (f+ + f-)        per blob (0 if uncharged)
    delta_g = mean over blobs of (sigma_b - sigma)^2
    kappa_g = delta_g / delta_max,g

    where sigma = NCPR^2 / FCR of the whole sequence and delta_max,g is the
    largest delta over the maximally segregated arrangements of the same
    composition (a positive block and a negative block with the neutral residues
    distributed over the three surrounding gaps). kappa is the mean of kappa_g
    over ``blob_sizes`` (5 and 6 in [2]); 0 means well mixed and 1 fully
    segregated. Charges: K, R = +1; D, E = -1; H = 0.

    kappa is reported as 0 (not NaN) when the sequence has no charges or only one
    charge type, where oppositely charged residues cannot be segregated. A blob
    size longer than the sequence is reduced to the sequence length.

    Parameters
    ==========
    :sequence: str
        protein sequence.
    :blob_sizes: tuple[int, ...] (default=(5, 6))
        blob sizes g to average over.

    Returns
    =======
    :kappa_df: pd.DataFrame
        pandas DataFrame of shape 1 x (1 + len(blob_sizes)), columns ``Kappa`` and
        ``Kappa_{g}``.
    """
    #validate input protein sequence
    sequence = _validate_sequence(sequence)
    blob_sizes = [_validate_int("blob size", g, 1) for g in blob_sizes]
    if not blob_sizes:
        raise ValueError("At least one blob size must be given.")

    q = _charge_vector(sequence)
    pos, neg = (q > 0).astype(float), (q < 0).astype(float)
    n_pos, n_neg, L = int(pos.sum()), int(neg.sum()), len(sequence)
    values = []
    for g in blob_sizes:
        if n_pos == 0 or n_neg == 0:
            values.append(0.0)
            continue
        g = min(g, L)
        sigma = (n_pos - n_neg) ** 2 / (L * (n_pos + n_neg))
        d_max = _delta_max(n_pos, n_neg, L, g, sigma)
        values.append(min(float(_delta(pos, neg, g, sigma)) / d_max, 1.0) if d_max else 0.0)
    return _to_df([float(np.mean(values))] + values, ["Kappa"] + [f"Kappa_{g}" for g in blob_sizes])

########################### Sequence Charge Decoration ####################################

def sequence_charge_decoration(sequence: str) -> pd.DataFrame:
    """
    Calculate the Sequence Charge Decoration (SCD) of Sawle & Ghosh [3]:

    SCD = (1/L) sum_{m=2}^{L} sum_{n=1}^{m-1} q_m q_n (m - n)^{1/2}

    with q = +1 for K/R, -1 for D/E and 0 otherwise. More negative values
    indicate more segregated opposite charges; uncharged sequences give 0.

    Parameters
    ==========
    :sequence: str
        protein sequence.

    Returns
    =======
    :scd_df: pd.DataFrame
        pandas DataFrame of shape 1 x 1, column ``SCD``.
    """
    #validate input protein sequence
    sequence = _validate_sequence(sequence)

    q = _charge_vector(sequence)
    idx = np.flatnonzero(q)
    #only charged residues contribute; pairwise over them is O(n_charged^2)
    m, n = np.triu_indices(len(idx), k=1)
    scd = np.sum(q[idx[m]] * q[idx[n]] * np.sqrt(idx[n] - idx[m])) / len(sequence)
    return _to_df([scd], ["SCD"])

#################################### Charge Density #######################################

def charge_density(sequence: str, ph: float = 7.4) -> pd.DataFrame:
    """
    Calculate the charge density of a protein sequence [4]: the net charge at
    ``ph`` (from ``charge_distribution``) divided by the molecular weight (from
    ``molecular_weight``), in charges per Da.

    ChargeDensity = NetCharge(pH) / MW

    Parameters
    ==========
    :sequence: str
        protein sequence.
    :ph: float (default=7.4)
        pH at which the net charge is evaluated.

    Returns
    =======
    :charge_density_df: pd.DataFrame
        pandas DataFrame of shape 1 x 1, column ``ChargeDensity``. The value is
        rounded to 6 dp because it is typically of order 1e-4.
    """
    #validate input protein sequence
    sequence = _validate_sequence(sequence)

    net = charge_distribution(sequence, ph=ph)["NetCharge"].iloc[0]
    mw = molecular_weight(sequence)["MolecularWeight"].iloc[0]
    return pd.DataFrame([[round(net / mw, 6)]], columns=["ChargeDensity"])
