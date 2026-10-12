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
[5] Kyte, J., & Doolittle, R. F. (1982). A simple method for displaying the
    hydropathic character of a protein. J Mol Biol, 157(1), 105-132.
    https://doi.org/10.1016/0022-2836(82)90515-0
[6] Gasteiger, E. et al. (2005). Protein identification and analysis tools on the
    ExPASy server. In The Proteomics Protocols Handbook (pp. 571-607). Humana
    Press. https://doi.org/10.1385/1-59259-890-0:571
[7] Wang, G. (2020). The antimicrobial peptide database provides a platform for
    decoding the design principles of naturally occurring antimicrobial peptides.
    Protein Sci, 29(1), 8-18. https://doi.org/10.1002/pro.3702
[8] Bachmair, A., Finley, D., & Varshavsky, A. (1986). In vivo half-life of a
    protein is a function of its amino-terminal residue. Science, 234(4773),
    179-186. https://doi.org/10.1126/science.3018930
[9] Gonda, D. K. et al. (1989). Universality and structure of the N-end rule.
    J Biol Chem, 264(28), 16700-16712. https://doi.org/10.1016/S0021-9258(19)84762-2
[10] Tobias, J. W. et al. (1991). The N-end rule in bacteria. Science, 254(5036),
    1374-1377. https://doi.org/10.1126/science.1962196
[11] Zeldovich, K. B., Berezovsky, I. N., & Shakhnovich, E. I. (2007). Protein and
    DNA sequence determinants of thermophilic adaptation. PLoS Comput Biol, 3(1),
    e5. https://doi.org/10.1371/journal.pcbi.0030005
[12] Suhre, K., & Claverie, J.-M. (2003). Genomic correlates of hyperthermostability,
    an update. J Biol Chem, 278(19), 17198-17202. https://doi.org/10.1074/jbc.M301327200
[13] Vihinen, M., Torkkila, E., & Riikonen, P. (1994). Accuracy of protein
    flexibility predictions. Proteins, 19(2), 141-149.
    https://doi.org/10.1002/prot.340190207
"""

import numpy as np
import pandas as pd

from ._constants import _validate_sequence, _validate_int
from ._helpers import _charge_vector, _to_df, _property_matrix, _sliding_windows, _profile_stats, _load_json_table
from .composition import charge_distribution, molecular_weight

__all__ = [
    "charge_fractions", "ncpr", "kappa", "sequence_charge_decoration", "charge_density",
    "hydrophobicity_profile", "hydrophobic_ratio", "half_life", "ivywrel_index", "cvp_bias", "flexibility",
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
        pandas DataFrame of shape 1 x 4, columns ``Charge_FCR``, ``Charge_FPlus``,
        ``Charge_FMinus`` and ``Charge_Sigma``.
    """
    #validate input protein sequence
    sequence = _validate_sequence(sequence)

    q = _charge_vector(sequence)
    f_plus, f_minus = np.mean(q > 0), np.mean(q < 0)
    fcr, ncpr_ = f_plus + f_minus, f_plus - f_minus
    sigma = ncpr_ ** 2 / fcr if fcr else 0.0
    #prefixed so "FCR" cannot clash with the tripeptide-composition column of the same name
    return _to_df([fcr, f_plus, f_minus, sigma], ["Charge_FCR", "Charge_FPlus", "Charge_FMinus", "Charge_Sigma"])

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
    gaps around them (left a, between b, right c).

    Blobs lying entirely in a neutral flank contribute (0 - sigma)^2 whatever the flank's
    length, so for a given b only flank lengths up to g - 1 are distinct. Trying every b
    with those a/c splits covers every distinct delta in O(n0 x g) candidates instead of
    the O(n0^2) of a full enumeration.
    """
    n0 = length - n_pos - n_neg
    candidates = set()
    for b in range(n0 + 1):
        rest = n0 - b
        for k in range(min(g, rest + 1)):
            candidates.add((k, b))
            candidates.add((rest - k, b))
    candidates = np.array(sorted(candidates))
    idx = np.arange(length)
    best = 0.0
    #chunk the (candidate x residue) arrays to bound memory on long sequences
    for chunk in np.array_split(candidates, max(1, len(candidates) // 1000)):
        a, b = chunk[:, :1], chunk[:, 1:]
        pos = ((idx >= a) & (idx < a + n_pos)).astype(float)
        neg = ((idx >= a + n_pos + b) & (idx < a + n_pos + b + n_neg)).astype(float)
        best = max(best, float(_delta(pos, neg, g, sigma).max()))
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
    distributed over the three surrounding gaps; every distinct arrangement is
    evaluated). kappa is the mean of kappa_g
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
        pandas DataFrame of shape 1 x 1, column ``Charge_SCD``.
    """
    #validate input protein sequence
    sequence = _validate_sequence(sequence)

    q = _charge_vector(sequence)
    idx = np.flatnonzero(q)
    #only charged residues contribute; pairwise over them is O(n_charged^2)
    m, n = np.triu_indices(len(idx), k=1)
    scd = np.sum(q[idx[m]] * q[idx[n]] * np.sqrt(idx[n] - idx[m])) / len(sequence)
    #prefixed so "SCD" cannot clash with the tripeptide-composition column of the same name
    return _to_df([scd], ["Charge_SCD"])

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

############################### Hydrophobicity Profile ####################################

def hydrophobicity_profile(sequence: str, scale: str = "kyte_doolittle", window: int = 9, edge: float = 1.0,
    threshold: float | None = None) -> pd.DataFrame:
    """
    Summarise a sliding-window hydrophobicity profile [5, 6]. Each residue is
    replaced by its scale value and a window of length ``window`` is slid along
    the sequence; ``edge < 1`` weights the window ends linearly down to ``edge``
    (ExPASy ProtScale semantics). The profile is reduced to its mean, max, min
    and standard deviation, plus the fraction of windows above ``threshold``
    when one is given (e.g. 1.6 for Kyte-Doolittle transmembrane segments).
    Windows longer than the sequence use the whole sequence.

    Parameters
    ==========
    :sequence: str
        protein sequence.
    :scale: str (default="kyte_doolittle")
        single-value scale: "kyte_doolittle", "hopp_woods", "eisenberg",
        "wimley_white", "engelman", "domcut", or any AAindex1 accession code.
    :window: int (default=9)
        window length.
    :edge: float (default=1.0)
        relative weight of the window ends, 0-1.
    :threshold: float | None (default=None)
        if set, add the fraction of windows with a value above it.

    Returns
    =======
    :profile_df: pd.DataFrame
        pandas DataFrame of shape 1 x 4 (1 x 5 with ``threshold``), columns
        ``HP_{scale}_w{window}_{mean,max,min,std[,frac_above]}``.
    """
    #validate input protein sequence
    sequence = _validate_sequence(sequence)
    window = _validate_int("window", window, 1)
    if not 0 <= edge <= 1:
        raise ValueError(f"edge must be between 0 and 1, got {edge}.")

    values = _property_matrix(sequence, scale)
    if values.shape[1] != 1:
        raise ValueError(f"Scale '{scale}' is multi-dimensional; use a single-value scale.")
    profile = _sliding_windows(values[:, 0], window, edge)
    stats = _profile_stats(profile, f"HP_{scale}_w{window}", threshold=threshold)
    return _to_df(list(stats.values()), list(stats))

################################ Hydrophobic Ratio ########################################

#hydrophobic residues as defined by the APD [7]
_apd_hydrophobic = frozenset("LIVFMACW")

def hydrophobic_ratio(sequence: str) -> pd.DataFrame:
    """
    Calculate the hydrophobic ratio of the Antimicrobial Peptide Database (APD)
    [7]: the percentage of residues that are L, I, V, F, M, A, C or W.

    HydrophobicRatio = N(LIVFMACW) / L x 100

    Parameters
    ==========
    :sequence: str
        protein sequence.

    Returns
    =======
    :hydrophobic_ratio_df: pd.DataFrame
        pandas DataFrame of shape 1 x 1, column ``HydrophobicRatio`` (percentage).
    """
    #validate input protein sequence
    sequence = _validate_sequence(sequence)
    return _to_df([sum(aa in _apd_hydrophobic for aa in sequence) / len(sequence) * 100], ["HydrophobicRatio"])

#################################### Half Life ############################################

def half_life(sequence: str) -> pd.DataFrame:
    """
    Estimate the in vivo half-life of a protein from its N-terminal residue by
    the N-end rule [8, 9, 10], as tabulated by ExPASy ProtParam [6]: mammalian
    reticulocytes (in vitro), yeast and E. coli (in vivo).

    Values are in hours. Lower bounds in the source (">20 h" in yeast, ">10 h"
    in E. coli, ">20 h" for mammalian Pro) are reported as the bound. The E. coli
    value for an N-terminal Pro is unknown in the source and returned as NaN.
    Not applicable to N-terminally modified proteins.

    Parameters
    ==========
    :sequence: str
        protein sequence.

    Returns
    =======
    :half_life_df: pd.DataFrame
        pandas DataFrame of shape 1 x 3, columns ``HalfLife_Mammal_h``,
        ``HalfLife_Yeast_h`` and ``HalfLife_Ecoli_h``.
    """
    #validate input protein sequence
    sequence = _validate_sequence(sequence)

    row = _load_json_table("residue_tables.json", "n_end_half_life")["values"][sequence[0]]
    values = [np.nan if row[k] is None else row[k] for k in ("mammal", "yeast", "ecoli")]
    return _to_df(values, ["HalfLife_Mammal_h", "HalfLife_Yeast_h", "HalfLife_Ecoli_h"])

################################## IVYWREL Index ##########################################

def ivywrel_index(sequence: str) -> pd.DataFrame:
    """
    Calculate the IVYWREL index of Zeldovich et al. [11]: the fraction of residues
    that are I, V, Y, W, R, E or L, which correlates with optimal growth
    temperature.

    IVYWREL = N(IVYWREL) / L

    Parameters
    ==========
    :sequence: str
        protein sequence.

    Returns
    =======
    :ivywrel_df: pd.DataFrame
        pandas DataFrame of shape 1 x 1, column ``IVYWREL`` (fraction).
    """
    #validate input protein sequence
    sequence = _validate_sequence(sequence)
    return _to_df([sum(aa in "IVYWREL" for aa in sequence) / len(sequence)], ["IVYWREL"])

#################################### CvP Bias #############################################

def cvp_bias(sequence: str) -> pd.DataFrame:
    """
    Calculate the charged-versus-polar (CvP) bias of Suhre & Claverie [12], an
    indicator of hyperthermostability:

    CvP = %(D + E + K + R) - %(N + Q + S + T)

    Parameters
    ==========
    :sequence: str
        protein sequence.

    Returns
    =======
    :cvp_df: pd.DataFrame
        pandas DataFrame of shape 1 x 1, column ``CvP_Bias`` (percentage points).
    """
    #validate input protein sequence
    sequence = _validate_sequence(sequence)
    charged = sum(aa in "DEKR" for aa in sequence)
    polar = sum(aa in "NQST" for aa in sequence)
    return _to_df([(charged - polar) / len(sequence) * 100], ["CvP_Bias"])

################################### Flexibility ###########################################

def flexibility(sequence: str, window: int = 9) -> pd.DataFrame:
    """
    Summarise the flexibility profile of Vihinen et al. [13]. Average normalised
    B-values (Vihinen's "average" parameter set, AAindex VINM940101) are averaged
    over a window of 9 whose end weights fall linearly to 0.25 (weights 0.25,
    0.4375, 0.625, 0.8125, 1, ..., 0.25; sum 5.25), and the profile is reduced
    to its mean, max, min and standard deviation.

    The neighbour-dependent parameter sets of [13] (residues with 0, 1 or 2 rigid
    neighbours) are not used. Compared with Biopython's
    ``ProteinAnalysis.flexibility``, which uses the same values and weights,
    protpy centres each window on its middle residue and includes the last
    window; Biopython takes the residue after the centre and omits the last window.

    Parameters
    ==========
    :sequence: str
        protein sequence.
    :window: int (default=9)
        window length (the parameters were optimised for 9).

    Returns
    =======
    :flexibility_df: pd.DataFrame
        pandas DataFrame of shape 1 x 4, columns ``Flex_{mean,max,min,std}``.
    """
    #validate input protein sequence
    sequence = _validate_sequence(sequence)
    window = _validate_int("window", window, 1)

    profile = _sliding_windows(_property_matrix(sequence, "vihinen_flexibility")[:, 0], window, edge=0.25)
    stats = _profile_stats(profile, "Flex")
    return _to_df(list(stats.values()), list(stats))
