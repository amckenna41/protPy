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
    Calculate the Charge Fractions (FCR) of protein sequence. The fraction of
    charged residues (FCR), the fractions of positive and negative residues and
    the overall charge asymmetry (sigma) describe the charge content of a protein
    sequence, and are key parameters of the conformations of intrinsically
    disordered proteins [1, 2]. They are calculated as:

    FCR = (n+ + n-) / L

    FPlus = n+ / L, FMinus = n- / L

    Sigma = NCPR^2 / FCR

    where n+ is the number of positive residues (K, R), n- the number of negative
    residues (D, E), L is the length of the sequence and NCPR = (n+ - n-) / L. His
    is treated as neutral and sigma is 0 when there are no charged residues.

    Parameters
    ==========
    :sequence: str
        protein sequence.

    Returns
    =======
    :charge_fractions_df: pd.DataFrame
        pandas dataframe of charge fractions for protein sequence. Dataframe will
        be of the shape 1 x 4, where 4 is the number of features calculated from
        the descriptor and 1 is the input sequence. Column names are Charge_FCR,
        Charge_FPlus, Charge_FMinus and Charge_Sigma.
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
    Calculate the Net Charge Per Residue (NCPR) of protein sequence. NCPR is the
    net charge of a protein sequence divided by its length [1], calculated as:

    NCPR = (n+ - n-) / L

    where n+ is the number of positive residues (K, R), n- the number of negative
    residues (D, E) and L is the length of the sequence.

    Parameters
    ==========
    :sequence: str
        protein sequence.

    Returns
    =======
    :ncpr_df: pd.DataFrame
        pandas dataframe of NCPR for protein sequence. Dataframe will be of the
        shape 1 x 1, with the column NCPR.
    """
    #validate input protein sequence
    sequence = _validate_sequence(sequence)
    return _to_df([_charge_vector(sequence).mean()], ["NCPR"])

#################################### Kappa ################################################

def _blob_sigma(pos: np.ndarray, neg: np.ndarray, g: int) -> np.ndarray:
    """
    Calculate the charge asymmetry of every overlapping blob (window) of size g
    along one or more sequences, as used in the kappa calculation [2], where:

    sigma_b = (f+ - f-)^2 / (f+ + f-)

    and f+ and f- are the fractions of positive and negative residues in the blob.
    Uncharged blobs have a charge asymmetry of 0.

    Parameters
    ==========
    :pos: np.ndarray
        1 for positive residues and 0 otherwise; one row per sequence when 2-D.
    :neg: np.ndarray
        1 for negative residues and 0 otherwise, same shape as pos.
    :g: int
        blob size.

    Returns
    =======
    :blob_sigma: np.ndarray
        charge asymmetry of each blob, with one row per sequence when 2-D.
    """
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
    """
    Calculate the mean squared deviation (delta) of the blob charge asymmetries
    from the charge asymmetry of the whole sequence, as used in the kappa
    calculation [2].

    Parameters
    ==========
    :pos: np.ndarray
        1 for positive residues and 0 otherwise; one row per sequence when 2-D.
    :neg: np.ndarray
        1 for negative residues and 0 otherwise, same shape as pos.
    :g: int
        blob size.
    :sigma: float
        charge asymmetry of the whole sequence.

    Returns
    =======
    :delta: np.ndarray
        delta of each sequence (a single value for a 1-D input).
    """
    return np.mean((_blob_sigma(pos, neg, g) - sigma) ** 2, axis=-1)


def _delta_max(n_pos: int, n_neg: int, length: int, g: int, sigma: float) -> float:
    """
    Calculate the maximum delta of a sequence composition, used to normalise kappa
    [2]. The maximum is taken over the most segregated arrangements of the
    composition: a block of positive residues and a block of negative residues,
    with the neutral residues split across the three gaps around them. Blobs lying
    entirely in a neutral flank contribute the same value whatever the length of
    the flank, so for each gap between the blocks only flank lengths up to g - 1
    need to be tried. This evaluates every distinct arrangement in O(n0 * g)
    candidates, where n0 is the number of neutral residues.

    Parameters
    ==========
    :n_pos: int
        number of positive residues.
    :n_neg: int
        number of negative residues.
    :length: int
        length of the sequence.
    :g: int
        blob size.
    :sigma: float
        charge asymmetry of the whole sequence.

    Returns
    =======
    :best: float
        maximum delta over the segregated arrangements.
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
    Calculate the Kappa (Kappa) charge patterning parameter of protein sequence.
    Kappa describes how well oppositely charged residues are mixed or segregated
    along a protein sequence, which strongly affects the conformations of
    intrinsically disordered proteins [2]. The sequence is split into overlapping
    blobs of g residues and kappa is calculated as:

    sigma_b = (f+ - f-)^2 / (f+ + f-)

    delta_g = mean over all blobs of (sigma_b - sigma)^2

    kappa_g = delta_g / delta_max,g

    where sigma_b is the charge asymmetry of a blob, f+ and f- are its fractions of
    positive (K, R) and negative (D, E) residues, sigma = NCPR^2 / FCR is the charge
    asymmetry of the whole sequence and delta_max,g is the largest delta over the
    most segregated arrangements of the same composition. The final kappa is the
    mean of kappa_g over the blob sizes (5 and 6 by default). Kappa ranges from 0
    (well mixed) to 1 (fully segregated). It is reported as 0 when the sequence has
    no charged residues or only one type of charge, and blobs longer than the
    sequence are reduced to its length.

    Parameters
    ==========
    :sequence: str
        protein sequence.
    :blob_sizes: tuple (default=(5, 6))
        blob sizes g to average kappa over.

    Returns
    =======
    :kappa_df: pd.DataFrame
        pandas dataframe of kappa values for protein sequence. Dataframe will be of
        the shape 1 x (1 + number of blob sizes), where the first column (Kappa) is
        the mean and the remaining columns (Kappa_g) are the values for each blob
        size. By default the shape will be 1 x 3.
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
    Calculate the Sequence Charge Decoration (SCD) of protein sequence. SCD
    describes the patterning of charges along a protein sequence, weighting each
    pair of charges by the square root of their separation [3]. It is calculated as:

    SCD = (1/L) * sum_{m=2}^{L} sum_{n=1}^{m-1} q_m * q_n * (m - n)^(1/2)

    where q_m and q_n are the charges of the residues at positions m and n (+1 for
    K and R, -1 for D and E, 0 otherwise) and L is the length of the sequence. More
    negative values indicate more segregated opposite charges and uncharged
    sequences have an SCD of 0.

    Parameters
    ==========
    :sequence: str
        protein sequence.

    Returns
    =======
    :scd_df: pd.DataFrame
        pandas dataframe of SCD for protein sequence. Dataframe will be of the
        shape 1 x 1, with the column Charge_SCD.
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
    Calculate the Charge Density of protein sequence. Charge density is the net
    charge of a protein sequence at a given pH divided by its molecular weight [4],
    calculated as:

    ChargeDensity = NetCharge(pH) / MW

    where NetCharge(pH) is calculated using the charge_distribution function and MW
    using the molecular_weight function.

    Parameters
    ==========
    :sequence: str
        protein sequence.
    :ph: float (default=7.4)
        pH at which to calculate the net charge.

    Returns
    =======
    :charge_density_df: pd.DataFrame
        pandas dataframe of charge density (charges per Da) for protein sequence.
        Dataframe will be of the shape 1 x 1, with the column ChargeDensity. The
        value is rounded to 6 decimal places because it is typically of order 1e-4.
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
    Calculate the Hydrophobicity Profile (HP) of protein sequence. Each amino acid
    in the sequence is replaced by its value on a hydrophobicity scale and the
    values are averaged over a sliding window, giving a profile along the sequence
    [5, 6]. The profile is summarised by its mean, maximum, minimum and standard
    deviation, and optionally the fraction of windows above a threshold, e.g. 1.6
    for transmembrane segments with the Kyte-Doolittle scale [5]. The window value
    is calculated as:

    HP(i) = sum_j w_j * H(s_i+j) / sum_j w_j

    where H is the hydrophobicity scale and w_j are the window weights. With
    edge=1.0 all weights are 1 (a plain mean); with edge < 1 the weights fall
    linearly from 1 at the centre of the window to edge at its ends, as in ExPASy
    ProtScale. Windows longer than the sequence use the whole sequence.

    Parameters
    ==========
    :sequence: str
        protein sequence.
    :scale: str (default="kyte_doolittle")
        single-value hydrophobicity scale: "kyte_doolittle", "hopp_woods",
        "eisenberg", "wimley_white", "engelman", "domcut" or any AAindex
        accession number.
    :window: int (default=9)
        length of the sliding window.
    :edge: float (default=1.0)
        relative weight of the ends of the window, between 0 and 1.
    :threshold: float | None (default=None)
        if set, also calculate the fraction of windows with a value above it.

    Returns
    =======
    :profile_df: pd.DataFrame
        pandas dataframe of hydrophobicity profile statistics for protein sequence.
        Dataframe will be of the shape 1 x 4 (1 x 5 with a threshold), with the
        columns HP_<scale>_w<window>_mean, _max, _min, _std and _frac_above.
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
    Calculate the Hydrophobic Ratio of protein sequence. The hydrophobic ratio is
    the percentage of hydrophobic residues in a protein sequence, as defined by
    the Antimicrobial Peptide Database (APD) [7], calculated as:

    HydrophobicRatio = N(LIVFMACW) / L * 100

    where N(LIVFMACW) is the number of L, I, V, F, M, A, C and W residues and L
    is the length of the sequence.

    Parameters
    ==========
    :sequence: str
        protein sequence.

    Returns
    =======
    :hydrophobic_ratio_df: pd.DataFrame
        pandas dataframe of hydrophobic ratio for protein sequence. Dataframe will
        be of the shape 1 x 1, with the column HydrophobicRatio.
    """
    #validate input protein sequence
    sequence = _validate_sequence(sequence)
    return _to_df([sum(aa in _apd_hydrophobic for aa in sequence) / len(sequence) * 100], ["HydrophobicRatio"])

#################################### Half Life ############################################

def half_life(sequence: str) -> pd.DataFrame:
    """
    Calculate the estimated Half Life of protein sequence. The half-life is a
    prediction of the time it takes for half of the amount of a protein in a cell
    to disappear after its synthesis. It is estimated from the N-terminal residue
    using the N-end rule [8, 9, 10], for 3 model systems: mammalian reticulocytes
    (in vitro), yeast and E. coli (in vivo), using the values tabulated by ExPASy
    ProtParam [6]. Values given as lower bounds in ProtParam (">20 hours" for yeast,
    ">10 hours" for E. coli) are reported as the bound. The E. coli value for an
    N-terminal proline is unknown in ProtParam and returned as NaN. The estimate is
    not applicable to N-terminally modified proteins.

    Parameters
    ==========
    :sequence: str
        protein sequence.

    Returns
    =======
    :half_life_df: pd.DataFrame
        pandas dataframe of half-lives (hours) for protein sequence. Dataframe will
        be of the shape 1 x 3, with the columns HalfLife_Mammal_h, HalfLife_Yeast_h
        and HalfLife_Ecoli_h.
    """
    #validate input protein sequence
    sequence = _validate_sequence(sequence)

    row = _load_json_table("residue_tables.json", "n_end_half_life")["values"][sequence[0]]
    values = [np.nan if row[k] is None else row[k] for k in ("mammal", "yeast", "ecoli")]
    return _to_df(values, ["HalfLife_Mammal_h", "HalfLife_Yeast_h", "HalfLife_Ecoli_h"])

################################## IVYWREL Index ##########################################

def ivywrel_index(sequence: str) -> pd.DataFrame:
    """
    Calculate the IVYWREL Index of protein sequence. The IVYWREL index is the
    fraction of I, V, Y, W, R, E and L residues in a protein sequence, which
    correlates with the optimal growth temperature of the organism [11]. It is
    calculated as:

    IVYWREL = N(IVYWREL) / L

    where N(IVYWREL) is the number of I, V, Y, W, R, E and L residues and L is the
    length of the sequence.

    Parameters
    ==========
    :sequence: str
        protein sequence.

    Returns
    =======
    :ivywrel_df: pd.DataFrame
        pandas dataframe of IVYWREL index for protein sequence. Dataframe will be
        of the shape 1 x 1, with the column IVYWREL.
    """
    #validate input protein sequence
    sequence = _validate_sequence(sequence)
    return _to_df([sum(aa in "IVYWREL" for aa in sequence) / len(sequence)], ["IVYWREL"])

#################################### CvP Bias #############################################

def cvp_bias(sequence: str) -> pd.DataFrame:
    """
    Calculate the Charged versus Polar Bias (CvP) of protein sequence. The CvP bias
    is the difference between the percentages of charged and polar residues in a
    protein sequence, an indicator of hyperthermostability [12], calculated as:

    CvP = %(D + E + K + R) - %(N + Q + S + T)

    where %(X) is the percentage of residues of type X in the sequence.

    Parameters
    ==========
    :sequence: str
        protein sequence.

    Returns
    =======
    :cvp_df: pd.DataFrame
        pandas dataframe of CvP bias for protein sequence. Dataframe will be of the
        shape 1 x 1, with the column CvP_Bias.
    """
    #validate input protein sequence
    sequence = _validate_sequence(sequence)
    charged = sum(aa in "DEKR" for aa in sequence)
    polar = sum(aa in "NQST" for aa in sequence)
    return _to_df([(charged - polar) / len(sequence) * 100], ["CvP_Bias"])

################################### Flexibility ###########################################

def flexibility(sequence: str, window: int = 9) -> pd.DataFrame:
    """
    Calculate the Flexibility (Flex) profile of protein sequence. Each amino acid
    in the sequence is replaced by its normalised flexibility parameter (B-value)
    and the values are averaged over a sliding window of 9 residues whose end
    weights fall linearly to 0.25 (weights 0.25, 0.4375, 0.625, 0.8125, 1, ...,
    0.25), as in Vihinen et al. [13]. The profile is summarised by its mean,
    maximum, minimum and standard deviation. The average B-values of Vihinen et
    al. are used; their parameter sets for residues with 0, 1 or 2 rigid
    neighbours are not. Compared with Biopython's ProteinAnalysis.flexibility,
    which uses the same values and weights, each window is centred on its middle
    residue and the last window is included.

    Parameters
    ==========
    :sequence: str
        protein sequence.
    :window: int (default=9)
        length of the sliding window (the parameters were optimised for 9).

    Returns
    =======
    :flexibility_df: pd.DataFrame
        pandas dataframe of flexibility profile statistics for protein sequence.
        Dataframe will be of the shape 1 x 4, with the columns Flex_mean, Flex_max,
        Flex_min and Flex_std.
    """
    #validate input protein sequence
    sequence = _validate_sequence(sequence)
    window = _validate_int("window", window, 1)

    profile = _sliding_windows(_property_matrix(sequence, "vihinen_flexibility")[:, 0], window, edge=0.25)
    stats = _profile_stats(profile, "Flex")
    return _to_df(list(stats.values()), list(stats))
