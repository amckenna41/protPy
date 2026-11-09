################################################################################
#############                  Physicochemical                    ##############
################################################################################

from __future__ import annotations

"""
Charge and charge-patterning descriptors (charge fractions and NCPR, kappa,
SCD), the hydrophobicity profile engine and thermostability indices.

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
[14] von Heijne, G. (1986). The distribution of positively charged residues in
    bacterial inner membrane proteins correlates with the trans-membrane
    topology. EMBO J, 5(11), 3021-3027.
    https://doi.org/10.1002/j.1460-2075.1986.tb04601.x
"""

from functools import lru_cache

import numpy as np
import pandas as pd

from ._constants import _validate_sequence, _validate_int
from ._helpers import (_charge_vector, _to_df, _property_matrix, _sliding_windows, _profile_stats,
                       _load_json_table, _runs, _warn_deprecated)
from .composition import charge_distribution

__all__ = [
    "charge_fractions", "kappa", "sequence_charge_decoration", "hydrophobicity_profile",
    "thermostability_indices",
    #deprecated, removed in protpy 2.0
    "ncpr", "charge_density", "hydrophobic_ratio", "half_life", "ivywrel_index", "cvp_bias", "flexibility",
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

    NCPR = (n+ - n-) / L

    Sigma = NCPR^2 / FCR

    where n+ is the number of positive residues (K, R), n- the number of negative
    residues (D, E), L is the length of the sequence and NCPR is the net charge per
    residue [1]. His is treated as neutral and sigma is 0 when there are no charged
    residues. (Only two of the values are independent; they are all given because
    each is a standard localCIDER parameter.)

    Parameters
    ==========
    :sequence: str
        protein sequence.

    Returns
    =======
    :charge_fractions_df: pd.DataFrame
        pandas dataframe of charge fractions for protein sequence. Dataframe will
        be of the shape 1 x 5, where 5 is the number of features calculated from
        the descriptor and 1 is the input sequence. Column names are Charge_FCR,
        Charge_FPlus, Charge_FMinus, Charge_NCPR and Charge_Sigma.
    """
    #validate input protein sequence
    sequence = _validate_sequence(sequence)

    q = _charge_vector(sequence)
    f_plus, f_minus = np.mean(q > 0), np.mean(q < 0)
    fcr, ncpr_ = f_plus + f_minus, f_plus - f_minus
    sigma = ncpr_ ** 2 / fcr if fcr else 0.0
    #prefixed so "FCR" cannot clash with the tripeptide-composition column of the same name
    return _to_df([fcr, f_plus, f_minus, ncpr_, sigma],
                  ["Charge_FCR", "Charge_FPlus", "Charge_FMinus", "Charge_NCPR", "Charge_Sigma"])

#################################### NCPR #################################################

def ncpr(sequence: str) -> pd.DataFrame:
    """
    Calculate the Net Charge Per Residue (NCPR) of protein sequence. NCPR is the
    net charge of a protein sequence divided by its length [1], calculated as:

    NCPR = (n+ - n-) / L

    where n+ is the number of positive residues (K, R), n- the number of negative
    residues (D, E) and L is the length of the sequence.

    Deprecated since v1.6.0, removed in 2.0: NCPR is the Charge_NCPR column of
    ``charge_fractions``.

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
    _warn_deprecated("ncpr", "charge_fractions (column Charge_NCPR)")
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


@lru_cache(maxsize=4096)
def _delta_max(n_pos: int, n_neg: int, length: int, g: int) -> float:
    """
    Calculate the maximum delta of a sequence composition, used to normalise kappa
    [2]. The maximum is taken over the most segregated arrangements of the
    composition: a block of positive residues and a block of negative residues,
    with the neutral residues split across the three gaps around them. Blobs lying
    entirely in a neutral flank contribute the same value whatever the length of
    the flank, so for each gap between the blocks only flank lengths up to g - 1
    need to be tried, giving O(n0 * g) candidate arrangements, where n0 is the
    number of neutral residues.

    Each arrangement is scored in O(g) rather than O(L): only blobs that straddle
    one of the four block boundaries have a composition other than all-neutral
    (charge asymmetry 0) or all-charged (charge asymmetry 1), so delta is the sum
    over the at most 4(g - 1) straddling blobs plus the counts of pure blobs times
    their constant values. The result is memoised, as it depends only on the
    composition and the blob size.

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

    Returns
    =======
    :best: float
        maximum delta over the segregated arrangements.
    """
    sigma = (n_pos - n_neg) ** 2 / (length * (n_pos + n_neg))
    n0 = length - n_pos - n_neg
    candidates = set()
    for b in range(n0 + 1):
        rest = n0 - b
        for k in range(min(g, rest + 1)):
            candidates.add((k, b))
            candidates.add((rest - k, b))
    a, b = np.array(sorted(candidates)).T
    c = n0 - a - b
    n_blobs = length - g + 1
    #blob starts that straddle a block boundary
    boundaries = np.stack([a, a + n_pos, a + n_pos + b, a + n_pos + b + n_neg], axis=1)
    starts = (boundaries[:, :, None] + np.arange(1 - g, 0)).reshape(len(a), -1)
    starts = np.sort(np.where((starts >= 0) & (starts < n_blobs), starts, -1), axis=1)
    mixed = starts >= 0
    mixed[:, 1:] &= starts[:, 1:] != starts[:, :-1]
    ends = starts + g
    n_p = np.clip(np.minimum(ends, (a + n_pos)[:, None]) - np.maximum(starts, a[:, None]), 0, None)
    lo = (a + n_pos + b)[:, None]
    n_m = np.clip(np.minimum(ends, lo + n_neg) - np.maximum(starts, lo), 0, None)
    total = n_p + n_m
    blob_sigma = np.divide((n_p - n_m) ** 2, total * g, out=np.zeros(total.shape), where=total > 0)
    mixed_sum = np.where(mixed, (blob_sigma - sigma) ** 2, 0.0).sum(axis=1)
    #blobs lying inside a single block
    pure_neutral = sum(np.maximum(0, r - g + 1) for r in (a, b, c))
    pure_charged = max(0, n_pos - g + 1) + max(0, n_neg - g + 1)
    delta = (mixed_sum + pure_neutral * sigma ** 2 + pure_charged * (1 - sigma) ** 2) / n_blobs
    return float(delta.max())


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
    most segregated arrangements of the same composition. Kappa is the mean of
    kappa_g over the blob sizes (5 and 6 by default) and ranges from 0 (well
    mixed) to 1 (fully segregated). Sequences with only one type of charge are
    scored the same way (as in localCIDER); kappa is NaN for a sequence with no
    charged residues, where it is undefined. Blobs longer than the sequence are
    reduced to its length.

    Parameters
    ==========
    :sequence: str
        protein sequence.
    :blob_sizes: tuple (default=(5, 6))
        blob sizes g to average kappa over.

    Returns
    =======
    :kappa_df: pd.DataFrame
        pandas dataframe of kappa for protein sequence. Dataframe will be of the
        shape 1 x 1, with the column Kappa. (Before v1.6.0 the values for each blob
        size were also returned; they correlate at r >= 0.986 with Kappa. Call
        kappa(sequence, blob_sizes=(g,)) for a single blob size.)
    """
    #validate input protein sequence
    sequence = _validate_sequence(sequence)
    blob_sizes = [_validate_int("blob size", g, 1) for g in blob_sizes]
    if not blob_sizes:
        raise ValueError("At least one blob size must be given.")

    q = _charge_vector(sequence)
    pos, neg = (q > 0).astype(float), (q < 0).astype(float)
    n_pos, n_neg, L = int(pos.sum()), int(neg.sum()), len(sequence)
    if n_pos + n_neg == 0:
        return _to_df([np.nan], ["Kappa"])
    values = []
    for g in blob_sizes:
        g = min(g, L)
        sigma = (n_pos - n_neg) ** 2 / (L * (n_pos + n_neg))
        d_max = _delta_max(n_pos, n_neg, L, g)
        values.append(min(float(_delta(pos, neg, g, sigma)) / d_max, 1.0) if d_max else 0.0)
    return _to_df([float(np.mean(values))], ["Kappa"])

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

    Deprecated since v1.6.0, removed in 2.0: use
    ``charge_distribution(sequence, ph, normalize="mass")`` (column
    NetCharge_perDa); charge density correlates at r = 0.996 with NCPR.

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
        Dataframe will be of the shape 1 x 1, with the column ChargeDensity.
    """
    _warn_deprecated("charge_density", 'charge_distribution(sequence, ph, normalize="mass")')
    net = charge_distribution(sequence, ph=ph, normalize="mass")["NetCharge_perDa"].iloc[0]
    return _to_df([net], ["ChargeDensity"])

############################### Hydrophobicity Profile ####################################

def hydrophobicity_profile(sequence: str, scale: str = "kyte_doolittle", window: int = 9, edge: float = 1.0,
    threshold: float | None = None, segments: bool = False) -> pd.DataFrame:
    """
    Calculate the Hydrophobicity Profile (HP) of protein sequence. Each amino acid
    in the sequence is replaced by its value on a hydrophobicity scale (or any
    other single-value scale) and the values are averaged over a sliding window,
    giving a profile along the sequence [5, 6]. The profile is summarised by its
    mean, maximum, minimum and standard deviation, and optionally the fraction of
    windows above a threshold, e.g. 1.6 for transmembrane segments with the
    Kyte-Doolittle scale [5]. The window value is calculated as:

    HP(i) = sum_j w_j * H(s_i+j) / sum_j w_j

    where H is the scale and w_j are the window weights. With edge=1.0 all
    weights are 1 (a plain mean); with edge < 1 the weights fall linearly from 1 at
    the centre of the window to edge at its ends, as in ExPASy ProtScale. Windows
    longer than the sequence use the whole sequence.

    With ``segments=True`` the windows above ``threshold`` are merged into
    segments (the residues covered by at least one such window), and the number of
    segments, the fraction of residues inside a segment and the positive-inside
    bias (K + R in the 15 residues before the first segment minus K + R in the 15
    residues after it, 0 without a segment) [14] are added. With the Kyte-Doolittle
    scale, window=19 and threshold=1.6 these are hydropathy-based transmembrane
    segments; this is a heuristic, not a topology predictor (signal peptides are
    also counted, and multi-pass proteins are often under-counted).

    The scale can be any protpy scale (e.g. "vihinen_flexibility" with window=9
    and edge=0.25 gives the Vihinen flexibility profile [13], "parker" or "janin"
    with window=7 the Parker and Janin profiles) or an AAindex1 accession number.

    Parameters
    ==========
    :sequence: str
        protein sequence.
    :scale: str (default="kyte_doolittle")
        single-value scale: "kyte_doolittle", "hopp_woods", "eisenberg_jmb1984",
        "eisenberg_arb1984", "wimley_white", "engelman", "domcut",
        "vihinen_flexibility", "parker", "janin", "top_idp", "eiip" or any AAindex
        accession number.
    :window: int (default=9)
        length of the sliding window.
    :edge: float (default=1.0)
        relative weight of the ends of the window, between 0 and 1.
    :threshold: float | None (default=None)
        if set, also calculate the fraction of windows with a value above it.
    :segments: bool (default=False)
        also calculate the segment features; needs a threshold.

    Returns
    =======
    :profile_df: pd.DataFrame
        pandas dataframe of hydrophobicity profile statistics for protein sequence.
        Dataframe will be of the shape 1 x 4 (1 x 5 with a threshold, 1 x 8 with
        segments), with the columns HP_<scale>_w<window>_mean, _max, _min, _std,
        _frac_above, _n_segments, _frac_in_segments and _positive_inside_bias.
    """
    #validate input protein sequence
    sequence = _validate_sequence(sequence)
    window = _validate_int("window", window, 1)
    if not 0 <= edge <= 1:
        raise ValueError(f"edge must be between 0 and 1, got {edge}.")
    if segments and threshold is None:
        raise ValueError("segments=True needs a threshold, e.g. threshold=1.6 for Kyte-Doolittle with window=19.")

    values = _property_matrix(sequence, scale)
    if values.shape[1] != 1:
        raise ValueError(f"Scale '{scale}' is multi-dimensional; use a single-value scale.")
    profile = _sliding_windows(values[:, 0], window, edge)
    prefix = f"HP_{scale}_w{window}"
    stats = _profile_stats(profile, prefix, threshold=threshold)
    if segments:
        width = min(window, len(sequence))
        covered = np.zeros(len(sequence), dtype=bool)
        for start in np.flatnonzero(profile > threshold):
            covered[start:start + width] = True
        runs = _runs(covered)
        bias = 0
        if runs:
            s, e = runs[0]
            bias = sum(aa in "KR" for aa in sequence[max(0, s - 15):s]) - sum(aa in "KR" for aa in sequence[e:e + 15])
        stats.update({f"{prefix}_n_segments": len(runs), f"{prefix}_frac_in_segments": float(covered.mean()),
                      f"{prefix}_positive_inside_bias": bias})
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

    Deprecated since v1.6.0, removed in 2.0: it is a fixed linear map of the amino
    acid composition (and close to the hydrophobic group of
    ``kmer_composition(k=1, alphabet=3)``).

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
    _warn_deprecated("hydrophobic_ratio", "kmer_composition(sequence, k=1, alphabet=3)")
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

    Deprecated since v1.6.0, removed in 2.0: it depends only on the first residue
    (which is Met in most proteins, so it is nearly constant) and the three
    columns are perfectly correlated.

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
    _warn_deprecated("half_life")
    #validate input protein sequence
    sequence = _validate_sequence(sequence)

    row = _load_json_table("residue_tables.json", "n_end_half_life")["values"][sequence[0]]
    values = [np.nan if row[k] is None else row[k] for k in ("mammal", "yeast", "ecoli")]
    return _to_df(values, ["HalfLife_Mammal_h", "HalfLife_Yeast_h", "HalfLife_Ecoli_h"])

############################## Thermostability Indices ####################################

def thermostability_indices(sequence: str) -> pd.DataFrame:
    """
    Calculate the Thermostability Indices of protein sequence: the IVYWREL index
    and the charged versus polar (CvP) bias, composition-based indicators of
    protein thermostability [11, 12]. They are calculated as:

    IVYWREL = N(IVYWREL) / L

    CvP = %(D + E + K + R) - %(N + Q + S + T)

    where N(IVYWREL) is the number of I, V, Y, W, R, E and L residues, %(X) is the
    percentage of residues of type X and L is the length of the sequence. IVYWREL
    was derived as a proteome-wide predictor of optimal growth temperature and
    does not reliably rank orthologous proteins; CvP bias was validated on
    hyperthermophile proteomes. Both are fixed linear maps of the amino acid
    composition.

    Parameters
    ==========
    :sequence: str
        protein sequence.

    Returns
    =======
    :thermostability_df: pd.DataFrame
        pandas dataframe of the thermostability indices for protein sequence.
        Dataframe will be of the shape 1 x 2, with the columns IVYWREL and CvP_Bias.
    """
    #validate input protein sequence
    sequence = _validate_sequence(sequence)
    L = len(sequence)
    charged = sum(aa in "DEKR" for aa in sequence)
    polar = sum(aa in "NQST" for aa in sequence)
    return _to_df([sum(aa in "IVYWREL" for aa in sequence) / L, (charged - polar) / L * 100], ["IVYWREL", "CvP_Bias"])

def ivywrel_index(sequence: str) -> pd.DataFrame:
    """
    Calculate the IVYWREL Index of protein sequence, the fraction of I, V, Y, W, R,
    E and L residues [11].

    Deprecated since v1.6.0, removed in 2.0: use ``thermostability_indices``
    (column IVYWREL).

    Parameters
    ==========
    :sequence: str
        protein sequence.

    Returns
    =======
    :ivywrel_df: pd.DataFrame
        pandas dataframe of shape 1 x 1, with the column IVYWREL.
    """
    _warn_deprecated("ivywrel_index", "thermostability_indices (column IVYWREL)")
    return thermostability_indices(sequence)[["IVYWREL"]]

def cvp_bias(sequence: str) -> pd.DataFrame:
    """
    Calculate the Charged versus Polar Bias (CvP) of protein sequence,
    %(D + E + K + R) - %(N + Q + S + T) [12].

    Deprecated since v1.6.0, removed in 2.0: use ``thermostability_indices``
    (column CvP_Bias).

    Parameters
    ==========
    :sequence: str
        protein sequence.

    Returns
    =======
    :cvp_df: pd.DataFrame
        pandas dataframe of shape 1 x 1, with the column CvP_Bias.
    """
    _warn_deprecated("cvp_bias", "thermostability_indices (column CvP_Bias)")
    return thermostability_indices(sequence)[["CvP_Bias"]]

################################### Flexibility ###########################################

def flexibility(sequence: str, window: int = 9) -> pd.DataFrame:
    """
    Calculate the Flexibility (Flex) profile of protein sequence. Each amino acid
    in the sequence is replaced by its normalised flexibility parameter (B-value)
    and the values are averaged over a sliding window of 9 residues whose end
    weights fall linearly to 0.25, as in Vihinen et al. [13], using the average
    B-values (AAindex VINM940101, as in Biopython). The profile is summarised by
    its mean, maximum, minimum and standard deviation.

    Deprecated since v1.6.0, removed in 2.0: use
    ``hydrophobicity_profile(sequence, scale="vihinen_flexibility", window=9, edge=0.25)``,
    which gives the same values.

    Parameters
    ==========
    :sequence: str
        protein sequence.
    :window: int (default=9)
        length of the sliding window (the parameters were optimised for 9).

    Returns
    =======
    :flexibility_df: pd.DataFrame
        pandas dataframe of shape 1 x 4, with the columns Flex_mean, Flex_max,
        Flex_min and Flex_std.
    """
    _warn_deprecated("flexibility", 'hydrophobicity_profile(sequence, scale="vihinen_flexibility", window=9, edge=0.25)')
    df = hydrophobicity_profile(sequence, scale="vihinen_flexibility", window=window, edge=0.25)
    df.columns = ["Flex_mean", "Flex_max", "Flex_min", "Flex_std"]
    return df
