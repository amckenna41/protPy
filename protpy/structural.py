################################################################################
#############            Structural & Functional Heuristics       ##############
################################################################################

from __future__ import annotations

"""
Sequence-based structural and functional heuristics: aggregation hot spots
(AGGRESCAN), charge-hydropathy disorder (Uversky, FoldIndex), disorder
propensity, transmembrane and signal-peptide features, motif counts,
prion-like domains (PLAAC), antigenicity profiles and cysteine features.

References
==========
[1] Conchillo-Solé, O. et al. (2007). AGGRESCAN: a server for the prediction and
    evaluation of "hot spots" of aggregation in polypeptides. BMC Bioinformatics,
    8, 65. https://doi.org/10.1186/1471-2105-8-65
[2] Uversky, V. N., Gillespie, J. R., & Fink, A. L. (2000). Why are "natively
    unfolded" proteins unstructured under physiologic conditions? Proteins,
    41(3), 415-427. https://doi.org/10.1002/1097-0134(20001115)41:3<415::AID-PROT130>3.0.CO;2-7
[3] Prilusky, J. et al. (2005). FoldIndex: a simple tool to predict whether a
    given protein sequence is intrinsically unfolded. Bioinformatics, 21(16),
    3435-3438. https://doi.org/10.1093/bioinformatics/bti537
[4] Dunker, A. K. et al. (2001). Intrinsically disordered protein. J Mol Graph
    Model, 19(1), 26-59. https://doi.org/10.1016/S1093-3263(00)00138-8
[5] Campen, A. et al. (2008). TOP-IDP-scale: a new amino acid scale measuring
    propensity for intrinsic disorder. Protein Pept Lett, 15(9), 956-963.
    https://doi.org/10.2174/092986608785849164
[6] Kyte, J., & Doolittle, R. F. (1982). A simple method for displaying the
    hydropathic character of a protein. J Mol Biol, 157(1), 105-132.
    https://doi.org/10.1016/0022-2836(82)90515-0
[7] von Heijne, G. (1986). The distribution of positively charged residues in
    bacterial inner membrane proteins correlates with the trans-membrane
    topology. EMBO J, 5(11), 3021-3027.
    https://doi.org/10.1002/j.1460-2075.1986.tb04601.x
[8] von Heijne, G. (1983). Patterns of amino acids near signal-sequence cleavage
    sites. Eur J Biochem, 133(1), 17-21.
    https://doi.org/10.1111/j.1432-1033.1983.tb07424.x
[9] Sigrist, C. J. A. et al. (2013). New and continuing developments at PROSITE.
    Nucleic Acids Res, 41(D1), D344-D347. https://doi.org/10.1093/nar/gks1067
[10] Olsen, J. V., Ong, S.-E., & Mann, M. (2004). Trypsin cleaves exclusively
    C-terminal to arginine and lysine residues. Mol Cell Proteomics, 3(6),
    608-614. https://doi.org/10.1074/mcp.T400003-MCP200
[11] Lancaster, A. K., Nutter-Upham, A., Lindquist, S., & King, O. D. (2014).
    PLAAC: a web and command-line application to identify proteins with
    prion-like amino acid composition. Bioinformatics, 30(17), 2501-2502.
    https://doi.org/10.1093/bioinformatics/btu310
[12] Parker, J. M. R., Guo, D., & Hodges, R. S. (1986). New hydrophilicity scale
    derived from high-performance liquid chromatography peptide retention data.
    Biochemistry, 25(19), 5425-5432. https://doi.org/10.1021/bi00367a013
[13] Janin, J. (1979). Surface and inside volumes in globular proteins. Nature,
    277, 491-492. https://doi.org/10.1038/277491a0
[14] Fomenko, D. E., & Gladyshev, V. N. (2003). Identity and functions of CxxC-
    derived motifs. Biochemistry, 42(38), 11214-11225. https://doi.org/10.1021/bi034459s
[15] Das, R. K., & Pappu, R. V. (2013). Conformations of intrinsically disordered
    proteins are influenced by linear sequence distributions of oppositely
    charged residues. PNAS, 110(33), 13392-13397.
    https://doi.org/10.1073/pnas.1304749110
"""

import math
import re

import numpy as np
import pandas as pd

from ._constants import _validate_sequence, _validate_int
from ._helpers import (_charge_vector, _load_json_table, _property_matrix, _profile_stats, _runs,
                       _sliding_windows, _to_df)
from .composition import kyte_doolittle_

__all__ = [
    "aggrescan", "uversky_charge_hydropathy", "fold_index", "disorder_propensity",
    "transmembrane_features", "signal_peptide_features", "functional_site_counts", "plaac",
    "antigenicity_surface", "cysteine_features",
]

################################### AGGRESCAN #############################################

#np.trapz was renamed np.trapezoid in numpy 2.0
_trapz = getattr(np, "trapezoid", None) or np.trapz

#hot-spot threshold: SwissProt-frequency-weighted mean a3v [1]
_AGGRESCAN_HST = -0.02

def aggrescan(sequence: str, window: int = 5) -> pd.DataFrame:
    """
    Calculate AGGRESCAN (AGG) aggregation descriptors of protein sequence.
    AGGRESCAN predicts aggregation "hot spots", short stretches that promote
    protein aggregation, from experimental aggregation propensities (a3v) of the
    20 amino acids [1]. The a3v values are averaged over a sliding window, giving
    the aggregation profile (a4v), with each window value assigned to its central
    residue and the off-centre N- and C-terminal residues given the first and last
    window values. A hot spot is a run of 5 or more residues with an a4v above the
    hot-spot threshold (HST) of -0.02 and containing no proline. Areas are
    calculated by trapezoidal integration of the profile, with HST as the zero axis.

    The 12 output values are: the mean a3v of the sequence (a3vSA); the number of
    hot spots (nHS) and the number per 100 residues (NnHS); the area of the profile
    above HST (AAT) and per residue (AATr); the total area of the hot spots (THSA)
    and per residue (THSAr); the total area of the profile (TA); the sum of a4v
    (a4vSS) and per 100 residues (Na4vSS); the mean a4v in the hot spots (a4vAHS,
    0 if there are none); and the largest single hot-spot area (HSA_max).

    AGGRESCAN chooses a window of 5, 7, 9 or 11 residues from the sequence length
    and gives the two extreme residues specific values; neither rule is given in
    the published supplementary material, so the window is a parameter and the
    extreme residues are treated like the other off-centre residues.

    Parameters
    ==========
    :sequence: str
        protein sequence.
    :window: int (default=5)
        odd sliding window length (AGGRESCAN uses 5, 7, 9 or 11); reduced to the
        largest odd length that fits a short sequence.

    Returns
    =======
    :aggrescan_df: pd.DataFrame
        pandas dataframe of AGGRESCAN descriptors for protein sequence. Dataframe
        will be of the shape 1 x 12, where 12 is the number of features calculated
        from the descriptor and 1 is the input sequence. Column names start
        with AGG, e.g. AGG_nHS and AGG_THSA.
    """
    #validate input protein sequence
    sequence = _validate_sequence(sequence)
    window = _validate_int("window", window, 1)
    if window % 2 == 0:
        raise ValueError(f"window must be odd, got {window}.")
    L = len(sequence)
    window = min(window, L if L % 2 else L - 1)

    a3v = _property_matrix(sequence, _load_json_table("structural_params.json", "aggrescan_a3v")["values"])[:, 0]
    centres = _sliding_windows(a3v, window)
    half = window // 2
    a4v = np.concatenate([np.full(half, centres[0]), centres, np.full(half, centres[-1])])
    above = a4v - _AGGRESCAN_HST

    #hot spots: >= 5 consecutive residues above HST, broken by proline
    is_hs = (above > 0) & np.array([aa != "P" for aa in sequence])
    hot_spots = _runs(is_hs, min_len=5)
    hsa = [_trapz(above[s:e]) for s, e in hot_spots]
    aat = float(_trapz(np.clip(above, 0, None)))
    thsa = float(sum(hsa))
    values = [a3v.mean(), len(hot_spots), len(hot_spots) / L * 100, aat, aat / L, thsa, thsa / L,
              float(_trapz(above)), a4v.sum(), a4v.sum() / L * 100,
              float(np.mean([a4v[s:e].mean() for s, e in hot_spots])) if hot_spots else 0.0,
              max(hsa) if hsa else 0.0]
    columns = ["AGG_a3vSA", "AGG_nHS", "AGG_NnHS", "AGG_AAT", "AGG_AATr", "AGG_THSA", "AGG_THSAr",
               "AGG_TA", "AGG_a4vSS", "AGG_Na4vSS", "AGG_a4vAHS", "AGG_HSA_max"]
    return _to_df(values, columns)

######################### Uversky Charge-Hydropathy / FoldIndex ###########################

def _mean_hydropathy(sequence: str) -> float:
    """
    Calculate the mean hydropathy of a protein sequence using the Kyte-Doolittle
    scale rescaled to the range 0-1, (KD + 4.5) / 9, as used by Uversky et al. [2].

    Parameters
    ==========
    :sequence: str
        protein sequence.

    Returns
    =======
    :mean_h: float
        mean rescaled hydropathy of the sequence.
    """
    return sum((kyte_doolittle_[aa] + 4.5) / 9 for aa in sequence) / len(sequence)


def uversky_charge_hydropathy(sequence: str) -> pd.DataFrame:
    """
    Calculate the Uversky Charge-Hydropathy (UV) descriptors of protein sequence.
    Uversky et al. showed that natively unfolded proteins can be separated from
    folded proteins by their mean net charge and mean hydropathy [2]. These are
    calculated as:

    <H> = mean of (KD + 4.5) / 9

    <R> = abs(n+ - n-) / L

    HQ = 2.785 <H> - <R> - 1.151

    where KD is the Kyte-Doolittle hydropathy of each residue, n+ and n- are the
    numbers of positive (K, R) and negative (D, E) residues and L is the length
    of the sequence. HQ is the distance from the boundary <R> = 2.785 <H> - 1.151
    (as restated in [15]); the protein is predicted to be unfolded when HQ < 0.

    Parameters
    ==========
    :sequence: str
        protein sequence.

    Returns
    =======
    :uversky_df: pd.DataFrame
        pandas dataframe of charge-hydropathy descriptors for protein sequence.
        Dataframe will be of the shape 1 x 4, with the columns UV_MeanH, UV_MeanR,
        UV_HQ and UV_Unfolded (1 if predicted unfolded, else 0).
    """
    #validate input protein sequence
    sequence = _validate_sequence(sequence)
    h = _mean_hydropathy(sequence)
    r = abs(_charge_vector(sequence).mean())
    hq = 2.785 * h - r - 1.151
    return _to_df([h, r, hq, int(hq < 0)], ["UV_MeanH", "UV_MeanR", "UV_HQ", "UV_Unfolded"])


def fold_index(sequence: str, window: int = 51) -> pd.DataFrame:
    """
    Calculate the FoldIndex (FI) of protein sequence. FoldIndex applies the
    Uversky charge-hydropathy boundary over a sliding window to find disordered
    regions of a protein sequence [3]. For each window it is calculated as:

    FI = 2.785 <H> - abs(<R>) - 1.151

    where <H> is the mean Kyte-Doolittle hydropathy of the window rescaled to 0-1
    and <R> is the mean net charge of the window (K, R positive; D, E negative).
    Residues covered by at least one window with FI < 0 are counted as disordered.
    A window longer than the sequence uses the whole sequence, so FI then equals
    UV_HQ from uversky_charge_hydropathy.

    Parameters
    ==========
    :sequence: str
        protein sequence.
    :window: int (default=51)
        length of the sliding window.

    Returns
    =======
    :fold_index_df: pd.DataFrame
        pandas dataframe of FoldIndex descriptors for protein sequence. Dataframe
        will be of the shape 1 x 4, with the columns FI_mean, FI_min,
        FI_frac_disordered and FI_longest_disordered (in residues).
    """
    #validate input protein sequence
    sequence = _validate_sequence(sequence)
    window = min(_validate_int("window", window, 1), len(sequence))

    h = _sliding_windows([(kyte_doolittle_[aa] + 4.5) / 9 for aa in sequence], window)
    r = np.abs(_sliding_windows(_charge_vector(sequence), window))
    fi = 2.785 * h - r - 1.151
    disordered = np.zeros(len(sequence), dtype=bool)
    for start in np.flatnonzero(fi < 0):
        disordered[start:start + window] = True
    runs = _runs(disordered)
    longest = max((e - s for s, e in runs), default=0)
    return _to_df([fi.mean(), fi.min(), disordered.mean(), longest],
                  ["FI_mean", "FI_min", "FI_frac_disordered", "FI_longest_disordered"])

################################ Disorder Propensity ######################################

def disorder_propensity(sequence: str) -> pd.DataFrame:
    """
    Calculate the Disorder Propensity (DP) of protein sequence. Intrinsically
    disordered proteins are enriched in disorder-promoting residues (A, R, G, Q,
    S, P, E, K) and depleted in order-promoting residues (W, C, F, I, Y, V, L, N)
    [4]. The descriptor gives the fraction of each, together with the mean value of
    the TOP-IDP scale, a scale of amino acid propensity for intrinsic disorder [5],
    calculated as:

    DP_TOP_IDP_mean = sum_i TOP-IDP(s_i) / L

    where TOP-IDP(s_i) is the TOP-IDP value of the amino acid at position i and L
    is the length of the sequence. Higher values indicate a higher propensity for
    disorder.

    Parameters
    ==========
    :sequence: str
        protein sequence.

    Returns
    =======
    :disorder_df: pd.DataFrame
        pandas dataframe of disorder propensity descriptors for protein sequence.
        Dataframe will be of the shape 1 x 3, with the columns
        DP_frac_disorder_promoting, DP_frac_order_promoting and DP_TOP_IDP_mean.
    """
    #validate input protein sequence
    sequence = _validate_sequence(sequence)
    classes = _load_json_table("structural_params.json", "disorder_residue_classes")["values"]
    top_idp = _load_json_table("structural_params.json", "top_idp")["values"]
    L = len(sequence)
    return _to_df([sum(aa in classes["disorder_promoting"] for aa in sequence) / L,
                   sum(aa in classes["order_promoting"] for aa in sequence) / L,
                   sum(top_idp[aa] for aa in sequence) / L],
                  ["DP_frac_disorder_promoting", "DP_frac_order_promoting", "DP_TOP_IDP_mean"])

############################## Transmembrane Features #####################################

def transmembrane_features(sequence: str, window: int = 19, kd_threshold: float = 1.6) -> pd.DataFrame:
    """
    Calculate Transmembrane (TM) segment features of protein sequence. Putative
    transmembrane segments are found from the Kyte-Doolittle hydropathy profile:
    windows of 19 residues with a mean hydropathy above 1.6 indicate a likely
    transmembrane segment [6], and overlapping windows are merged into segments.
    The positive-inside rule states that the cytoplasmic side of a transmembrane
    segment is enriched in positive residues [7]. The 4 output values are: the
    number of segments; the maximum window hydropathy; the fraction of residues
    inside a segment; and the positive-inside bias, the number of K and R residues
    in the 15 residues before the first segment minus the number in the 15 residues
    after it (0 if there is no segment). This is a hydropathy heuristic and not a
    topology predictor.

    Parameters
    ==========
    :sequence: str
        protein sequence.
    :window: int (default=19)
        length of the hydropathy window.
    :kd_threshold: float (default=1.6)
        mean hydropathy above which a window is part of a segment.

    Returns
    =======
    :tm_df: pd.DataFrame
        pandas dataframe of transmembrane features for protein sequence. Dataframe
        will be of the shape 1 x 4, with the columns TM_n_segments, TM_max_KD,
        TM_frac_in_TM and TM_positive_inside_bias.
    """
    #validate input protein sequence
    sequence = _validate_sequence(sequence)
    window = min(_validate_int("window", window, 1), len(sequence))

    profile = _sliding_windows([kyte_doolittle_[aa] for aa in sequence], window)
    in_tm = np.zeros(len(sequence), dtype=bool)
    for start in np.flatnonzero(profile > kd_threshold):
        in_tm[start:start + window] = True
    segments = _runs(in_tm)
    bias = 0
    if segments:
        s, e = segments[0]
        bias = sum(aa in "KR" for aa in sequence[max(0, s - 15):s]) - sum(aa in "KR" for aa in sequence[e:e + 15])
    return _to_df([len(segments), profile.max(), in_tm.mean(), bias],
                  ["TM_n_segments", "TM_max_KD", "TM_frac_in_TM", "TM_positive_inside_bias"])

############################## Signal Peptide Features ####################################

def signal_peptide_features(sequence: str, n_term: int = 30) -> pd.DataFrame:
    """
    Calculate Signal Peptide (SP) features of protein sequence. Signal peptides
    have a positively charged N-terminal n-region followed by a hydrophobic
    h-region [8]. The 2 output values are: the net charge of the first 5 residues
    (K and R +1, D and E -1); and the maximum mean Kyte-Doolittle hydropathy over
    8-residue windows within the first n_term residues. The von Heijne
    cleavage-site weight matrix is not included.

    Parameters
    ==========
    :sequence: str
        protein sequence.
    :n_term: int (default=30)
        length of the N-terminal region to examine.

    Returns
    =======
    :sp_df: pd.DataFrame
        pandas dataframe of signal peptide features for protein sequence.
        Dataframe will be of the shape 1 x 2, with the columns SP_nregion_charge
        and SP_hregion_maxKD.
    """
    #validate input protein sequence
    sequence = _validate_sequence(sequence)
    n_term = _validate_int("n_term", n_term, 1)

    region = sequence[:n_term]
    charge = _charge_vector(sequence[:5]).sum()
    h = _sliding_windows([kyte_doolittle_[aa] for aa in region], 8).max()
    return _to_df([int(charge), h], ["SP_nregion_charge", "SP_hregion_maxKD"])

############################## Functional Site Counts #####################################

#PROSITE patterns [9] and the trypsin rule [10]; lookaheads count overlapping matches
functional_motifs_ = {
    "N_glycosylation": r"(?=N[^P][ST][^P])",   #PS00001 N-{P}-[ST]-{P}
    "ER_retention": r"[KRHQSA][DENQ]EL$",      #PS00014 [KRHQSA]-[DENQ]-E-L>
    "PTS1": r"[STAGCN][RKH][LIVMAFY]$",        #PS00342 [STAGCN]-[RKH]-[LIVMAFY]>
    "tryptic_site": r"[KR](?!P)",              #trypsin: after K/R, not before P
}

def functional_site_counts(sequence: str) -> pd.DataFrame:
    """
    Calculate Functional Site (FS) counts of protein sequence. The number of each
    functional sequence motif, and its density per 100 residues, is counted for:
    N-glycosylation sites (PROSITE PS00001, N-{P}-[ST]-{P}), the C-terminal
    endoplasmic reticulum retention signal (PS00014, [KRHQSA]-[DENQ]-E-L at the
    C-terminus), the C-terminal peroxisomal targeting signal 1 (PS00342,
    [STAGCN]-[RKH]-[LIVMAFY] at the C-terminus) [9] and tryptic cleavage sites
    (after K or R, not before P) [10]. Overlapping matches are counted. These motifs
    complement the default motifs of motif_composition.

    Parameters
    ==========
    :sequence: str
        protein sequence.

    Returns
    =======
    :sites_df: pd.DataFrame
        pandas dataframe of functional site counts for protein sequence. Dataframe
        will be of the shape 1 x 8, where 8 is the number of features calculated
        from the descriptor (4 counts and 4 densities) and 1 is the input sequence.
        Column names are FS_<motif> and FS_<motif>_per100, e.g. FS_N_glycosylation.
    """
    #validate input protein sequence
    sequence = _validate_sequence(sequence)
    counts = {name: len(re.findall(pattern, sequence)) for name, pattern in functional_motifs_.items()}
    values = list(counts.values()) + [c / len(sequence) * 100 for c in counts.values()]
    columns = [f"FS_{n}" for n in counts] + [f"FS_{n}_per100" for n in counts]
    return _to_df(values, columns)

###################################### PLAAC ##############################################

def _plaac_model() -> tuple[dict, dict, np.ndarray, np.ndarray]:
    """
    Get the parameters of the PLAAC prion-like domain model [11]: the amino acid
    frequencies of S. cerevisiae prion domains (foreground) and of the
    S. cerevisiae proteome (background), and the transition and initial
    probabilities of the two-state hidden Markov model. Ported from the PLAAC
    source code (MIT License).

    Returns
    =======
    :bg: dict
        background amino acid frequencies.
    :fg: dict
        prion-domain amino acid frequencies.
    :log_t: np.ndarray
        log transition probabilities of the hidden Markov model (2 x 2).
    :log_i: np.ndarray
        log initial probabilities of the hidden Markov model.
    """
    p = _load_json_table("structural_params.json", "plaac")["values"]
    eps = p["epsilon"]

    def _normalise(freqs: dict) -> dict:
        #as PLAAC: normalise, add epsilon for X and *, normalise again
        """
        Normalise amino acid frequencies as PLAAC does: normalise to sum to 1, then
        add a small probability for the unknown residue (X) and stop (*) symbols and
        normalise again.

        Parameters
        ==========
        :freqs: dict
            amino acid frequencies.

        Returns
        =======
        :freqs: dict
            normalised amino acid frequencies.
        """
        total = sum(freqs.values())
        freqs = {aa: v / total for aa, v in freqs.items()}
        total = 1 + 2 * eps
        return {aa: v / total for aa, v in freqs.items()}

    return (_normalise(p["bg"]), _normalise(p["fg"]), np.log(np.array(p["transition"])),
            np.log(np.array(p["initial"])))


def _viterbi(sequence: str, bg: dict, fg: dict, log_t: np.ndarray, log_i: np.ndarray) -> np.ndarray:
    """
    Calculate the most likely sequence of states of the PLAAC two-state hidden
    Markov model (background or prion-like) for a protein sequence using the
    Viterbi algorithm [11].

    Parameters
    ==========
    :sequence: str
        protein sequence.
    :bg: dict
        background amino acid frequencies.
    :fg: dict
        prion-domain amino acid frequencies.
    :log_t: np.ndarray
        log transition probabilities.
    :log_i: np.ndarray
        log initial probabilities.

    Returns
    =======
    :path: np.ndarray
        state of each residue, 0 for background and 1 for prion-like.
    """
    emit = np.log(np.array([[bg[aa], fg[aa]] for aa in sequence]))
    score = log_i + emit[0]
    back = np.zeros((len(sequence), 2), dtype=int)
    for i in range(1, len(sequence)):
        cand = score[:, None] + log_t
        back[i] = np.argmax(cand, axis=0)
        score = cand[back[i], [0, 1]] + emit[i]
    path = np.zeros(len(sequence), dtype=int)
    path[-1] = int(np.argmax(score))
    for i in range(len(sequence) - 1, 0, -1):
        path[i - 1] = back[i, path[i]]
    return path


def _max_window_sum(values: np.ndarray, width: int) -> float:
    """
    Calculate the maximum sum of values over any window of a given width.

    Parameters
    ==========
    :values: np.ndarray
        values along the sequence.
    :width: int
        window width.

    Returns
    =======
    :max_sum: float
        maximum window sum.
    """
    csum = np.concatenate([[0.0], np.cumsum(values)])
    return float((csum[width:] - csum[:-width]).max())


def plaac(sequence: str, core_length: int = 60) -> pd.DataFrame:
    """
    Calculate the Prion-Like Amino Acid Composition (PLAAC) descriptors of protein
    sequence. PLAAC identifies prion-like domains, whose composition resembles
    that of yeast prion domains, using the amino acid frequencies of prion domains
    and of the S. cerevisiae proteome [11]. Each residue has a log-likelihood ratio:

    LLR(s_i) = ln(f_prion(s_i) / f_background(s_i))

    The 4 output values are: the maximum sum of LLR over any window of core_length
    residues (PLAAC_LLR); the maximum sum of LLR over a window of core_length
    residues lying entirely inside a prion-like region of the two-state hidden
    Markov model Viterbi parse (PLAAC_COREscore); the length of the prion-like
    region containing that core (PLAAC_PrD_length); and whether the sequence has a
    prion-like region of at least core_length residues (PLAAC_has_PrD). As in
    PLAAC, PLAAC_LLR and PLAAC_COREscore are NaN when they are undefined (a
    sequence shorter than core_length, or no prion-like region of that length).
    The parameters are ported from the PLAAC source code (MIT License), using the
    S. cerevisiae background (alpha = 1).

    Parameters
    ==========
    :sequence: str
        protein sequence.
    :core_length: int (default=60)
        minimum length of a prion-like domain core.

    Returns
    =======
    :plaac_df: pd.DataFrame
        pandas dataframe of PLAAC descriptors for protein sequence. Dataframe will
        be of the shape 1 x 4, with the columns PLAAC_LLR, PLAAC_COREscore,
        PLAAC_PrD_length and PLAAC_has_PrD.
    """
    #validate input protein sequence
    sequence = _validate_sequence(sequence)
    core_length = _validate_int("core_length", core_length, 1)

    bg, fg, log_t, log_i = _plaac_model()
    llr = np.array([math.log(fg[aa] / bg[aa]) for aa in sequence])
    window_llr = _max_window_sum(llr, core_length) if len(sequence) >= core_length else np.nan

    core, prd_length = np.nan, 0
    for s, e in _runs(_viterbi(sequence, bg, fg, log_t, log_i) == 1, min_len=core_length):
        score = _max_window_sum(llr[s:e], core_length)
        if np.isnan(core) or score > core:
            core, prd_length = score, e - s
    return _to_df([window_llr, core, prd_length, int(prd_length > 0)],
                  ["PLAAC_LLR", "PLAAC_COREscore", "PLAAC_PrD_length", "PLAAC_has_PrD"])

############################### Antigenicity / Surface ####################################

def antigenicity_surface(sequence: str, window: int = 7) -> pd.DataFrame:
    """
    Calculate the Antigenicity and Surface (Antigenicity) profiles of protein
    sequence. Antigenic sites tend to be hydrophilic and exposed on the protein
    surface. Each amino acid in the sequence is replaced by its value on the
    Parker hydrophilicity scale, derived from HPLC peptide retention times [12],
    and on the Janin scale of free energy of transfer from the inside to the
    outside of a globular protein [13]. Each is averaged over a sliding window of 7
    residues, as used by Parker et al., and the profile is summarised by its mean,
    maximum, minimum and standard deviation.

    Parameters
    ==========
    :sequence: str
        protein sequence.
    :window: int (default=7)
        length of the sliding window.

    Returns
    =======
    :antigenicity_df: pd.DataFrame
        pandas dataframe of antigenicity and surface profile statistics for protein
        sequence. Dataframe will be of the shape 1 x 8, with the columns
        PARKER_mean, PARKER_max, PARKER_min, PARKER_std and the same for JANIN.
    """
    #validate input protein sequence
    sequence = _validate_sequence(sequence)
    window = _validate_int("window", window, 1)
    stats = {}
    for name in ("parker", "janin"):
        profile = _sliding_windows(_property_matrix(sequence, name)[:, 0], window)
        stats.update(_profile_stats(profile, name.upper()))
    return _to_df(list(stats.values()), list(stats))

################################ Cysteine Features ########################################

def cysteine_features(sequence: str) -> pd.DataFrame:
    """
    Calculate Cysteine (Cys) features of protein sequence. Cysteines form
    disulfide bonds and functional motifs such as the CxxC motif of thiol
    oxidoreductases [14]. The 6 output values are: the number of cysteines; their
    fraction; whether the number is even; the maximum number of disulfide bonds
    (the number of cysteines // 2); the number of overlapping C-x-x-C motifs; and
    the mean gap between consecutive cysteines (0 if there are fewer than 2).

    Parameters
    ==========
    :sequence: str
        protein sequence.

    Returns
    =======
    :cys_df: pd.DataFrame
        pandas dataframe of cysteine features for protein sequence. Dataframe will
        be of the shape 1 x 6, with the columns Cys_count, Cys_frac, Cys_even,
        Cys_max_disulfides, Cys_CxxC and Cys_mean_spacing.
    """
    #validate input protein sequence
    sequence = _validate_sequence(sequence)
    positions = [i for i, aa in enumerate(sequence) if aa == "C"]
    n = len(positions)
    spacing = float(np.mean(np.diff(positions))) if n > 1 else 0.0
    return _to_df([n, n / len(sequence), int(n % 2 == 0), n // 2,
                   len(re.findall(r"(?=C..C)", sequence)), spacing],
                  ["Cys_count", "Cys_frac", "Cys_even", "Cys_max_disulfides", "Cys_CxxC", "Cys_mean_spacing"])
