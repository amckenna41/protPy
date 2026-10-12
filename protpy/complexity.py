################################################################################
#############                Sequence Complexity                  ##############
################################################################################

from __future__ import annotations

"""
Sequence complexity descriptors: Wootton-Federhen (SEG-like) complexity,
linguistic complexity and Lempel-Ziv complexity. ``shannon_entropy`` lives in
the composition module and is re-exported here.

References
==========
[1] Wootton, J. C., & Federhen, S. (1993). Statistics of local complexity in amino
    acid sequences and sequence databases. Comput Chem, 17(2), 149-163.
    https://doi.org/10.1016/0097-8485(93)85006-X
[2] Wootton, J. C., & Federhen, S. (1996). Analysis of compositionally biased
    regions in sequence databases. Methods Enzymol, 266, 554-571.
    https://doi.org/10.1016/S0076-6879(96)66035-2
[3] Troyanskaya, O. G., Arbell, O., Koren, Y., Landau, G. M., & Bolshoy, A.
    (2002). Sequence complexity profiles of prokaryotic genomic sequences: a fast
    algorithm for calculating linguistic complexity. Bioinformatics, 18(5),
    679-688. https://doi.org/10.1093/bioinformatics/18.5.679
[4] Lempel, A., & Ziv, J. (1976). On the complexity of finite sequences. IEEE
    Trans Inf Theory, 22(1), 75-81. https://doi.org/10.1109/TIT.1976.1055501
[5] Diao, Y. et al. (2008). Using pseudo amino acid composition to predict
    transmembrane regions in protein: cellular automata and Lempel-Ziv
    complexity. Amino Acids, 34(1), 111-117. https://doi.org/10.1007/s00726-007-0550-z
"""

import math
from collections import Counter

import numpy as np
import pandas as pd

from ._constants import _validate_sequence, _validate_int
from ._helpers import _to_df
from .composition import shannon_entropy

__all__ = [
    "wootton_federhen_complexity", "linguistic_complexity", "lempel_ziv_complexity", "shannon_entropy",
]

########################### Wootton-Federhen Complexity ###################################

def wootton_federhen_complexity(sequence: str, window: int = 12, k2_trigger: float = 2.2) -> pd.DataFrame:
    """
    Calculate the Wootton-Federhen Complexity (WF) of protein sequence. The
    Wootton-Federhen measures describe the local compositional complexity of a
    protein sequence over a sliding window, and are the basis of the SEG program
    for finding low-complexity regions [1, 2]. For each window of W residues
    they are calculated as:

    K1 = (1/W) * log_20(W! / prod_i n_i!)

    K2 = -sum_i (n_i / W) * log2(n_i / W)

    where K1 is the complexity state, K2 is the Shannon entropy of the window in
    bits and n_i is the number of residues of amino acid i in the window. The
    low-complexity fraction is the fraction of residues covered by at least one
    window with K2 <= k2_trigger. This is the trigger stage of SEG only; SEG's
    extension and merging of low-complexity segments are not performed. A window
    longer than the sequence uses the whole sequence.

    Parameters
    ==========
    :sequence: str
        protein sequence.
    :window: int (default=12)
        length of the sliding window W (the SEG default).
    :k2_trigger: float (default=2.2)
        K2 value in bits at or below which a window is low complexity (the SEG
        default).

    Returns
    =======
    :wf_df: pd.DataFrame
        pandas dataframe of Wootton-Federhen complexity for protein sequence.
        Dataframe will be of the shape 1 x 4, with the columns WF_K1_mean,
        WF_K2_mean, WF_K2_min and WF_LowComplexity_frac.
    """
    #validate input protein sequence
    sequence = _validate_sequence(sequence)
    window = min(_validate_int("window", window, 1), len(sequence))

    k1, k2 = [], []
    covered = np.zeros(len(sequence), dtype=bool)
    counts = Counter(sequence[:window])
    for start in range(len(sequence) - window + 1):
        if start:
            #slide the window by one residue
            counts[sequence[start - 1]] -= 1
            counts[sequence[start + window - 1]] += 1
        n = [c for c in counts.values() if c]
        log_states = math.lgamma(window + 1) - sum(math.lgamma(c + 1) for c in n)
        k1.append(log_states / math.log(20) / window)
        entropy = 0.0 - sum(c / window * math.log2(c / window) for c in n)
        k2.append(entropy)
        if entropy <= k2_trigger:
            covered[start:start + window] = True
    return _to_df([np.mean(k1), np.mean(k2), min(k2), covered.mean()],
                  ["WF_K1_mean", "WF_K2_mean", "WF_K2_min", "WF_LowComplexity_frac"])

############################### Linguistic Complexity #####################################

def linguistic_complexity(sequence: str, max_k: int | None = None) -> pd.DataFrame:
    """
    Calculate the Linguistic Complexity (LC) of protein sequence. Linguistic
    complexity describes how repetitive a sequence is, as the number of distinct
    subwords (k-mers) it contains relative to the maximum possible number for a
    sequence of its length [3]. It is calculated as:

    LC = sum_{k=1}^{K} U_k / sum_{k=1}^{K} min(20^k, L - k + 1)

    where U_k is the number of distinct k-mers in the sequence, L is the length of
    the sequence and K is the longest subword length (by default L, i.e. all
    subwords). Repetitive sequences have a low LC and random sequences an LC
    close to 1.

    Parameters
    ==========
    :sequence: str
        protein sequence.
    :max_k: int | None (default=None)
        longest subword length K; None uses the length of the sequence.

    Returns
    =======
    :lc_df: pd.DataFrame
        pandas dataframe of LC for protein sequence. Dataframe will be of the
        shape 1 x 1, with the column LinguisticComplexity.
    """
    #validate input protein sequence
    sequence = _validate_sequence(sequence)
    L = len(sequence)
    max_k = L if max_k is None else min(_validate_int("max_k", max_k, 1), L)

    observed = possible = 0
    for k in range(1, max_k + 1):
        n_kmers = L - k + 1
        observed += len({sequence[i:i + k] for i in range(n_kmers)})
        #20^k overtakes L - k + 1 quickly; avoid huge integers
        possible += n_kmers if k > math.log(n_kmers, 20) else min(20 ** k, n_kmers)
    return _to_df([observed / possible], ["LinguisticComplexity"])

############################### Lempel-Ziv Complexity #####################################

def _lz76(sequence: str) -> int:
    """
    Calculate the number of components in the Lempel-Ziv (LZ76) parsing of a
    sequence [4], using the Kaspar-Schuster algorithm. The sequence is read from
    left to right and each component is the shortest substring that cannot be
    copied from the preceding part of the sequence.

    Parameters
    ==========
    :sequence: str
        protein sequence.

    Returns
    =======
    :c: int
        number of components c(n) in the LZ76 parsing.
    """
    n = len(sequence)
    if n == 1:
        return 1
    c, i, k, l, k_max = 1, 0, 1, 1, 1
    while True:
        if sequence[i + k - 1] == sequence[l + k - 1]:
            k += 1
            if l + k > n:
                c += 1
                break
        else:
            k_max = max(k, k_max)
            i += 1
            if i == l:
                #no earlier copy: the component ends here
                c += 1
                l += k_max
                if l + 1 > n:
                    break
                i, k, k_max = 0, 1, 1
            else:
                k = 1
    return c


def lempel_ziv_complexity(sequence: str) -> pd.DataFrame:
    """
    Calculate the Lempel-Ziv Complexity (LZ) of protein sequence. The Lempel-Ziv
    (LZ76) complexity is the number of new patterns, c(n), found when a sequence is
    read from left to right, each pattern being the shortest substring that cannot
    be copied from the preceding part of the sequence [4, 5]. As c(n) grows with the
    length of the sequence, a normalised value is also calculated as:

    LZ_Normalized = c(n) * log_20(n) / n

    where n is the length of the sequence and 20 is the size of the amino acid
    alphabet. The normalised value is 0 for a single residue.

    Parameters
    ==========
    :sequence: str
        protein sequence.

    Returns
    =======
    :lz_df: pd.DataFrame
        pandas dataframe of LZ complexity for protein sequence. Dataframe will be
        of the shape 1 x 2, with the columns LZ_Complexity and LZ_Normalized.
    """
    #validate input protein sequence
    sequence = _validate_sequence(sequence)
    n = len(sequence)
    c = _lz76(sequence)
    return _to_df([c, c * math.log(n, 20) / n if n > 1 else 0.0], ["LZ_Complexity", "LZ_Normalized"])
