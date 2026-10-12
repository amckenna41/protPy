################################################################################
############                    Conjoint Triad                     #############
################################################################################

from __future__ import annotations

from collections import Counter
import pandas as pd

"""
References
==========
[1] J. Shen et al., “Predicting protein-protein interactions based only on sequences
    information,” Proc. Natl. Acad. Sci. U. S. A., vol. 104, no. 11, pp. 4337–4341, 2007.
    https://doi.org/10.1073/pnas.0607879104
[2] Chen, Z. et al. (2018). iFeature: a Python package and web server for features
    extraction and selection from protein and peptide sequences. Bioinformatics,
    34(14), 2499-2502. https://doi.org/10.1093/bioinformatics/bty140
"""
#list of amino acids
from ._constants import amino_acids, _validate_sequence, _validate_int
    
#amino acid triads
aa_triads = {
    1: ["A", "G", "V"],
    2: ["I", "L", "F", "P"],
    3: ["Y", "M", "T", "S"],
    4: ["H", "N", "Q", "W"],
    5: ["R", "K"],
    6: ["D", "E"],
    7: ["C"],
    }

def conjoint_triad(sequence: str, gap: int = 0) -> pd.DataFrame:
    """
    Calculate Conjoint Triad features (CTriad) for a protein sequence. The descriptor
    mainly considers neighbor relationships in protein sequences by encoding
    each protein sequence using the triad (continuous three amino acids)
    frequency distribution extracted from a 7-letter reduced alphabet [1]. This
    descriptor calculates 343 different features (7x7x7), with the output
    being of shape 1 x 343 for a sequence.

    With ``gap > 0`` the k-spaced conjoint triad (KSCTriad) is returned: for each
    g = 0..gap, triads (s_i, s_{i+g+1}, s_{i+2g+2}) are counted, giving an output
    of shape 1 x 343·(gap+1) with columns suffixed ``_gap{g}`` [2]. Values are
    raw counts, as for the default descriptor.

    Parameters
    ==========
    :sequence: str
        protein sequence.
    :gap: int (default=0)
        maximum residue spacing for KSCTriad; 0 gives the original conjoint triad.

    Returns
    =======
    :conjoint_triad_df: pd.Dataframe
        pandas Dataframe of CTriad descriptor values for all protein sequences. Dataframe
        will be of the shape 1 x 343, where 343 is the number of features calculated 
        from the descriptor and 1 is the input sequence.
    """
    #validate input protein sequence
    sequence = _validate_sequence(sequence)
    gap = _validate_int("gap", gap)

    conjoint_triad = {}
    _aa_triads = {}

    #expand aa_triads into form {AminoAcid: triad}
    for i in aa_triads:
        for j in aa_triads[i]:
            _aa_triads[j] = i

    #get protein number for each triad
    protein_num = sequence
    for i in _aa_triads:
        protein_num = protein_num.replace(i, str(_aa_triads[i]))

    for g in range(gap + 1):
        #count overlapping (g-spaced) triads with a sliding window (str.count misses overlaps, e.g. "111" in "11111")
        step = g + 1
        triad_counts = Counter(protein_num[i] + protein_num[i+step] + protein_num[i+2*step]
            for i in range(len(protein_num) - 2*step))
        suffix = f"_gap{g}" if gap > 0 else ""

        #calculate 7 triads for amino acids
        for i in range(1, 8):
            for j in range(1, 8):
                for k in range(1, 8):
                    temp = str(i) + str(j) + str(k)
                    conjoint_triad["conj_triad_" + temp + suffix] = triad_counts[temp]
    
    #transform descriptor values to dataframe
    conjoint_triad_df = pd.DataFrame([list(conjoint_triad.values())], columns=list(conjoint_triad.keys()))

    return conjoint_triad_df