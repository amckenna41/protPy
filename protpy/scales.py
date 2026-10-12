################################################################################
#############                  Amino Acid Scales                  ##############
################################################################################

from __future__ import annotations

"""
Amino-acid factor scales (z-scales, VHSE, Atchley factors, T-/ST-scales,
ProtFP, FASGAI, BLOSUM indices, MS-WHIM) and the auto/cross-covariance and
physicochemical distance transformations built on them.

References
==========
[1] Sandberg, M. et al. (1998). New chemical descriptors relevant for the design
    of biologically active peptides. A multivariate characterization of 87 amino
    acids. J Med Chem, 41(14), 2481-2491. https://doi.org/10.1021/jm9700575
[2] Mei, H. et al. (2005). A new set of amino acid descriptors and its application
    in peptide QSARs. Biopolymers, 80(6), 775-786. https://doi.org/10.1002/bip.20296
[3] Atchley, W. R. et al. (2005). Solving the protein sequence metric problem.
    PNAS, 102(18), 6395-6400. https://doi.org/10.1073/pnas.0408677102
[4] Tian, F., Zhou, P., & Li, Z. (2007). T-scale as a novel vector of topological
    descriptors for amino acids and its application in QSARs of peptides. J Mol
    Struct, 830, 106-115. https://doi.org/10.1016/j.molstruc.2006.07.004
[5] Yang, L. et al. (2010). ST-scale as a novel amino acid descriptor and its
    application in QSAM of peptides and analogues. Amino Acids, 38(3), 805-816.
    https://doi.org/10.1007/s00726-009-0287-y
[6] van Westen, G. J. P. et al. (2013). Benchmarking of protein descriptor sets in
    proteochemometric modeling (part 1): comparative study of 13 amino acid
    descriptor sets. J Cheminform, 5, 41. https://doi.org/10.1186/1758-2946-5-41
[7] Liang, G., & Li, Z. (2007). Factor analysis scale of generalized amino acid
    information as the source of a new set of descriptors for elucidating the
    structure and activity relationships of cationic antimicrobial peptides. QSAR
    Comb Sci, 26(6), 754-763. https://doi.org/10.1002/qsar.200630145
[8] Georgiev, A. G. (2009). Interpretable numerical descriptors of amino acid
    space. J Comput Biol, 16(5), 703-723. https://doi.org/10.1089/cmb.2008.0173
[9] Zaliani, A., & Gancia, E. (1999). MS-WHIM scores for amino acids: a new
    3D-description for peptide QSAR and QSPR studies. J Chem Inf Comput Sci,
    39(3), 525-533. https://doi.org/10.1021/ci980211b
[10] Wold, S. et al. (1993). DNA and peptide sequences and chemical processes
    multivariately modelled by principal component analysis and partial
    least-squares projections to latent structures. Anal Chim Acta, 277(2),
    239-253. https://doi.org/10.1016/0003-2670(93)80437-P
[11] Lapinsh, M. et al. (2002). Development of proteo-chemometrics: a novel
    technology for the analysis of drug-receptor interactions. Mol Pharmacol,
    61(6), 1465. https://doi.org/10.1124/mol.61.6.1465
[12] Liu, B. et al. (2012). Using amino acid physicochemical distance
    transformation for fast protein remote homology detection. PLoS ONE, 7(9),
    e46633. https://doi.org/10.1371/journal.pone.0046633

The values of [1, 2, 4, 5, 7, 8, 9] were transcribed from the CC-BY descriptor
component published with [6]; see protpy/data/aa_scales.json for provenance.
"""

import itertools

import numpy as np
import pandas as pd

from ._constants import _validate_sequence, _validate_int
from ._helpers import _property_matrix, _read_json, _to_df

__all__ = [
    "aa_scale_descriptors", "z_scales_5", "vhse", "atchley_factors", "t_scales", "st_scales",
    "protfp", "fasgai", "blosum_indices", "ms_whim", "auto_cross_covariance",
    "physicochemical_distance_transformation",
]

#factor-scale families: name -> column abbreviation
factor_scales_ = {
    "z_scales_5": "Z5",
    "vhse": "VHSE",
    "atchley_factors": "AF",
    "t_scales": "T",
    "st_scales": "ST",
    "protfp": "ProtFP",
    "fasgai": "FASGAI",
    "blosum_indices": "BLOSUM",
    "ms_whim": "MSWHIM",
}

#default AAindex properties for PDT (same 8 as the autocorrelation descriptors)
_default_pdt_properties = ["CIDH920105", "BHAR880101", "CHAM820101", "CHAM820102",
    "CHOC760101", "BIGC670101", "CHAM810101", "DAYM780201"]

################################## Factor Scale Descriptors ###############################

def aa_scale_descriptors(sequence: str, scale: str, domains: int = 1) -> pd.DataFrame:
    """
    Calculate the mean of an amino-acid factor scale over a protein sequence. Each
    residue is replaced by its d-dimensional scale vector and each dimension is
    averaged:

    S_j = (1/L) sum_i z_j(s_i),  j = 1..d

    With ``domains > 1`` the sequence is split into ``domains`` near-equal
    contiguous segments and the per-segment means are concatenated.

    Parameters
    ==========
    :sequence: str
        protein sequence.
    :scale: str
        factor-scale family: "z_scales_5", "vhse", "atchley_factors", "t_scales",
        "st_scales", "protfp", "fasgai", "blosum_indices" or "ms_whim".
    :domains: int (default=1)
        number of segments, at most the sequence length.

    Returns
    =======
    :scale_df: pd.DataFrame
        pandas DataFrame of shape 1 x d·domains, columns ``{ABBR}_{j}`` or
        ``{ABBR}_{j}_d{k}`` when ``domains > 1``.
    """
    #validate input protein sequence
    sequence = _validate_sequence(sequence)
    if scale not in factor_scales_:
        raise ValueError(f"Invalid scale '{scale}'. Valid scales are: {list(factor_scales_)}.")
    domains = _validate_int("domains", domains, 1)
    if domains > len(sequence):
        raise ValueError(f"domains ({domains}) cannot exceed the sequence length ({len(sequence)}).")

    abbr = factor_scales_[scale]
    matrix = _property_matrix(sequence, scale)
    values, columns = [], []
    for k, segment in enumerate(np.array_split(matrix, domains), start=1):
        values += list(segment.mean(axis=0))
        suffix = f"_d{k}" if domains > 1 else ""
        columns += [f"{abbr}_{j}{suffix}" for j in range(1, matrix.shape[1] + 1)]
    return _to_df(values, columns)


def _scale_wrapper(scale: str, ref: str, dim: int):
    def wrapper(sequence: str, domains: int = 1) -> pd.DataFrame:
        return aa_scale_descriptors(sequence, scale, domains)
    wrapper.__name__ = scale
    wrapper.__qualname__ = scale
    wrapper.__doc__ = f"""
    Calculate the mean {scale} factor-scale descriptor of a protein sequence {ref}.
    See ``aa_scale_descriptors`` for the formula.

    Parameters
    ==========
    :sequence: str
        protein sequence.
    :domains: int (default=1)
        number of near-equal segments whose means are concatenated.

    Returns
    =======
    :scale_df: pd.DataFrame
        pandas DataFrame of shape 1 x {dim}·domains, columns ``{factor_scales_[scale]}_{{j}}``.
    """
    return wrapper

z_scales_5 = _scale_wrapper("z_scales_5", "(Sandberg 5 z-scales) [1]", 5)
vhse = _scale_wrapper("vhse", "(principal components score Vectors of Hydrophobic, Steric and Electronic properties) [2]", 8)
atchley_factors = _scale_wrapper("atchley_factors", "(Atchley factors I-V) [3]", 5)
t_scales = _scale_wrapper("t_scales", "(topological T-scales) [4]", 5)
st_scales = _scale_wrapper("st_scales", "(structural topology ST-scales) [5]", 8)
protfp = _scale_wrapper("protfp", "(ProtFP PCA, 8 components) [6]", 8)
fasgai = _scale_wrapper("fasgai", "(factor analysis scales of generalized amino acid information) [7]", 6)
blosum_indices = _scale_wrapper("blosum_indices", "(BLOSUM62-derived indices) [8]", 10)
ms_whim = _scale_wrapper("ms_whim", "(MS-WHIM scores) [9]", 3)

############################## Auto/Cross Covariance ######################################

def auto_cross_covariance(sequence: str, scale: str | list[str] = "z_scales_5", lag: int = 5,
    mode: str = "acc") -> pd.DataFrame:
    """
    Calculate the auto- and cross-covariance (ACC) transform of an amino-acid
    scale [10, 11]. For lag l and scale dimensions j, k:

    ACC_{j,k}(l) = sum_{i=1}^{L-l} z_{j,i} x z_{k,i+l} / (L - l)

    ``mode="acc"`` returns every (j, k) pair, ``"ac"`` only j = k (auto
    covariance) and ``"cc"`` only j != k (cross covariance). ``scale`` may be a
    protpy factor-scale family (used as published) or an AAindex1 code / list
    of codes, which are standardised over the 20 amino acids.

    Parameters
    ==========
    :sequence: str
        protein sequence.
    :scale: str | list[str] (default="z_scales_5")
        factor-scale family, single-value scale name or AAindex1 code(s).
    :lag: int (default=5)
        maximum lag, less than the sequence length.
    :mode: str (default="acc")
        "acc", "ac" or "cc".

    Returns
    =======
    :acc_df: pd.DataFrame
        pandas DataFrame of shape 1 x d²·lag ("acc"), d·lag ("ac") or
        d(d-1)·lag ("cc"), columns ``ACC_{j}_{k}_lag{l}``.
    """
    #validate input protein sequence
    sequence = _validate_sequence(sequence)
    lag = _validate_int("lag", lag, 1)
    if lag >= len(sequence):
        raise ValueError(f"lag ({lag}) must be less than the sequence length ({len(sequence)}).")
    if mode not in ("acc", "ac", "cc"):
        raise ValueError(f"Invalid mode '{mode}'. Valid modes are: ['acc', 'ac', 'cc'].")

    #published scale tables are used as-is; raw AAindex properties are standardised
    standardize = not (isinstance(scale, str) and scale in _read_json("aa_scales.json"))
    z = _property_matrix(sequence, scale, standardize=standardize)
    d = z.shape[1]
    pairs = [(j, k) for j, k in itertools.product(range(d), repeat=2)
             if mode == "acc" or (mode == "ac") == (j == k)]
    values, columns = [], []
    for l in range(1, lag + 1):
        #cov[j, k] = mean over i of z[i, j] * z[i + l, k]
        #einsum instead of matmul: macOS Accelerate BLAS raises spurious float warnings
        cov = np.einsum("ij,ik->jk", z[:-l], z[l:]) / (len(sequence) - l)
        values += [cov[j, k] for j, k in pairs]
        columns += [f"ACC_{j + 1}_{k + 1}_lag{l}" for j, k in pairs]
    return _to_df(values, columns)

########################## Physicochemical Distance Transformation ########################

def physicochemical_distance_transformation(sequence: str, lag: int = 3,
    properties: list[str] | None = None) -> pd.DataFrame:
    """
    Calculate the Physicochemical Distance Transformation (PDT) [12]. For each
    AAindex property P_j (standardised over the 20 amino acids) and lag l:

    PDT_j(l) = sum_{i=1}^{L-l} (P_j(s_i) - P_j(s_{i+l}))^2 / (L - l)

    Parameters
    ==========
    :sequence: str
        protein sequence.
    :lag: int (default=3)
        maximum lag, less than the sequence length.
    :properties: list[str] | None (default=None)
        AAindex1 accession codes; defaults to the 8 properties used by the
        autocorrelation descriptors.

    Returns
    =======
    :pdt_df: pd.DataFrame
        pandas DataFrame of shape 1 x len(properties)·lag, columns ``PDT_{code}_{l}``.
    """
    #validate input protein sequence
    sequence = _validate_sequence(sequence)
    lag = _validate_int("lag", lag, 1)
    if lag >= len(sequence):
        raise ValueError(f"lag ({lag}) must be less than the sequence length ({len(sequence)}).")
    properties = list(properties) if properties is not None else _default_pdt_properties
    if not properties:
        raise ValueError("At least one property must be given.")

    p = _property_matrix(sequence, properties, standardize=True)
    values, columns = [], []
    for j, code in enumerate(properties):
        for l in range(1, lag + 1):
            values.append(np.mean((p[:-l, j] - p[l:, j]) ** 2))
            columns.append(f"PDT_{code}_{l}")
    return _to_df(values, columns)
