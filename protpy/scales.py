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
    Calculate an amino acid factor-scale descriptor of protein sequence. Factor
    scales describe each amino acid with a small number of values (dimensions)
    derived from many physicochemical properties, typically by principal component
    or factor analysis [1-9]. Each amino acid in the sequence is replaced by its
    d-dimensional scale vector and the mean of each dimension is taken over the
    sequence, calculated as:

    S_j = (1/L) * sum_i z_j(s_i)

    where S_j is the descriptor value for dimension j (j = 1,2,..,d), z_j(s_i) is
    the value of dimension j for the amino acid at position i and L is the length
    of the sequence. With domains > 1 the sequence is split into that many
    near-equal contiguous segments and the means of each segment are concatenated.
    The individual scales are also available as their own functions, e.g. vhse().

    Parameters
    ==========
    :sequence: str
        protein sequence.
    :scale: str
        factor-scale family to use: "z_scales_5", "vhse", "atchley_factors",
        "t_scales", "st_scales", "protfp", "fasgai", "blosum_indices" or "ms_whim".
    :domains: int (default=1)
        number of near-equal segments of the sequence to calculate the means over,
        at most the length of the sequence.

    Returns
    =======
    :scale_df: pd.DataFrame
        pandas dataframe of factor-scale values for protein sequence. Dataframe
        will be of the shape 1 x d * domains, where d is the number of dimensions
        of the scale, d * domains is the number of features calculated from the
        descriptor and 1 is the input sequence. Column names are of the form
        ABBR_j, or ABBR_j_dk when domains > 1, e.g. VHSE_1 or VHSE_1_d2.
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


def z_scales_5(sequence: str, domains: int = 1) -> pd.DataFrame:
    """
    Calculate Sandberg 5 z-scales (Z5) of protein sequence. The z-scales are the
    principal components of a set of physicochemical properties measured for 87
    coded and non-coded amino acids [1]. z1 mainly describes lipophilicity, z2
    steric bulk and polarizability, z3 polarity and charge, and z4 and z5
    electronegativity, heat of formation, electrophilicity and hardness.

    Each amino acid in the sequence is replaced by its 5 z-scales and the mean of
    each is taken over the sequence, calculated as:

    S_j = (1/L) * sum_i z_j(s_i)

    where S_j is the descriptor value for scale j, z_j(s_i) is the value of scale
    j for the amino acid at position i and L is the length of the sequence. With
    domains > 1 the sequence is split into that many near-equal segments and the
    means of each segment are concatenated.

    Parameters
    ==========
    :sequence: str
        protein sequence.
    :domains: int (default=1)
        number of near-equal segments of the sequence to calculate the means over.

    Returns
    =======
    :z_scales_5_df: pd.DataFrame
        pandas dataframe of Z5 values for protein sequence. Dataframe will be
        of the shape 1 x 5 * domains, where 5 * domains is the number of
        features calculated from the descriptor (5 scales for each domain) and
        1 is the input sequence. Column names are of the form Z5_j, or
        Z5_j_dk when domains > 1.
    """
    return aa_scale_descriptors(sequence, "z_scales_5", domains)

def vhse(sequence: str, domains: int = 1) -> pd.DataFrame:
    """
    Calculate Vectors of Hydrophobic, Steric and Electronic properties (VHSE) of
    protein sequence. The VHSE scales are principal components score vectors derived
    from physicochemical variables of the 20 coded amino acids [2]. VHSE1-2 describe
    hydrophobic properties, VHSE3-4 steric properties and VHSE5-8 electronic
    properties.

    Each amino acid in the sequence is replaced by its 8 VHSE scales and the mean of
    each is taken over the sequence, calculated as:

    S_j = (1/L) * sum_i z_j(s_i)

    where S_j is the descriptor value for scale j, z_j(s_i) is the value of scale
    j for the amino acid at position i and L is the length of the sequence. With
    domains > 1 the sequence is split into that many near-equal segments and the
    means of each segment are concatenated.

    Parameters
    ==========
    :sequence: str
        protein sequence.
    :domains: int (default=1)
        number of near-equal segments of the sequence to calculate the means over.

    Returns
    =======
    :vhse_df: pd.DataFrame
        pandas dataframe of VHSE values for protein sequence. Dataframe will be
        of the shape 1 x 8 * domains, where 8 * domains is the number of
        features calculated from the descriptor (8 scales for each domain) and
        1 is the input sequence. Column names are of the form VHSE_j, or
        VHSE_j_dk when domains > 1.
    """
    return aa_scale_descriptors(sequence, "vhse", domains)

def atchley_factors(sequence: str, domains: int = 1) -> pd.DataFrame:
    """
    Calculate Atchley Factors (AF) of protein sequence. The Atchley factors are 5
    factors obtained from a factor analysis of 54 amino acid attributes [3]. Factor
    I mainly reflects polarity, accessibility and hydrophobicity, factor II
    propensity for secondary structure, factor III molecular size, factor IV
    relative amino acid composition and factor V electrostatic charge.

    Each amino acid in the sequence is replaced by its 5 Atchley factors and the
    mean of each is taken over the sequence, calculated as:

    S_j = (1/L) * sum_i z_j(s_i)

    where S_j is the descriptor value for scale j, z_j(s_i) is the value of scale
    j for the amino acid at position i and L is the length of the sequence. With
    domains > 1 the sequence is split into that many near-equal segments and the
    means of each segment are concatenated.

    Parameters
    ==========
    :sequence: str
        protein sequence.
    :domains: int (default=1)
        number of near-equal segments of the sequence to calculate the means over.

    Returns
    =======
    :atchley_factors_df: pd.DataFrame
        pandas dataframe of AF values for protein sequence. Dataframe will be
        of the shape 1 x 5 * domains, where 5 * domains is the number of
        features calculated from the descriptor (5 scales for each domain) and
        1 is the input sequence. Column names are of the form AF_j, or
        AF_j_dk when domains > 1.
    """
    return aa_scale_descriptors(sequence, "atchley_factors", domains)

def t_scales(sequence: str, domains: int = 1) -> pd.DataFrame:
    """
    Calculate Topological Scales (T-scales) of protein sequence. The T-scales are
    principal components of topological descriptors of amino acids, describing their
    structural and topological characteristics [4].

    Each amino acid in the sequence is replaced by its 5 T-scales and the mean of
    each is taken over the sequence, calculated as:

    S_j = (1/L) * sum_i z_j(s_i)

    where S_j is the descriptor value for scale j, z_j(s_i) is the value of scale
    j for the amino acid at position i and L is the length of the sequence. With
    domains > 1 the sequence is split into that many near-equal segments and the
    means of each segment are concatenated.

    Parameters
    ==========
    :sequence: str
        protein sequence.
    :domains: int (default=1)
        number of near-equal segments of the sequence to calculate the means over.

    Returns
    =======
    :t_scales_df: pd.DataFrame
        pandas dataframe of T values for protein sequence. Dataframe will be
        of the shape 1 x 5 * domains, where 5 * domains is the number of
        features calculated from the descriptor (5 scales for each domain) and
        1 is the input sequence. Column names are of the form T_j, or
        T_j_dk when domains > 1.
    """
    return aa_scale_descriptors(sequence, "t_scales", domains)

def st_scales(sequence: str, domains: int = 1) -> pd.DataFrame:
    """
    Calculate Structural Topology Scales (ST-scales) of protein sequence. The ST-
    scales are principal components of 3D structural and topological properties of
    amino acids [5].

    Each amino acid in the sequence is replaced by its 8 ST-scales and the mean of
    each is taken over the sequence, calculated as:

    S_j = (1/L) * sum_i z_j(s_i)

    where S_j is the descriptor value for scale j, z_j(s_i) is the value of scale
    j for the amino acid at position i and L is the length of the sequence. With
    domains > 1 the sequence is split into that many near-equal segments and the
    means of each segment are concatenated.

    Parameters
    ==========
    :sequence: str
        protein sequence.
    :domains: int (default=1)
        number of near-equal segments of the sequence to calculate the means over.

    Returns
    =======
    :st_scales_df: pd.DataFrame
        pandas dataframe of ST values for protein sequence. Dataframe will be
        of the shape 1 x 8 * domains, where 8 * domains is the number of
        features calculated from the descriptor (8 scales for each domain) and
        1 is the input sequence. Column names are of the form ST_j, or
        ST_j_dk when domains > 1.
    """
    return aa_scale_descriptors(sequence, "st_scales", domains)

def protfp(sequence: str, domains: int = 1) -> pd.DataFrame:
    """
    Calculate ProtFP (ProtFP) of protein sequence. ProtFP is a set of 8 principal
    components calculated from a selection of physicochemical properties in the
    AAindex database [6].

    Each amino acid in the sequence is replaced by its 8 ProtFP components and the
    mean of each is taken over the sequence, calculated as:

    S_j = (1/L) * sum_i z_j(s_i)

    where S_j is the descriptor value for scale j, z_j(s_i) is the value of scale
    j for the amino acid at position i and L is the length of the sequence. With
    domains > 1 the sequence is split into that many near-equal segments and the
    means of each segment are concatenated.

    Parameters
    ==========
    :sequence: str
        protein sequence.
    :domains: int (default=1)
        number of near-equal segments of the sequence to calculate the means over.

    Returns
    =======
    :protfp_df: pd.DataFrame
        pandas dataframe of ProtFP values for protein sequence. Dataframe will be
        of the shape 1 x 8 * domains, where 8 * domains is the number of
        features calculated from the descriptor (8 scales for each domain) and
        1 is the input sequence. Column names are of the form ProtFP_j, or
        ProtFP_j_dk when domains > 1.
    """
    return aa_scale_descriptors(sequence, "protfp", domains)

def fasgai(sequence: str, domains: int = 1) -> pd.DataFrame:
    """
    Calculate Factor Analysis Scales of Generalized Amino Acid Information (FASGAI)
    of protein sequence. The FASGAI scales are 6 factors obtained from a factor
    analysis of amino acid properties, describing hydrophobicity, alpha and turn
    propensities, bulky properties, compositional characteristics, local flexibility
    and electronic properties [7].

    Each amino acid in the sequence is replaced by its 6 FASGAI factors and the mean
    of each is taken over the sequence, calculated as:

    S_j = (1/L) * sum_i z_j(s_i)

    where S_j is the descriptor value for scale j, z_j(s_i) is the value of scale
    j for the amino acid at position i and L is the length of the sequence. With
    domains > 1 the sequence is split into that many near-equal segments and the
    means of each segment are concatenated.

    Parameters
    ==========
    :sequence: str
        protein sequence.
    :domains: int (default=1)
        number of near-equal segments of the sequence to calculate the means over.

    Returns
    =======
    :fasgai_df: pd.DataFrame
        pandas dataframe of FASGAI values for protein sequence. Dataframe will be
        of the shape 1 x 6 * domains, where 6 * domains is the number of
        features calculated from the descriptor (6 scales for each domain) and
        1 is the input sequence. Column names are of the form FASGAI_j, or
        FASGAI_j_dk when domains > 1.
    """
    return aa_scale_descriptors(sequence, "fasgai", domains)

def blosum_indices(sequence: str, domains: int = 1) -> pd.DataFrame:
    """
    Calculate BLOSUM Indices (BLOSUM) of protein sequence. The BLOSUM indices are 10
    factors derived from the BLOSUM62 substitution matrix by principal component
    analysis and VARIMAX rotation, so they describe the evolutionary similarity of
    amino acids [8].

    Each amino acid in the sequence is replaced by its 10 BLOSUM indices and the
    mean of each is taken over the sequence, calculated as:

    S_j = (1/L) * sum_i z_j(s_i)

    where S_j is the descriptor value for scale j, z_j(s_i) is the value of scale
    j for the amino acid at position i and L is the length of the sequence. With
    domains > 1 the sequence is split into that many near-equal segments and the
    means of each segment are concatenated.

    Parameters
    ==========
    :sequence: str
        protein sequence.
    :domains: int (default=1)
        number of near-equal segments of the sequence to calculate the means over.

    Returns
    =======
    :blosum_indices_df: pd.DataFrame
        pandas dataframe of BLOSUM values for protein sequence. Dataframe will be
        of the shape 1 x 10 * domains, where 10 * domains is the number of
        features calculated from the descriptor (10 scales for each domain) and
        1 is the input sequence. Column names are of the form BLOSUM_j, or
        BLOSUM_j_dk when domains > 1.
    """
    return aa_scale_descriptors(sequence, "blosum_indices", domains)

def ms_whim(sequence: str, domains: int = 1) -> pd.DataFrame:
    """
    Calculate Molecular Surface Weighted Holistic Invariant Molecular Scores (MS-
    WHIM) of protein sequence. The MS-WHIM scores are 3 principal components of 3D
    descriptors calculated from the molecular electrostatic potential on the surface
    of each amino acid [9].

    Each amino acid in the sequence is replaced by its 3 MS-WHIM scores and the mean
    of each is taken over the sequence, calculated as:

    S_j = (1/L) * sum_i z_j(s_i)

    where S_j is the descriptor value for scale j, z_j(s_i) is the value of scale
    j for the amino acid at position i and L is the length of the sequence. With
    domains > 1 the sequence is split into that many near-equal segments and the
    means of each segment are concatenated.

    Parameters
    ==========
    :sequence: str
        protein sequence.
    :domains: int (default=1)
        number of near-equal segments of the sequence to calculate the means over.

    Returns
    =======
    :ms_whim_df: pd.DataFrame
        pandas dataframe of MSWHIM values for protein sequence. Dataframe will be
        of the shape 1 x 3 * domains, where 3 * domains is the number of
        features calculated from the descriptor (3 scales for each domain) and
        1 is the input sequence. Column names are of the form MSWHIM_j, or
        MSWHIM_j_dk when domains > 1.
    """
    return aa_scale_descriptors(sequence, "ms_whim", domains)

############################## Auto/Cross Covariance ######################################

def auto_cross_covariance(sequence: str, scale: str | list[str] = "z_scales_5", lag: int = 5,
    mode: str = "acc") -> pd.DataFrame:
    """
    Calculate Auto Cross Covariance (ACC) of protein sequence. ACC transforms the
    amino acid scale values of a sequence into a fixed-length vector that
    describes the covariance between the scale dimensions of residues separated
    by l positions, for each lag l in the range 1..lag [10, 11]. It is calculated as:

    ACC_jk(l) = sum_{i=1}^{L-l} z_j(s_i) * z_k(s_i+l) / (L - l)

    where ACC_jk(l) is the covariance between dimensions j and k at lag l, z_j(s_i)
    is the value of dimension j for the amino acid at position i and L is the length
    of the sequence. When j = k this is the auto covariance and when j != k the cross
    covariance. The scale can be one of the factor-scale families, which are used as
    published, or one or more AAindex accession numbers, whose values are first
    standardised over the 20 amino acids.

    Parameters
    ==========
    :sequence: str
        protein sequence.
    :scale: str | list[str] (default="z_scales_5")
        factor-scale family, single-value scale name, AAindex accession number or
        list of AAindex accession numbers.
    :lag: int (default=5)
        maximum lag, must be less than the length of the sequence.
    :mode: str (default="acc")
        which covariances to return: "acc" (all dimension pairs), "ac" (auto
        covariance only, j = k) or "cc" (cross covariance only, j != k).

    Returns
    =======
    :acc_df: pd.DataFrame
        pandas dataframe of ACC values for protein sequence. Dataframe will be of
        the shape 1 x N, where N is the number of features calculated from the
        descriptor and 1 is the input sequence; N is d^2 * lag for "acc", d * lag
        for "ac" and d(d-1) * lag for "cc", where d is the number of dimensions of
        the scale. By default (5 z-scales, lag=5) the shape will be 1 x 125. Column
        names are of the form ACC_j_k_lagl, e.g. ACC_1_2_lag3.
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
    Calculate Physicochemical Distance Transformation (PDT) of protein sequence.
    PDT describes how much a physicochemical property changes between residues
    separated by l positions, for each lag l in the range 1..lag [12]. For each
    AAindex property it is calculated as:

    PDT_j(l) = sum_{i=1}^{L-l} (P_j(s_i) - P_j(s_i+l))^2 / (L - l)

    where PDT_j(l) is the descriptor value for property j at lag l, P_j(s_i) is the
    value of property j for the amino acid at position i, standardised over the 20
    amino acids, and L is the length of the sequence. By default the same 8 AAindex
    properties as the autocorrelation descriptors are used.

    Parameters
    ==========
    :sequence: str
        protein sequence.
    :lag: int (default=3)
        maximum lag, must be less than the length of the sequence.
    :properties: list[str] | None (default=None)
        list of AAindex accession numbers for the properties to use. Defaults to
        the 8 properties used by the autocorrelation descriptors.

    Returns
    =======
    :pdt_df: pd.DataFrame
        pandas dataframe of PDT values for protein sequence. Dataframe will be of
        the shape 1 x N, where N is the number of features calculated from the
        descriptor (number of properties * lag) and 1 is the input sequence. By
        default the shape will be 1 x 24. Column names are of the form
        PDT_<accession>_l, e.g. PDT_CIDH920105_1.
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
