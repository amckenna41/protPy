################################################################################
#############                   Composition                       ##############
################################################################################

from __future__ import annotations

import re
import itertools
from collections import Counter
import numpy as np
import pandas as pd
import math
from aaindex import aaindex1

"""
References
==========
[1] Gromiha, M. M. (2010). Protein Sequence Analysis. In M. M. Gromiha (Ed.), 
    Protein Bioinformatics (pp. 29–62). Elsevier.
[2] Hua, S. and Sun, Z. (2001) Support vector machine approach for protein
    subcellular localization prediction. Bioinformatics, 17, 721-728.
[3] Grassmann, J., Reczko, M., Suhai, S. and Edler, L. (1999) Protein fold
    class prediction: new methods of statistical classification. Proc Int Conf
    Intell Syst Mol Biol, 106-112.
[4] Kuo-Chen Chou. Prediction of Protein Cellular Attributes Using Pseudo-Amino Acid
    Composition. PROTEINS: Structure, Function, and Genetics, 2001, 43: 246-255.
[5] Kuo-Chen Chou, Using amphiphilic pseudo amino acid composition to predict enzyme 
    subfamily classes, Bioinformatics, Volume 21, Issue 1, January 2005, Pages 10–19, 
    https://doi.org/10.1093/bioinformatics/bth466
[6] http://www.csbio.sjtu.edu.cn/bioinf/PseAAC/ParaValue.htm.
[7] Reczko, M. and Bohr, H. (1994) The DEF data base of sequence based protein
    fold class predictions. Nucleic Acids Res, 22, 3616-3619.
[8] Kyte, J., & Doolittle, R. F. (1982). A simple method for displaying the
    hydropathic character of a protein. Journal of Molecular Biology, 157(1),
    105-132. https://doi.org/10.1016/0022-2836(82)90515-0
[9] Lobry, J. R., & Gautier, C. (1994). Hydrophobicity, expressivity and
    aromaticity are the major trends of amino-acid usage in 999 Escherichia coli
    chromosome-encoded genes. Nucleic Acids Research, 22(15), 3174-3180.
[10] Guruprasad, K., Reddy, B. V. B., & Pandit, M. W. (1990). Correlation
    between stability of a protein and its dipeptide composition: a novel
    approach for predicting in vivo stability of a protein from its primary
    sequence. Protein Engineering, 4(2), 155-161.
[11] Bjellqvist, B., et al. (1993). The focusing positions of polypeptides in
    immobilized pH gradients can be predicted from their amino acid sequences.
    Electrophoresis, 14(1), 1023-1031.
[12] Gasteiger, E., et al. (2005). Protein Identification and Analysis Tools on
    the ExPASy Server. In J. M. Walker (Ed.), The Proteomics Protocols Handbook
    (pp. 571-607). Humana Press.
[13] Cameselle-Teijeiro, J. (1979). The Henderson-Hasselbalch equation. 
    Biochemical Education, 7(3), 69-70.
[14] Chou, P. Y., & Fasman, G. D. (1974). Conformational parameters for amino
    acids in helical, beta-sheet, and random coil regions calculated from
    proteins. Biochemistry, 13(2), 211-222.
[15] Prosite database of protein families and domains:
    https://prosite.expasy.org/
[16] Ikai, A. J. (1980). Thermostability and aliphatic index of globular
    proteins. Journal of Biochemistry, 88(6), 1895-1898.
[17] Gasteiger, E., et al. (2005). Protein Identification and Analysis Tools
    on the ExPASy Server. In J. M. Walker (Ed.), The Proteomics Protocols
    Handbook (pp. 571-607). Humana Press. (extinction coefficient method)
[18] Boman, H. G. (2003). Antibacterial peptides: basic facts and emerging
    concepts. Journal of Internal Medicine, 254(3), 197-215.
[19] Eisenberg, D., Weiss, R. M., & Terwilliger, T. C. (1982). The helical
    hydrophobic moment: a measure of the amphiphilicity of a helix. Nature,
    299, 371-374. https://doi.org/10.1038/299371a0
[20] Chen, K., Kurgan, L. A., & Ruan, J. (2007). Prediction of flexible/rigid
    regions from protein sequences using k-spaced amino acid pairs. BMC
    Structural Biology, 7, 25. https://doi.org/10.1186/1472-6807-7-25
[21] Chen, Z. et al. (2018). iFeature: a Python package and web server for
    features extraction and selection from protein and peptide sequences.
    Bioinformatics, 34(14), 2499-2502. https://doi.org/10.1093/bioinformatics/bty140
[22] Saravanan, V., & Gautham, N. (2015). Harnessing computational biology for
    exact linear B-cell epitope prediction: a novel amino acid composition-based
    feature descriptor. OMICS, 19(10), 648-658. https://doi.org/10.1089/omi.2015.0095
[23] Wei, L. et al. (2018). ACPred-FL: a sequence-based predictor using effective
    feature representation to improve the prediction of anti-cancer peptides.
    Bioinformatics, 34(23), 4007-4016. https://doi.org/10.1093/bioinformatics/bty451
[24] Liu, B. et al. (2014). iDNA-Prot|dis: identifying DNA-binding proteins by
    incorporating amino acid distance-pairs and reduced alphabet profile into
    the general pseudo amino acid composition. PLoS ONE, 9(9), e106691.
    https://doi.org/10.1371/journal.pone.0106691
[25] Cannata, N. et al. (2002). Simplifying amino acid alphabets by means of a
    branch and bound algorithm and substitution matrices. Bioinformatics, 18(8),
    1102-1108. https://doi.org/10.1093/bioinformatics/18.8.1102
[26] Yu, C. et al. (2013). Protein space: a natural method for realizing the
    nature of protein universe. Journal of Theoretical Biology, 318, 197-204.
    https://doi.org/10.1016/j.jtbi.2012.11.005
[27] Pande, A. et al. (2023). Computing wide range of protein/peptide features
    from their sequence and structure. Journal of Computational Biology, 30(2),
    204-222. https://doi.org/10.1089/cmb.2022.0241
[28] Kumar, M., Verma, R., & Raghava, G. P. S. (2006). Prediction of
    mitochondrial proteins using support vector machine and hidden Markov
    model. Journal of Biological Chemistry, 281(9), 5357-5363.
    https://doi.org/10.1074/jbc.M511061200
[29] Lata, S., Sharma, B. K., & Raghava, G. P. S. (2007). Analysis and
    prediction of antibacterial peptides. BMC Bioinformatics, 8, 263.
    https://doi.org/10.1186/1471-2105-8-263
"""
import warnings

#list of amino acids
from ._constants import amino_acids, _validate_sequence, _validate_int, _validate_weight
from ._helpers import (_load_json_table, _read_json, _gapped_pair_counts, _kmer_counts, _to_df, _property_matrix,
                       _warn_deprecated, _check_lag, _runs)
#shannon_entropy moved to the complexity module in v1.6.0; re-exported here for compatibility
from .complexity import shannon_entropy  # noqa: F401 (re-export)

#public descriptors of this module; the correlation helpers used by PseAAC/APAAC
#(sequence_order_correlation_factor, amphiphilic_sequence_order_correlation_factor,
#normalize_property) are internal and no longer exported from the protpy namespace
__all__ = [
    "amino_acid_composition", "dipeptide_composition", "tripeptide_composition", "gravy", "aromaticity",
    "instability_index", "isoelectric_point", "molecular_weight", "charge_distribution",
    "kmer_composition", "motif_composition", "aliphatic_index", "extinction_coefficient", "boman_index",
    "hydrophobic_moment", "pseudo_amino_acid_composition", "amphiphilic_pseudo_amino_acid_composition",
    "cksaap", "cksaagp", "dde", "asdc", "distance_pair", "natural_vector", "atomic_composition",
    #deprecated, removed in protpy 2.0
    "reduced_alphabet_composition", "hydrophobic_polar_charged_composition", "secondary_structure_propensity",
    "amino_acid_pair_composition", "aggregation_propensity", "correlation_function",
    "amphiphilic_correllation_function",
]

#physicochemical values for all amino acids for hydrophobicity, hydrophilicty and residue mass
#values taken from http://www.csbio.sjtu.edu.cn/bioinf/PseAAC/ParaValue.htm
#property values differ from their equivalent in the aaindex therefore hard-coding them below
#hydrophobicity: the normalised consensus scale of Eisenberg et al. (1984), as listed by ExPASy
#ProtScale; it is the "eisenberg_jmb1984" table in protpy/data/aa_scales.json
hydrophobicity_ = dict(_read_json("aa_scales.json")["eisenberg_jmb1984"]["values"])

hydrophilicity_ = {
        "A": -0.5, "C": -1.0, "D": 3.0, "E": 3.0, "F": -2.5,  "G": 0.0, 
        "H": -0.5, "I": -1.8, "K": 3.0, "L": -1.8, "M": -1.3, "N": 0.2, 
        "P": 0.0, "Q": 0.2, "R": 3.0, "S": 0.3, "T": -0.4, "V": -1.5, 
        "W": -3.4, "Y": -2.3
        }

residue_mass_ = {
        "A": 15.0, "C": 47.0, "D": 59.0, "E": 73.0, "F": 91.0, "G": 1.000, 
        "H": 82.0, "I": 57.0, "K": 73.0, "L": 57.0,  "M": 75.0, "N": 58.0, 
        "P": 42.0, "Q": 72.0, "R": 101.0, "S": 31.0, "T": 45.0, "V": 43.0, 
        "W": 130.0, "Y": 107.0
        }

################################# Amino Acid Composition ##################################

def amino_acid_composition(sequence: str, region: str | None = None, n: int = 25) -> pd.DataFrame:
    """
    Calculate Amino Acid Composition (AAComp) of protein sequence. AAComp
    describes the fraction of each amino acid type within a protein sequence,
    and is calculated as:

    AA_Comp(s) = AA(t)/N(s)

    where AA_Comp(s) is the AAComp of protein sequence s, AA(t) is the number
    of amino acid types t (where t = 1,2,..,20) and N(s) is the length of the
    sequences [1].

    Parameters
    ==========
    :sequence: str
        protein sequence.
    :region: str | None (default=None)
        restrict the calculation to part of the sequence: None (whole sequence),
        "nterm" (first ``n`` residues, columns prefixed ``NT_``), "cterm" (last
        ``n`` residues, ``CT_``) or "split" (N-terminal, middle and C-terminal
        blocks concatenated, prefixed ``NT_``, ``MID_``, ``CT_``; split
        composition [28, 29]).
    :n: int (default=25)
        terminal block length used by ``region``.

    Returns
    =======
    :amino_acid_composition_df: pd.DataFrame
        pandas dataframe of AAComp for protein sequence. Dataframe will
        be of the shape 1 x 20, where 20 is the number of features 
        calculated from the descriptor (for the 20 amino acids) and 1 is 
        the input sequence.
    """
    #validate input protein sequence
    sequence = _validate_sequence(sequence)
    if region is not None:
        return _regional(amino_acid_composition, sequence, region, n)

    aa_composition = {}

    #iterate through each amino acid, calculating each composition (percentage)
    for aa in amino_acids:
        aa_composition[aa] = sequence.count(aa) / len(sequence) * 100

    #transform values and columns to DataFrame
    amino_acid_composition_df = pd.DataFrame([list(aa_composition.values())], 
        columns=list(aa_composition.keys()))

    return amino_acid_composition_df

################################## Dipeptide Composition ##################################

def dipeptide_composition(sequence: str, region: str | None = None, n: int = 25,
    normalize: str | None = None) -> pd.DataFrame:
    """
    Calculate Dipeptide Composition (DPComp) for protein sequence.
    Dipeptide composition is the fraction of each dipeptide type within a
    protein sequence. With dipeptides being of length 2 and there being 20
    canonical amino acids this creates 20^2 different combinations, thus a
    400-Dimensional vector will be produced such that:

    DPComp(s,t) = AA(s,t) / (N - 1) * 100

    where DPComp(s,t) is the dipeptide composition of the protein sequence
    for amino acid type s and t (where s and t = 1,2,..,20), AA(s,t) is the number
    of dipeptides represented by amino acid type s and t and N - 1 is the total number
    of dipeptides [1]. Values are percentages.

    With ``normalize="dde"`` the Dipeptide Deviation from Expected Mean (DDE) is
    returned instead: each composition is compared with the composition expected
    from the number of codons that code for each amino acid [22]:

    Dc(XY) = N(XY) / (L - 1)

    Tm(XY) = (C_X / 61) * (C_Y / 61)

    Tv(XY) = Tm(XY) * (1 - Tm(XY)) / (L - 1)

    DDE(XY) = (Dc(XY) - Tm(XY)) / sqrt(Tv(XY))

    where C_X and C_Y are the number of codons for amino acids X and Y out of the
    61 sense codons. DDE is a per-dipeptide affine transform of DPComp scaled by
    sqrt(L - 1), so for sequences of equal length it carries the same information.
    The columns are then named DDE_XY. ``dde()`` is an alias of this option.

    Parameters
    ==========
    :sequence: str
        protein sequence.
    :region: str | None (default=None)
        restrict the calculation to part of the sequence: None (whole sequence),
        "nterm" (first ``n`` residues, columns prefixed ``NT_``), "cterm" (last
        ``n`` residues, ``CT_``) or "split" (N-terminal, middle and C-terminal
        blocks concatenated, prefixed ``NT_``, ``MID_``, ``CT_``; split
        composition [28, 29]).
    :n: int (default=25)
        terminal block length used by ``region``.
    :normalize: str | None (default=None)
        None for the percentage composition, or "dde" for the dipeptide deviation
        from expected mean (needs at least 2 residues).

    Returns
    =======
    :dipeptide_composition_df: pd.DataFrame
        pandas dataframe of dipeptide composition for protein sequence. Dataframe will
        be of the shape 1 x 400, where 400 is the number of features calculated
        from the descriptor (20^2 for the 20 canonical amino acids) and 1 is the
        input sequence.
    """
    #validate input protein sequence
    sequence = _validate_sequence(sequence)
    if normalize not in (None, "dde"):
        raise ValueError(f"Invalid normalize '{normalize}'. Valid values are: [None, 'dde'].")
    if region is not None:
        return _regional(lambda part: dipeptide_composition(part, normalize=normalize), sequence, region, n)

    #count overlapping dipeptides with a sliding window (str.count misses overlaps, e.g. "AA" in "AAA")
    total_dipeptides = len(sequence) - 1
    counts = Counter(sequence[i:i+2] for i in range(total_dipeptides))

    if normalize == "dde":
        if total_dipeptides < 1:
            raise ValueError("DDE needs a sequence of at least 2 residues.")
        codons = _load_json_table("residue_tables.json", "codon_counts")["values"]
        values, columns = [], []
        for x, y in itertools.product(amino_acids, repeat=2):
            tm = (codons[x] / 61) * (codons[y] / 61)
            tv = tm * (1 - tm) / total_dipeptides
            values.append((counts[x + y] / total_dipeptides - tm) / math.sqrt(tv))
            columns.append(f"DDE_{x}{y}")
        return _to_df(values, columns)

    dipep_composition = {}

    #iterate through each amino acid, calculating each dipeptide combination's composition
    for i in amino_acids:
        for j in amino_acids:
            dipep = i + j
            dipep_composition[dipep] = counts[dipep] / total_dipeptides * 100 if total_dipeptides > 0 else 0.0

    #transform values and columns to DataFrame
    dipeptide_composition_df = pd.DataFrame([list(dipep_composition.values())], 
        columns=list(dipep_composition.keys()))

    return dipeptide_composition_df

def _regional(func, sequence: str, region: str, n: int) -> pd.DataFrame:
    """
    Calculate a composition descriptor over part of a protein sequence. The
    N-terminal ("nterm"), C-terminal ("cterm") or N-terminal, middle and
    C-terminal ("split") regions are passed to the descriptor function and the
    resulting columns are prefixed with NT_, MID_ and CT_, respectively. Split
    composition was used for predicting mitochondrial proteins and antibacterial
    peptides [28, 29].

    Parameters
    ==========
    :func: function
        composition descriptor function to apply to each region, e.g.
        amino_acid_composition.
    :sequence: str
        validated protein sequence.
    :region: str
        region of the sequence to use: "nterm", "cterm" or "split".
    :n: int
        length of the N- and C-terminal regions.

    Returns
    =======
    :regional_df: pd.DataFrame
        pandas dataframe of the descriptor values for each region, with the
        region prefix added to each column. For "split" the dataframe will be
        of the shape 1 x 3N, where N is the number of features of the descriptor.
    """
    n = _validate_int("n", n, 1)
    parts = {"nterm": [("NT_", sequence[:n])], "cterm": [("CT_", sequence[-n:])],
             "split": [("NT_", sequence[:n]), ("MID_", sequence[n:-n]), ("CT_", sequence[-n:])]}
    if region not in parts:
        raise ValueError(f"Invalid region '{region}'. Valid regions are: None, {list(parts)}.")
    #the middle block of a split must be non-empty
    if region == "split" and len(sequence) <= 2 * n:
        raise ValueError(f"region='split' needs a sequence longer than 2n ({2 * n}), got length {len(sequence)}.")
    frames = []
    for prefix, part in parts[region]:
        df = func(part)
        frames.append(df.add_prefix(prefix))
    return pd.concat(frames, axis=1)

################################## Tripeptide Composition #################################

def tripeptide_composition(sequence: str) -> pd.DataFrame:
    """
    Calculate Tripeptide Composition (TPComp) of protein sequence.
    Tripeptide composition is the fraction of each tripeptide type within a
    protein sequence. With tripeptides being of length 3 and there being 20
    canonical amino acids this creates 20^3 different combinations, thus a
    8000-Dimensional vector will be produced such that:

    TPComp(s,t,u) = AA(s,t,u) / (N - 2) * 100

    where TPComp(s,t,u) is the tripeptide composition of the protein sequence
    for amino acid type s, t and u (where s, t and u = 1,2,..,20), AA(s,t,u) is
    the number of tripeptides represented by amino acid type s, t and u, and N - 2 is
    the total number of tripeptides [1].

    This is a permanent alias of ``kmer_composition(sequence, k=3)``. It is not
    part of any ``all_descriptors`` set (8,000 mostly near-zero columns); request
    it with ``include=["tripeptide_composition"]``.

    Parameters
    ==========
    :sequence: str
        protein sequence in str form.

    Returns
    =======
    :tripeptide_composition_df: pd.DataFrame
        pandas DataFrame of tripeptide composition for protein sequence. Dataframe will
        be of the shape 1 x 8000, where 8000 is the number of features calculated
        from the descriptor (20^3 for the 20 canonical amino acids) and 1 is the
        input sequence. Values are percentages (0–100), consistent with
        dipeptide_composition.
    """
    return kmer_composition(sequence, k=3)

################################## GRAVY ##################################################

#Kyte-Doolittle hydropathy scale for all 20 standard amino acids [8] (protpy/data/aa_scales.json)
kyte_doolittle_ = dict(_read_json("aa_scales.json")["kyte_doolittle"]["values"])

def gravy(sequence: str, scale: str = "kyte_doolittle") -> pd.DataFrame:
    """
    Calculate the Grand Average of Hydropathy (GRAVY) for a protein sequence.
    GRAVY is the sum of the Kyte-Doolittle hydropathy values of all residues
    divided by the sequence length [8]:

    GRAVY = (sum of hydropathy values) / N

    A positive GRAVY value indicates an overall hydrophobic protein; a negative
    value indicates an overall hydrophilic protein.

    Another single-value scale can be averaged with ``scale`` ("hopp_woods",
    "eisenberg", "wimley_white", "engelman", "domcut" or any AAindex1 accession
    code); the column is then named ``GRAVY_{scale}``.

    Parameters
    ==========
    :sequence: str
        protein sequence.
    :scale: str (default="kyte_doolittle")
        hydrophobicity scale to average.

    Returns
    =======
    :gravy_df: pd.DataFrame
        pandas DataFrame containing the GRAVY score for the input sequence.
        Dataframe will be of the shape 1 x 1.
    """
    #validate input protein sequence
    sequence = _validate_sequence(sequence)

    if scale != "kyte_doolittle":
        values = _property_matrix(sequence, scale)
        if values.shape[1] != 1:
            raise ValueError(f"Scale '{scale}' is multi-dimensional; use a single-value scale.")
        return _to_df([values.mean()], [f"GRAVY_{scale}"])

    #sum Kyte-Doolittle hydropathy values and divide by sequence length
    gravy_score = sum(kyte_doolittle_[aa] for aa in sequence) / len(sequence)

    #transform to DataFrame
    gravy_df = pd.DataFrame([[gravy_score]], columns=["GRAVY"])

    return gravy_df

################################## Aromaticity ############################################

def aromaticity(sequence: str) -> pd.DataFrame:
    """
    Calculate Aromaticity of a protein sequence. Aromaticity is the fraction of
    aromatic amino acids (Phe, Trp, Tyr) in the sequence [9]:

    Aromaticity = (F + W + Y) / N

    where F, W and Y are the counts of Phe, Trp and Tyr residues, and N is the
    total sequence length.

    Parameters
    ==========
    :sequence: str
        protein sequence.

    Returns
    =======
    :aromaticity_df: pd.DataFrame
        pandas DataFrame of shape 1 x 1 containing the aromaticity score.
    """
    #validate input protein sequence
    sequence = _validate_sequence(sequence)

    #count aromatic residues (F, W, Y) and divide by sequence length
    aromatic_count = sum(sequence.count(aa) for aa in ("F", "W", "Y"))
    aromaticity_score = aromatic_count / len(sequence)

    #transform to DataFrame
    aromaticity_df = pd.DataFrame([[aromaticity_score]], columns=["Aromaticity"])

    return aromaticity_df

############################### Instability Index #########################################

#DIWV (dipeptide instability weight values) of Guruprasad et al. (1990) [10] are in
#protpy/data/residue_tables.json ("diwv"); they match Biopython's ProtParamData.DIWV

def instability_index(sequence: str) -> pd.DataFrame:
    """
    Calculate the Instability Index (II) of a protein sequence. The instability
    index is a measure of in vivo stability, computed as a weighted sum of
    dipeptide instability weight values (DIWV) [10]:

    II = (10 / N) * sum(DIWV(x_i, x_{i+1})) for i in 1..N-1

    Proteins with II < 40 are predicted to be stable; II >= 40 indicates
    instability. The DIWV table matches Biopython's ProtParamData (fixed in
    v1.6.0; earlier versions had 25 wrong entries).

    Parameters
    ==========
    :sequence: str
        protein sequence.

    Returns
    =======
    :instability_index_df: pd.DataFrame
        pandas DataFrame of shape 1 x 1 containing the instability index.
    """
    #validate input protein sequence
    sequence = _validate_sequence(sequence)

    #sum DIWV values for all consecutive dipeptides
    diwv = _read_json("residue_tables.json")["diwv"]["values"]
    dipeptide_sum = sum(diwv[sequence[i:i+2]] for i in range(len(sequence) - 1))
    ii = (dipeptide_sum * 10.0) / len(sequence)

    #transform to DataFrame
    instability_index_df = pd.DataFrame([[ii]], columns=["InstabilityIndex"])

    return instability_index_df

################################# Isoelectric Point ######################################

def _pka_set(pka_scale: str) -> dict:
    """
    Get a set of pKa values used in the isoelectric point and charge
    calculations from the residue_tables.json data file.

    Parameters
    ==========
    :pka_scale: str
        name of the pKa set: "protpy", "emboss", "bjellqvist" or "ipc_protein".

    Returns
    =======
    :pka_set: dict
        dict of pKa values for each ionisable group.
    """
    sets = _load_json_table("residue_tables.json", "pka_sets")["values"]
    if pka_scale not in sets:
        raise ValueError(f"Invalid pka_scale '{pka_scale}'. Valid scales are: {list(sets)}.")
    return sets[pka_scale]

def isoelectric_point(sequence: str, pka_scale: str = "bjellqvist") -> pd.DataFrame:
    """
    Calculate the theoretical Isoelectric Point (pI) of a protein sequence.
    pI is the pH at which the net charge of the protein is zero. It is computed
    by bisecting the pH range 0-14 until positive and negative charges balance [11].

    ``pka_scale`` selects the pKa values: "bjellqvist" (default since v1.6.0;
    Bjellqvist 1993/1994, including the residue-specific N-terminal pKa values for
    A, M, S, P, T, V, E and C-terminal values for D and E, as used by ExPASy
    ProtParam and Biopython), "emboss" (EMBOSS Epk.dat), "ipc_protein" (Kozlowski
    2016) or "protpy" (the unsourced v1.0-v1.5 default, kept for reproducing old
    results). The same charge model (K, R, H, N-terminus positive; D, E, C, Y,
    C-terminus negative) is used for every scale.

    Parameters
    ==========
    :sequence: str
        protein sequence.
    :pka_scale: str (default="bjellqvist")
        pKa value set: "bjellqvist", "emboss", "ipc_protein" or "protpy".

    Returns
    =======
    :isoelectric_point_df: pd.DataFrame
        pandas DataFrame of shape 1 x 1 containing the pI value.
    """
    #validate input protein sequence
    sequence = _validate_sequence(sequence)

    pka = dict(_pka_set(pka_scale))
    #residue-specific terminal pKa values (Bjellqvist)
    pka["N_term"] = pka.get("N_term_by_residue", {}).get(sequence[0], pka["N_term"])
    pka["C_term"] = pka.get("C_term_by_residue", {}).get(sequence[-1], pka["C_term"])

    #count ionisable residues
    counts = {aa: sequence.count(aa) for aa in "DEHCYKR"}

    def _net_charge(ph: float) -> float:
        pos = (counts["K"] / (1 + 10 ** (ph - pka["K"])) +
               counts["R"] / (1 + 10 ** (ph - pka["R"])) +
               counts["H"] / (1 + 10 ** (ph - pka["H"])) +
               1.0 / (1 + 10 ** (ph - pka["N_term"])))
        neg = (counts["D"] / (1 + 10 ** (pka["D"] - ph)) +
               counts["E"] / (1 + 10 ** (pka["E"] - ph)) +
               counts["C"] / (1 + 10 ** (pka["C"] - ph)) +
               counts["Y"] / (1 + 10 ** (pka["Y"] - ph)) +
               1.0 / (1 + 10 ** (pka["C_term"] - ph)))
        return pos - neg

    #net charge decreases monotonically with pH, so bisect over the full 0-14 range
    low, high = 0.0, 14.0
    for _ in range(50):
        mid = (low + high) / 2
        if _net_charge(mid) > 0:
            low = mid
        else:
            high = mid
    ph = (low + high) / 2

    #transform to DataFrame
    isoelectric_point_df = pd.DataFrame([[ph]], columns=["IsoelectricPoint"])

    return isoelectric_point_df

############################### Molecular Weight ##########################################

#average residue molecular weights (Da) for each of the 20 standard amino acids [12]
residue_molecular_weight_ = {
    "A": 89.094,  "C": 121.154, "D": 133.104, "E": 147.130, "F": 165.192,
    "G": 75.067,  "H": 155.156, "I": 131.174, "K": 146.189, "L": 131.174,
    "M": 149.208, "N": 132.119, "P": 115.132, "Q": 146.146, "R": 174.203,
    "S": 105.093, "T": 119.119, "V": 117.148, "W": 204.228, "Y": 181.191
}

def molecular_weight(sequence: str) -> pd.DataFrame:
    """
    Calculate the Molecular Weight (MW) of a protein sequence. MW is computed
    as the sum of average residue masses minus one water molecule per peptide
    bond (18.015 Da) [12]:

    MW = sum(residue_mass_i) - (N - 1) * 18.015

    Parameters
    ==========
    :sequence: str
        protein sequence.

    Returns
    =======
    :molecular_weight_df: pd.DataFrame
        pandas DataFrame of shape 1 x 1 containing the molecular weight (Da).
    """
    #validate input protein sequence
    sequence = _validate_sequence(sequence)

    #sum residue masses and subtract one water per peptide bond
    mw = sum(residue_molecular_weight_[aa] for aa in sequence) - (len(sequence) - 1) * 18.015
    molecular_weight_df = pd.DataFrame([[mw]], columns=["MolecularWeight"])

    return molecular_weight_df

############################## Charge Distribution ########################################

def charge_distribution(sequence: str, ph: float = 7.4, pka_scale: str = "bjellqvist",
    normalize: str = "none") -> pd.DataFrame:
    """
    Calculate Charge Distribution of a protein sequence at a given pH. Positive
    charge is contributed by Lys (K), Arg (R) and His (H); negative charge by
    Asp (D) and Glu (E). Henderson-Hasselbalch equations are used [13]:

    positive = sum(count(aa) / (1 + 10^(pH - pKa))) for aa in {K, R, H}
    negative = sum(count(aa) / (1 + 10^(pKa - pH))) for aa in {D, E}
    net      = positive - negative

    The charges are totals over the sequence. ``normalize="length"`` divides them
    by the sequence length (charge per residue) and ``normalize="mass"`` by the
    molecular weight (charge per Da, the charge density used by modlAMP).

    This is a side-chain Henderson-Hasselbalch model at a given pH. It differs on
    purpose from ``charge_fractions`` (integer K/R vs D/E counts, the Das & Pappu
    convention) and from ``isoelectric_point`` (which also includes C, Y and the
    termini); the pKa values default to the same Bjellqvist set as
    ``isoelectric_point``.

    Parameters
    ==========
    :sequence: str
        protein sequence.
    :ph: float (default=7.4)
        pH at which to evaluate charges, between 0 and 14.
    :pka_scale: str (default="bjellqvist")
        pKa value set for K, R, H, D and E: "bjellqvist", "emboss", "ipc_protein"
        or "protpy" (see ``isoelectric_point``). Termini are not included.
    :normalize: str (default="none")
        "none" (total charge), "length" (per residue) or "mass" (per Da).

    Returns
    =======
    :charge_distribution_df: pd.DataFrame
        pandas DataFrame of shape 1 x 3 with columns
        ["PositiveCharge", "NegativeCharge", "NetCharge"], suffixed "_perResidue"
        or "_perDa" when normalized.
    """
    #validate input protein sequence
    sequence = _validate_sequence(sequence)
    if isinstance(ph, bool) or not isinstance(ph, (int, float)) or not 0 <= ph <= 14:
        raise ValueError(f"ph must be a number between 0 and 14, got {ph}.")
    suffixes = {"none": "", "length": "_perResidue", "mass": "_perDa"}
    if normalize not in suffixes:
        raise ValueError(f"Invalid normalize '{normalize}'. Valid values are: {list(suffixes)}.")

    #pKa values for charged residues
    pka_charged_ = _pka_set(pka_scale)

    #compute positive charge (K, R, H) using Henderson-Hasselbalch
    pos = sum(sequence.count(aa) / (1 + 10 ** (ph - pka_charged_[aa])) for aa in ("K", "R", "H"))
    #compute negative charge (D, E)
    neg = sum(sequence.count(aa) / (1 + 10 ** (pka_charged_[aa] - ph)) for aa in ("D", "E"))

    divisor = 1
    if normalize == "length":
        divisor = len(sequence)
    elif normalize == "mass":
        divisor = molecular_weight(sequence)["MolecularWeight"].iloc[0]
    suffix = suffixes[normalize]
    return _to_df([pos / divisor, neg / divisor, (pos - neg) / divisor],
                  [f"PositiveCharge{suffix}", f"NegativeCharge{suffix}", f"NetCharge{suffix}"])

####################### Hydrophobic/Polar/Charged Composition ############################

def hydrophobic_polar_charged_composition(sequence: str) -> pd.DataFrame:
    """
    Calculate Hydrophobic, Polar and Charged Composition of a protein sequence.
    Residues are grouped into three physicochemical classes [1]:

    - Hydrophobic (nonpolar): A, C, F, I, L, M, V, W, Y
    - Polar (uncharged):       G, N, Q, S, T
    - Charged:                 D, E, R, H, K, P

    Each value is expressed as a percentage of the total sequence length.

    Deprecated since v1.6.0, removed in 2.0: use
    ``kmer_composition(sequence, k=1, alphabet=3)``, which uses the same three
    groups. Before v1.6.0 proline was left out of every class, so the three
    values summed to 100 - P%.

    Parameters
    ==========
    :sequence: str
        protein sequence.

    Returns
    =======
    :hpc_df: pd.DataFrame
        pandas DataFrame of shape 1 x 3 with columns
        ["Hydrophobic", "Polar", "Charged"].
    """
    _warn_deprecated("hydrophobic_polar_charged_composition", "kmer_composition(sequence, k=1, alphabet=3)")
    hpc_df = kmer_composition(sequence, k=1, alphabet=3)
    hpc_df.columns = ["Hydrophobic", "Polar", "Charged"]
    return hpc_df

########################## Secondary Structure Propensity ################################

#Chou-Fasman propensity values for alpha-helix, beta-sheet and coil [14]
helix_propensity_ = {
    "A": 1.45, "C": 0.77, "D": 0.98, "E": 1.53, "F": 1.12, "G": 0.53, "H": 1.24,
    "I": 1.00, "K": 1.07, "L": 1.34, "M": 1.20, "N": 0.73, "P": 0.59, "Q": 1.17,
    "R": 0.79, "S": 0.79, "T": 0.82, "V": 1.14, "W": 1.14, "Y": 0.61
}
sheet_propensity_ = {
    "A": 0.97, "C": 1.30, "D": 0.80, "E": 0.26, "F": 1.28, "G": 0.81, "H": 0.71,
    "I": 1.60, "K": 0.74, "L": 1.22, "M": 1.67, "N": 0.65, "P": 0.62, "Q": 1.23,
    "R": 0.90, "S": 0.72, "T": 1.20, "V": 1.65, "W": 1.19, "Y": 1.29
}
coil_propensity_ = {
    "A": 0.72, "C": 1.17, "D": 1.41, "E": 0.80, "F": 0.58, "G": 1.64, "H": 1.00,
    "I": 0.51, "K": 1.01, "L": 0.57, "M": 0.67, "N": 1.68, "P": 1.91, "Q": 0.98,
    "R": 1.24, "S": 1.33, "T": 1.03, "V": 0.50, "W": 0.99, "Y": 1.29
}

def secondary_structure_propensity(sequence: str) -> pd.DataFrame:
    """
    Calculate Secondary Structure Propensity of a protein sequence. Each residue
    is assigned Chou-Fasman propensity values for alpha-helix (H), beta-sheet (E)
    and coil/loop (C) conformations [14]. The mean propensity for each secondary
    structure class across the whole sequence is returned.

    Deprecated since v1.6.0, removed in 2.0: the three values are fixed linear
    maps of the amino acid composition, and the tables match no AAindex record.
    Average a Chou-Fasman AAindex record instead, e.g.
    ``gravy(sequence, scale="CHOP780201")``.

    Parameters
    ==========
    :sequence: str
        protein sequence.

    Returns
    =======
    :ssp_df: pd.DataFrame
        pandas DataFrame of shape 1 x 3 with columns ["Helix", "Sheet", "Coil"].
    """
    _warn_deprecated("secondary_structure_propensity", 'gravy(sequence, scale="<AAindex code>")')
    #validate input protein sequence
    sequence = _validate_sequence(sequence)

    n = len(sequence)
    #average propensity per secondary structure class
    h = sum(helix_propensity_[aa] for aa in sequence) / n
    e = sum(sheet_propensity_[aa] for aa in sequence) / n
    c = sum(coil_propensity_[aa] for aa in sequence) / n

    return _to_df([h, e, c], ["Helix", "Sheet", "Coil"])

################################# k-mer Composition ######################################

#predefined reduced alphabets mapping 20 amino acids to fewer groups
#key = alphabet size, value = list of frozensets (each frozenset = one group)
reduced_alphabets_ = {
    2: [frozenset("ACFILMVWY"), frozenset("DEGNHKPQRST")],
    3: [frozenset("ACFILMVWY"), frozenset("GNQST"), frozenset("DERHKP")],
    4: [frozenset("ACFILMVW"), frozenset("GNQSTY"), frozenset("DERHK"), frozenset("P")],
    6: [
        frozenset("ACST"),   # small polar/aliphatic
        frozenset("FILMVWY"),# hydrophobic/aromatic
        frozenset("DE"),      # acidic
        frozenset("KRH"),     # basic
        frozenset("NQ"),      # amide
        frozenset("GP"),      # special
    ]
}

#names of the reduced alphabets in protpy/data/reduced_alphabets.json
_named_alphabets = ["ifeature5", "cp13", "cp14", "cp19", "cp20"]

def _alphabet(alphabet: int | str | list[str]) -> tuple[str, list[tuple[str, str]]]:
    """
    Get the groups of a reduced amino acid alphabet. The alphabet can be the size
    of one of protpy's reduced alphabets (2, 3, 4 or 6), the name of an alphabet in
    the reduced_alphabets.json data file (the 5-group iFeature alphabet "ifeature5"
    [21] and the cp13, cp14, cp19 and cp20 cluster profiles [24]) or a user-supplied
    list of residue strings that together contain each amino acid exactly once,
    e.g. ["AGV", "ILFP", "YMTS", "HNQW", "RK", "DE", "C"] (any PseKRAAC grouping can
    be given this way).

    Parameters
    ==========
    :alphabet: int | str | list[str]
        alphabet size, alphabet name or list of residue groups.

    Returns
    =======
    :prefix: str
        column prefix for the alphabet.
    :groups: list
        list of (label, residues) tuples, one per group of the alphabet.
    """
    if isinstance(alphabet, (list, tuple)):
        residues = "".join(alphabet).upper()
        if sorted(residues) != sorted(amino_acids):
            raise ValueError("A custom alphabet must be a list of residue strings that together contain each "
                             f"of the 20 amino acids exactly once, got {list(alphabet)}.")
        return "RA_custom", [(group.upper(), group.upper()) for group in alphabet]
    if isinstance(alphabet, str):
        if alphabet not in _named_alphabets:
            raise ValueError(f"Invalid alphabet '{alphabet}'. Valid alphabets are: {_named_alphabets}, "
                             f"a size in {list(reduced_alphabets_)} or a list of residue groups.")
        values = _load_json_table("reduced_alphabets.json", alphabet)["values"]
        groups = list(values.items()) if isinstance(values, dict) else [(group, group) for group in values]
        return f"RA_{alphabet}", groups
    if isinstance(alphabet, bool) or alphabet not in reduced_alphabets_:
        raise ValueError(f"Invalid alphabet size {alphabet!r}. Valid sizes are: {list(reduced_alphabets_)}.")
    return "ReducedAlphabet", [(str(idx + 1), "".join(sorted(group)))
                               for idx, group in enumerate(reduced_alphabets_[alphabet])]

def kmer_composition(sequence: str, k: int = 2, alphabet: int | str | list[str] | None = None) -> pd.DataFrame:
    """
    Calculate k-mer Composition of a protein sequence. A k-mer is a subsequence
    of length k. For each of the 20^k possible k-mers the percentage of the
    total number of overlapping k-mers in the sequence is calculated [1]:

    kmer(x1..xk) = N(x1..xk) / (L - k + 1) * 100

    k=1 is the amino acid composition, k=2 the dipeptide composition and k=3 the
    tripeptide composition.

    With ``alphabet`` the amino acids are first mapped to the groups of a reduced
    alphabet, giving the reduced alphabet composition of G^k group k-mers, where G
    is the number of groups. ``alphabet`` can be a size (2, 3, 4 or 6; columns
    ReducedAlphabet_*), a named alphabet ("ifeature5", "cp13", "cp14", "cp19",
    "cp20"; columns RA_<name>_*) or a list of residue strings (columns
    RA_custom_*). With ``alphabet="ifeature5"`` (aliphatic GAVLMI, aromatic FYW,
    positive KRH, negative DE, uncharged STCPNQ) this gives grouped amino acid
    composition (GAAC, k=1), grouped dipeptide composition (GDPC, k=2) and grouped
    tripeptide composition (GTPC, k=3), with those column prefixes [21].
    ``alphabet=3`` is the hydrophobic (ACFILMVWY), polar (GNQST) and charged
    (DERHKP) composition.

    Parameters
    ==========
    :sequence: str
        protein sequence.
    :k: int (default=2)
        length of the k-mers, 1-3 (k=4 would give 160,000 almost entirely zero
        columns).
    :alphabet: int | str | list[str] | None (default=None)
        reduced alphabet; None uses the 20 amino acids.

    Returns
    =======
    :kmer_composition_df: pd.DataFrame
        pandas DataFrame of shape 1 x 20^k (1 x G^k with an alphabet) containing
        k-mer percentages. Column names are the k-mers (e.g. "AC"), or the group
        labels joined by "." with an alphabet prefix (e.g. "GDPC_aliphatic.aromatic").
    """
    #validate input protein sequence
    sequence = _validate_sequence(sequence)

    #validate k; raise error for impractically large k
    k = _validate_int("k", k, 1)
    if k > 3:
        raise ValueError(f"k must be 1, 2 or 3, got {k} ({20**k:,} columns).")
    total = len(sequence) - k + 1

    if alphabet is None:
        #count every (overlapping) k-mer in one pass, divide by total k-mer count
        counts = Counter(sequence[i:i+k] for i in range(max(total, 0)))
        kmers = ["".join(p) for p in itertools.product(amino_acids, repeat=k)]
        return _to_df([counts[kmer] / total * 100 if total > 0 else 0.0 for kmer in kmers], kmers)

    prefix, groups = _alphabet(alphabet)
    if alphabet == "ifeature5":
        prefix = {1: "GAAC", 2: "GDPC", 3: "GTPC"}[k]
    mapping = {aa: label for label, residues in groups for aa in residues}
    counts = _kmer_counts(sequence, k, mapping)
    combos = list(itertools.product([label for label, _ in groups], repeat=k))
    values = [counts[c] / total * 100 if total > 0 else 0.0 for c in combos]
    return _to_df(values, [f"{prefix}_" + ".".join(c) for c in combos])

########################## Reduced Alphabet Composition ##################################

def reduced_alphabet_composition(sequence: str, alphabet_size: int = 6, alphabet: str | None = None,
    k: int = 1) -> pd.DataFrame:
    """
    Calculate Reduced Alphabet Composition of a protein sequence. The 20 standard
    amino acids are clustered into a smaller alphabet of physicochemically similar
    groups and the percentage of each group k-mer is returned [1].

    Deprecated since v1.6.0, removed in 2.0: use
    ``kmer_composition(sequence, k, alphabet=alphabet_size or alphabet)``, which
    gives the same output. Since v1.6.0 an unsupported ``alphabet_size`` raises a
    ValueError instead of silently using 6.

    Parameters
    ==========
    :sequence: str
        protein sequence.
    :alphabet_size: int (default=6)
        number of amino acid groups (2, 3, 4 or 6); ignored when ``alphabet`` is set.
    :alphabet: str | None (default=None)
        named alphabet: "ifeature5", "cp13", "cp14", "cp19" or "cp20".
    :k: int (default=1)
        group k-mer length, 1-3.

    Returns
    =======
    :reduced_alphabet_df: pd.DataFrame
        pandas DataFrame of shape 1 x G^k, where G is the number of groups. Values
        are percentages.
    """
    _warn_deprecated("reduced_alphabet_composition", "kmer_composition(sequence, k, alphabet=...)")
    return kmer_composition(sequence, k=k, alphabet=alphabet if alphabet is not None else alphabet_size)

############################## Motif Composition ##########################################

#default set of biologically relevant motifs with regex patterns [15]; KDEL is a
#C-terminal signal, so it is anchored to the end of the sequence (fixed in v1.6.0)
default_motifs_ = {
    "NxS/T_glycosylation": r"N[^P][ST]",
    "RGD_integrin":        r"RGD",
    "KDEL_retention":      r"KDEL$",
    "CxxC_zinc_finger":    r"C.{2}C",
    "CAAX_prenylation":    r"C[ACFILMVWY]{2}[ACEILMQV]$",
    "cAMP_PKA":            r"RR.S",
    "dileucine_sorting":   r"[DE]...L[LI]",
    "PEST_region":         r"P.{1,10}[ESD].{1,10}T"
}

#PROSITE patterns [15] and the trypsin rule, selectable by name in motif_composition
prosite_motifs_ = {
    "N_glycosylation": r"N[^P][ST][^P]",       #PS00001 N-{P}-[ST]-{P}
    "ER_retention":    r"[KRHQSA][DENQ]EL$",   #PS00014 [KRHQSA]-[DENQ]-E-L>
    "PTS1":            r"[STAGCN][RKH][LIVMAFY]$",  #PS00342 [STAGCN]-[RKH]-[LIVMAFY]>
    "tryptic_site":    r"[KR](?=[^P])",        #after K/R, not before P or at the C-terminus
}

#every built-in motif, by name
motif_library_ = {**default_motifs_, **prosite_motifs_}

def motif_composition(sequence: str, motifs: dict[str, str] | list[str] | None = None,
    per_100: bool = False) -> pd.DataFrame:
    """
    Calculate Motif Composition of a protein sequence. For each motif the number
    of (overlapping) regex matches within the sequence is counted [15]. If no
    motifs are supplied, a set of 8 biologically relevant default motifs is used.

    The built-in library also contains the PROSITE patterns N_glycosylation
    (PS00001, N-{P}-[ST]-{P}), ER_retention (PS00014, [KRHQSA]-[DENQ]-E-L at the
    C-terminus) and PTS1 (PS00342, [STAGCN]-[RKH]-[LIVMAFY] at the C-terminus),
    and tryptic_site (after K or R, not before P or at the C-terminus); select
    them by name, e.g. ``motifs=["N_glycosylation", "PTS1"]``.

    Parameters
    ==========
    :sequence: str
        protein sequence.
    :motifs: dict | list | None (default=None)
        dict mapping motif names to regex patterns, or a list whose entries are
        built-in motif names or regex patterns (a pattern's column is named by the
        pattern itself). None uses the 8 default motifs. An invalid pattern raises
        a ValueError. Patterns run with Python's ``re`` module, so a badly written
        pattern such as ``(A+)+B`` can take very long on a long sequence; do not
        pass patterns from untrusted sources.
    :per_100: bool (default=False)
        also return each count per 100 residues, in columns suffixed "_per100".

    Returns
    =======
    :motif_composition_df: pd.DataFrame
        pandas DataFrame of shape 1 x M (1 x 2M with ``per_100``), where M is the
        number of motifs. Each value is the integer count of matches for that motif.
    """
    #validate input protein sequence
    sequence = _validate_sequence(sequence)

    #use default motifs if none provided
    if motifs is None:
        motifs = default_motifs_
    elif isinstance(motifs, (list, tuple)):
        motifs = {m: motif_library_.get(m, m) for m in motifs}
    elif not isinstance(motifs, dict):
        raise TypeError(f"motifs must be a dict, a list or None, got {type(motifs)}.")

    #count overlapping matches for each motif pattern
    compiled = []
    for name, pattern in motifs.items():
        try:
            compiled.append(re.compile(r'(?=' + pattern + r')'))
        except (re.error, TypeError) as error:
            raise ValueError(f"Invalid regex pattern for motif '{name}': {pattern!r} ({error}).") from None
    counts = [len(regex.findall(sequence)) for regex in compiled]
    columns = list(motifs)
    if per_100:
        counts += [c / len(sequence) * 100 for c in counts]
        columns += [f"{name}_per100" for name in motifs]
    return _to_df(counts, columns)

######################### Amino Acid Pair Composition #####################################

def amino_acid_pair_composition(sequence: str) -> pd.DataFrame:
    """
    Calculate Amino Acid Pair Composition (PairComp) of a protein sequence.
    For each pair of consecutive residues (i, i+1), the fractional frequency
    of all 400 possible dipeptide types (20 x 20) is computed [1]:

    PairComp(s,t) = count(s,t) / (N - 1) * 100

    Unlike the standard dipeptide composition, results are annotated with the
    physicochemical class of each residue in the pair (Hydrophobic/Polar/Charged,
    the groups of ``kmer_composition(alphabet=3)``, so P is Charged).

    Deprecated since v1.6.0, removed in 2.0: the values are exactly the dipeptide
    composition; use ``dipeptide_composition``.

    Parameters
    ==========
    :sequence: str
        protein sequence.

    Returns
    =======
    :pair_composition_df: pd.DataFrame
        pandas DataFrame of shape 1 x 400, where each column label has the form
        'XY_Class1-Class2' (e.g. 'AL_Hydrophobic-Hydrophobic').
    """
    _warn_deprecated("amino_acid_pair_composition", "dipeptide_composition")
    df = dipeptide_composition(sequence)
    classes = {aa: name for name, (_, residues) in zip(("Hydrophobic", "Polar", "Charged"), _alphabet(3)[1])
               for aa in residues}
    df.columns = [f"{c}_{classes[c[0]]}-{classes[c[1]]}" for c in df.columns]
    return df

################################## Aliphatic Index ########################################

def aliphatic_index(sequence: str) -> pd.DataFrame:
    """
    Calculate the Aliphatic Index (AI) of a protein sequence. The aliphatic
    index is a measure of thermostability, defined as the relative volume
    occupied by aliphatic side chains (Ala, Val, Ile, Leu) [16]:

    AI = Xala + 2.9 * Xval + 3.9 * (Xile + Xleu)

    where X values are mole percentages of each residue. Higher values
    indicate greater thermostability.

    Parameters
    ==========
    :sequence: str
        protein sequence.

    Returns
    =======
    :aliphatic_index_df: pd.DataFrame
        pandas DataFrame of shape 1 x 1 containing the aliphatic index score.
    """
    #validate input protein sequence
    sequence = _validate_sequence(sequence)

    n = len(sequence)
    #compute mole percentages for aliphatic residues
    xa = sequence.count("A") / n * 100
    xv = sequence.count("V") / n * 100
    xi = sequence.count("I") / n * 100
    xl = sequence.count("L") / n * 100

    #apply Ikai (1980) aliphatic index formula
    ai = xa + 2.9 * xv + 3.9 * (xi + xl)

    aliphatic_index_df = pd.DataFrame([[ai]], columns=["AliphaticIndex"])

    return aliphatic_index_df

############################### Extinction Coefficient ####################################

def extinction_coefficient(sequence: str, absorbance: bool = False) -> pd.DataFrame:
    """
    Calculate the molar Extinction Coefficient of a protein sequence at 280 nm.
    The coefficient is derived from the number of Trp (W), Tyr (Y) and
    Cys (C) residues using the method of Gasteiger et al. (2005) [17]:

    EC = (W * 5500) + (Y * 1490) + (SS * 125)

    where SS is the number of disulfide bonds (Cys_count // 2) for the
    oxidised form, and 0 for the reduced form. Both are returned.

    With ``absorbance=True`` the A280 of a 1 g/L (0.1%) solution is added,
    Abs 0.1% = EC / MW, for both forms (Pace et al. 1995, Gill & von Hippel 1989).

    Parameters
    ==========
    :sequence: str
        protein sequence.
    :absorbance: bool (default=False)
        also return ``Abs01_Reduced`` and ``Abs01_Oxidized``.

    Returns
    =======
    :extinction_coefficient_df: pd.DataFrame
        pandas DataFrame of shape 1 x 2 (1 x 4 with ``absorbance``) with columns
        ``ExtCoeff_Reduced`` and ``ExtCoeff_Oxidized`` (M⁻¹cm⁻¹).
    """
    #validate input protein sequence
    sequence = _validate_sequence(sequence)

    w = sequence.count("W")
    y = sequence.count("Y")
    c = sequence.count("C")

    #reduced: no disulfide bonds
    ec_reduced = w * 5500 + y * 1490
    #oxidized: each pair of Cys contributes one disulfide bond (125 M⁻¹cm⁻¹)
    ec_oxidized = ec_reduced + (c // 2) * 125

    extinction_coefficient_df = pd.DataFrame(
        [[ec_reduced, ec_oxidized]],
        columns=["ExtCoeff_Reduced", "ExtCoeff_Oxidized"]
    )

    if absorbance:
        mw = molecular_weight(sequence)["MolecularWeight"].iloc[0]
        extinction_coefficient_df["Abs01_Reduced"] = ec_reduced / mw
        extinction_coefficient_df["Abs01_Oxidized"] = ec_oxidized / mw

    return extinction_coefficient_df

################################ Boman Index ##############################################

def boman_index(sequence: str) -> pd.DataFrame:
    """
    Calculate the Boman Index (potential protein interaction index) for
    a protein sequence. The index is the sum of the solubility values of
    all residues divided by the sequence length, negated [18]:

    BI = -sum(solubility(aa_i)) / N

    where the solubility values are those of Radzicka & Wolfenden (AAindex
    RADA880101), so BI = -gravy(sequence, scale="RADA880101").

    Higher values indicate a greater tendency to bind other proteins or
    cell membranes (e.g. antimicrobial peptides). Values >2.48 suggest
    high interaction potential.

    Parameters
    ==========
    :sequence: str
        protein sequence.

    Returns
    =======
    :boman_index_df: pd.DataFrame
        pandas DataFrame of shape 1 x 1 containing the Boman Index score.
    """
    #negated mean of the Radzicka-Wolfenden solubility values (AAindex RADA880101)
    return _to_df([-gravy(sequence, scale="RADA880101").iloc[0, 0]], ["BomanIndex"])

############################### Aggregation Propensity ####################################

def aggregation_propensity(sequence: str, window: int = 5, hydrophobicity_threshold: float = 2.0, charge_threshold: int = 0) -> pd.DataFrame:
    """
    Calculate the Aggregation Propensity of a protein sequence. Aggregation-
    prone regions (APRs) are identified using a sliding window: a window is
    classified as an APR when its mean Kyte-Doolittle hydrophobicity exceeds
    ``hydrophobicity_threshold`` and the number of charged residues (D, E, K,
    R) is at most ``charge_threshold`` [8].

    Parameters
    ==========
    :sequence: str
        protein sequence.
    :window: int (default=5)
        sliding window length in residues.
    :hydrophobicity_threshold: float (default=2.0)
        minimum mean Kyte-Doolittle hydrophobicity for an APR window.
    :charge_threshold: int (default=0)
        maximum number of charged residues (D, E, K, R) allowed in an APR.

    Deprecated since v1.6.0, removed in 2.0: this is an unpublished heuristic; use
    the published AGGRESCAN method, ``aggrescan``. Since v1.6.0 overlapping APR
    windows are merged, so ``AggregProneRegions`` counts regions, not windows.

    Returns
    =======
    :aggregation_propensity_df: pd.DataFrame
        pandas DataFrame of shape 1 x 2 with columns:
        ``AggregProneRegions`` — number of aggregation-prone regions (runs of
        residues covered by APR windows);
        ``AggregProneFraction`` — percentage of residues covered by any APR.
    """
    _warn_deprecated("aggregation_propensity", "aggrescan")
    #validate input protein sequence
    sequence = _validate_sequence(sequence)

    window = _validate_int("window", window, 1)

    n = len(sequence)
    charged = frozenset("DEKR")
    #track which residue positions are covered by at least one APR
    covered = np.zeros(n, dtype=bool)

    #score each window
    for i in range(n - window + 1):
        seg = sequence[i:i + window]
        avg_hydro = sum(kyte_doolittle_[aa] for aa in seg) / window
        n_charged = sum(1 for aa in seg if aa in charged)
        if avg_hydro >= hydrophobicity_threshold and n_charged <= charge_threshold:
            covered[i:i + window] = True

    return _to_df([len(_runs(covered)), covered.mean() * 100], ["AggregProneRegions", "AggregProneFraction"])

############################### Hydrophobic Moment ########################################

def hydrophobic_moment(sequence: str, window: int = 11, angle: float | str = 100,
    scale: str | None = None) -> pd.DataFrame:
    """
    Calculate the mean and maximum Hydrophobic Moment of a protein sequence.
    The hydrophobic moment measures the amphiphilicity of a helical segment —
    the directional asymmetry of hydrophobicity projected around a helical
    axis [19].

    For each window of length ``window``, the moment is:

        mu = sqrt( (sum H_i * sin(i * angle))^2 + (sum H_i * cos(i * angle))^2 ) / window

    where H_i is the hydrophobicity of residue i and angle is in degrees.
    Default angle of 100° corresponds to an alpha-helix; 160° is used for
    beta-strands (Eisenberg et al. 1984).

    By default H is the normalised consensus scale of Eisenberg, Schwarz, Komarony
    & Wall (1984, J Mol Biol 179:125-142; A = 0.62) as listed by ExPASy ProtScale,
    available by name as "eisenberg_jmb1984". ``scale`` selects another
    single-value scale, e.g. "eisenberg_arb1984" (AAindex EISD840101, Eisenberg 1984
    Annu Rev Biochem, A = 0.25; a different table, formerly called "eisenberg"),
    "kyte_doolittle" or any AAindex1 code.

    Parameters
    ==========
    :sequence: str
        protein sequence.
    :window: int (default=11)
        sliding window length in residues.
    :angle: float | str (default=100)
        residue rotation angle in degrees, or "helix" (100) / "sheet" (160).
    :scale: str | None (default=None)
        hydrophobicity scale; None uses "eisenberg_jmb1984".

    Returns
    =======
    :hydrophobic_moment_df: pd.DataFrame
        pandas DataFrame of shape 1 x 2 with columns
        ``HydrophobicMoment_Mean`` and ``HydrophobicMoment_Max``.
    """
    #validate input protein sequence
    sequence = _validate_sequence(sequence)

    window = _validate_int("window", window, 1)
    if isinstance(angle, str):
        if angle not in ("helix", "sheet"):
            raise ValueError(f"angle must be a number, 'helix' or 'sheet', got '{angle}'.")
        angle = {"helix": 100, "sheet": 160}[angle]
    if scale is None:
        h_scale = hydrophobicity_
    else:
        values = _property_matrix("".join(amino_acids), scale)
        if values.shape[1] != 1:
            raise ValueError(f"Scale '{scale}' is multi-dimensional; use a single-value scale.")
        h_scale = dict(zip(amino_acids, values[:, 0]))

    n = len(sequence)
    #use full sequence as single window when shorter than default window
    effective_window = min(window, n)
    angle_rad = angle * math.pi / 180.0
    moments = []

    #compute hydrophobic moment for each window position using Eisenberg scale
    for i in range(n - effective_window + 1):
        seg = sequence[i:i + effective_window]
        hsin = sum(h_scale[aa] * math.sin(j * angle_rad) for j, aa in enumerate(seg))
        hcos = sum(h_scale[aa] * math.cos(j * angle_rad) for j, aa in enumerate(seg))
        mu = math.sqrt(hsin ** 2 + hcos ** 2) / effective_window
        moments.append(mu)

    return _to_df([sum(moments) / len(moments), max(moments)], ["HydrophobicMoment_Mean", "HydrophobicMoment_Max"])

############################## Pseudo Amino Acid Composition ##############################

def _resolve_properties(properties, default: list[dict[str, float]]) -> list[dict[str, float]]:
    """
    Get the amino acid values of the properties used by PseAAC and APseAAC. Each
    property can be an AAindex1 accession number or a dict of amino acid values,
    so both descriptors accept either form.

    Parameters
    ==========
    :properties: str | dict | list | None
        AAindex1 accession number(s) and/or dicts of amino acid values; None or an
        empty list uses ``default``.
    :default: list[dict]
        default property dicts.

    Returns
    =======
    :property_values: list[dict]
        dict of amino acid values for each property.
    """
    if properties is None or (isinstance(properties, list) and not properties):
        return default
    if not isinstance(properties, list):
        properties = [properties]
    values = []
    for prop in properties:
        if isinstance(prop, dict):
            values.append(prop)
            continue
        if prop not in aaindex1.record_codes():
            raise ValueError(f"Accession number not found in aaindex1: {prop}.")
        #keep only the 20 canonical amino acids (drops the '-' gap entry)
        values.append({aa: v for aa, v in aaindex1[prop].values.items() if aa in amino_acids})
    return values

def _pseudo_composition(sequence: str, factors: np.ndarray, weight: float, prefix: str) -> pd.DataFrame:
    """
    Combine the amino acid composition with sequence-order correlation factors
    as in Chou's pseudo amino acid compositions [4, 5]:

    X_u = f_u / (1 + w * sum_j tau_j) * 100          for the 20 amino acids

    X_{20+j} = w * tau_j / (1 + w * sum_j tau_j) * 100   for each correlation factor

    where f_u is the fraction of amino acid u and tau_j the correlation factors.
    The values are percentages and sum to 100. Amphiphilic factors can be negative,
    so the denominator is checked: a RuntimeWarning is raised when it is below
    0.1 (the values are then dominated by noise) and a ValueError when it is not
    positive.

    Parameters
    ==========
    :sequence: str
        validated protein sequence.
    :factors: np.ndarray
        correlation factors tau_j.
    :weight: float
        weighting factor w.
    :prefix: str
        column prefix, e.g. "PAAC".

    Returns
    =======
    :pseudo_df: pd.DataFrame
        dataframe of shape 1 x (20 + number of factors).
    """
    denominator = 1 + weight * factors.sum()
    if denominator <= 0:
        raise ValueError(f"{prefix}: the normalisation term 1 + w * sum(tau) is {denominator:.3g} <= 0, so the "
                         "descriptor is undefined; use a smaller weight or lamda.")
    if denominator < 0.1:
        warnings.warn(f"{prefix}: the normalisation term 1 + w * sum(tau) is only {denominator:.3g}, so the values "
                      "are numerically unstable (this happens for short sequences with a large weight and lamda); "
                      "use a smaller weight or lamda.", RuntimeWarning, stacklevel=3)
    comp = np.array([sequence.count(aa) for aa in amino_acids]) / len(sequence)
    values = np.concatenate([comp, weight * factors]) / denominator * 100
    return _to_df(list(values), [f"{prefix}_{i + 1}" for i in range(len(values))])

def pseudo_amino_acid_composition(sequence: str, lamda: int = 30, weight: float = 0.05,
    properties: list[str | dict[str, float]] | None = None) -> pd.DataFrame:
    """
    Calculate Pseudo Amino Acid Composition (PAAComp) of protein sequence. PAAComp
    combines the amino acid composition with sequence-order correlation factors,
    as the amino acid composition alone carries no sequence-order information [4, 5].
    The correlation factor of rank k is:

    theta_k = (1/(L-k)) * sum_{i=1}^{L-k} (1/P) * sum_p (H_p(s_i) - H_p(s_i+k))^2

    where H_p are the P properties, each standardised to zero mean and unit
    standard deviation over the 20 amino acids, and L is the length of the
    sequence. The first 20 values are the weighted amino acid composition and the
    last lamda the weighted correlation factors theta_1..theta_lamda (see
    _pseudo_composition); the 20 + lamda values are percentages that sum to 100.
    By default the three properties of Chou (2001) are used: hydrophobicity,
    hydrophilicity and side-chain mass [6], with lamda=30 and weight=0.05, giving
    an output of 1 x 50.

    Parameters
    ==========
    :sequence: str
        protein sequence.
    :lamda: int (default=30)
        highest rank of correlation, less than the length of the sequence.
    :weight: float (default=0.05)
        weighting factor of the correlation factors.
    :properties: list[str | dict] | None (default=None)
        AAindex1 accession numbers and/or dicts of amino acid values; None uses
        hydrophobicity, hydrophilicity and side-chain mass.

    Returns
    =======
    :pseudo_amino_acid_composition_df: pd.Dataframe
        output dataframe of calculated pseudo amino acid composition descriptors 
        for input sequence. Dataframe will be of the shape 1 x N, where N is the 
        number of features calculated from the descriptor (20 + lamda) and 1 is 
        the input sequence. By default, the shape will be 1 x 50 (using default lamda=30).
    """
    #validate input protein sequence
    sequence = _validate_sequence(sequence)
    lamda = _check_lag("lamda", lamda, len(sequence), minimum=0)
    weight = _validate_weight(weight)
    props = _resolve_properties(properties, [hydrophobicity_, hydrophilicity_, residue_mass_])

    z = _property_matrix(sequence, props, standardize=True)
    theta = np.array([np.mean((z[:-k] - z[k:]) ** 2) for k in range(1, lamda + 1)])
    return _pseudo_composition(sequence, theta, weight, "PAAC")

def sequence_order_correlation_factor(sequence: str, k: int = 1, properties: list[dict[str, float]] | None = None) -> float:
    """
    Calculate the PseAAC sequence order correlation factor theta_k of a protein
    sequence for the given properties (see pseudo_amino_acid_composition).

    Deprecated since v1.6.0, removed in 2.0: protpy no longer uses it (PseAAC
    computes the factors directly) and it is not exported from the protpy
    namespace.

    Parameters
    ==========
    :sequence: str
        protein sequence.
    :k: int
        gap between amino acids in the sequence.
    :properties: list
        list of dicts of physicochemical/structural property values for amino acids.

    Returns
    =======
    :result: float
        correlation factor theta_k; 0 when there are no residue pairs.
    """
    _warn_deprecated("sequence_order_correlation_factor", "pseudo_amino_acid_composition")
    if not properties or len(sequence) <= k:
        return 0.0
    z = _property_matrix(sequence, list(properties), standardize=True)
    return float(np.mean((z[:-k] - z[k:]) ** 2))

def correlation_function(Ri: str = "S", Rj: str = "D", properties: list[dict[str, float]] | None = None) -> float:
    """
    Calculate the correlation between 2 input amino acids based on the 
    physicochemical/structural properties from the protein sequence.

    Deprecated since v1.6.0, removed in 2.0: it has no call sites in protpy.

    Parameters
    ==========
    :Ri: str
        1st amino acid.
    :Rj: str
        2nd amino acid.
    :properties: list | None (default=None)
        list of dicts of property values; None uses the PseAAC defaults
        (hydrophobicity, hydrophilicity and side-chain mass).

    Returns
    =======
    :correlation: float
        correlation value for the two input amino acids based on input properties.
    """
    _warn_deprecated("correlation_function")
    properties = properties or [hydrophobicity_, hydrophilicity_, residue_mass_]
    norm = [normalize_property(p) for p in properties]
    return sum((p[Ri] - p[Rj]) ** 2 for p in norm) / len(norm)

def normalize_property(properties: dict[str, float]) -> dict[str, float]:
    """
    Normalize physicochemical/structural property values using
    their mean and standard deviation. Used only by the deprecated
    correlation_function and amphiphilic_correllation_function; not exported
    from the protpy namespace and removed with them in 2.0.

    Parameters
    ==========
    :properties: dict
        dictionary of amino acid values for a physicochemical property.

    Returns
    =======
    :normalized_vals: dict
        dict of normalized property values.
    """
    vals = np.array(list(properties.values()), dtype=float)
    if not len(vals):
        return {}
    std = vals.std()
    return {aa: (v - vals.mean()) / std if std else 0.0 for aa, v in properties.items()}

######################## Amphiphilic Pseudo Amino Acid Composition ########################

def amphiphilic_pseudo_amino_acid_composition(sequence: str, lamda: int = 30, weight: float = 0.5,
    properties: list[str | dict[str, float]] | None = None) -> pd.DataFrame:
    """
    Calculate Amphiphilic Pseudo Amino Acid Composition (APAAComp) of protein
    sequence. APAAComp has the same form as the amino acid composition, but also
    describes the distribution of hydrophobic and hydrophilic amino acids along
    the chain [5]. For each property p and rank k the correlation factor is:

    tau_{k,p} = (1/(L-k)) * sum_{i=1}^{L-k} H_p(s_i) * H_p(s_i+k)

    where H_p is property p standardised to zero mean and unit standard deviation
    over the 20 amino acids. The first 20 values are the weighted amino acid
    composition and the next P * lamda the weighted correlation factors, ordered
    tau_{1,1}, tau_{1,2}, tau_{2,1}, ... (see _pseudo_composition); the values are
    percentages that sum to 100.

    The amphiphilic factors can be negative, so for short sequences with the
    default weight of 0.5 the normalisation term 1 + w * sum(tau) can approach 0
    and the values become unstable; a RuntimeWarning is raised when it is below
    0.1 and a ValueError when it is not positive. Use a smaller weight or lamda
    for short sequences.

    Parameters
    ==========
    :sequence: str
        protein sequence.
    :lamda: int (default=30)
        highest rank of correlation, less than the length of the sequence.
    :weight: float (default=0.5)
        weighting factor of the correlation factors.
    :properties: list[str | dict] | None (default=None)
        AAindex1 accession numbers and/or dicts of amino acid values; None uses
        hydrophobicity and hydrophilicity.

    Returns
    =======
    :amp_pseudo_amino_acid_composition_df: pd.Dataframe
        output dataframe of calculated amphiphilic pseudo amino acid composition
        descriptors for input sequence. Dataframe will be of the shape 1 x N, where
        N is 20 + n_props * lamda (n_props = number of properties supplied).
        With the default two properties and lamda=30 the shape is 1 x 80.
    """
    #validate input protein sequence
    sequence = _validate_sequence(sequence)
    lamda = _check_lag("lamda", lamda, len(sequence), minimum=0)
    weight = _validate_weight(weight)
    props = _resolve_properties(properties, [hydrophobicity_, hydrophilicity_])

    z = _property_matrix(sequence, props, standardize=True)
    #tau ordered [p0_lag1, p1_lag1, p0_lag2, ...]
    tau = np.array([np.mean(z[:-k] * z[k:], axis=0) for k in range(1, lamda + 1)]).reshape(-1)
    return _pseudo_composition(sequence, tau, weight, "APAAC")

def amphiphilic_sequence_order_correlation_factor(sequence: str, k: int = 1,
    properties: list[dict[str, float]] | None = None) -> list[float]:
    """
    Calculate the APseAAC sequence order correlation factors tau_{k,p} of a
    protein sequence, one per property (see
    amphiphilic_pseudo_amino_acid_composition).

    Deprecated since v1.6.0, removed in 2.0: protpy no longer uses it and it is
    not exported from the protpy namespace.

    Parameters
    ==========
    :sequence: str
        protein sequence.
    :k: int (default=1)
        gap between amino acids in the sequence.
    :properties: list[dict] | None (default=None → [``hydrophobicity_``, ``hydrophilicity_``])
        list of property dicts.

    Returns
    =======
    :correlation_factor: list[float]
        one correlation factor per property; 0 when there are no residue pairs.
    """
    _warn_deprecated("amphiphilic_sequence_order_correlation_factor", "amphiphilic_pseudo_amino_acid_composition")
    properties = properties or [hydrophobicity_, hydrophilicity_]
    if len(sequence) <= k:
        return [0.0] * len(properties)
    z = _property_matrix(sequence, list(properties), standardize=True)
    return [float(v) for v in np.mean(z[:-k] * z[k:], axis=0)]

def amphiphilic_correllation_function(Ri: str = "S", Rj: str = "D",
    normalized_props: list[dict[str, float]] | None = None) -> list[float]:
    """
    Calculate correlation value based on physicochemical property values
    for input amino acids Ri and Rj.

    Deprecated since v1.6.0, removed in 2.0: it has no call sites in protpy.

    Parameters
    ==========
    :Ri: str
        1st amino acid.
    :Rj: str
        2nd amino acid.
    :normalized_props: list[dict] | None
        pre-normalized property dicts. Defaults to normalized hydrophobicity
        and hydrophilicity when None.

    Returns
    =======
    :thetas: list[float]
        one correlation value per property.
    """
    _warn_deprecated("amphiphilic_correllation_function")
    if normalized_props is None:
        normalized_props = [normalize_property(hydrophobicity_),
                            normalize_property(hydrophilicity_)]
    return [prop[Ri] * prop[Rj] for prop in normalized_props]

####################### Composition of k-Spaced Amino Acid Pairs ######################

def cksaap(sequence: str, gap: int = 5, alphabet: int | str | list[str] | None = None) -> pd.DataFrame:
    """
    Calculate Composition of k-Spaced Amino Acid Pairs (CKSAAP) of protein
    sequence. CKSAAP describes the frequency of each amino acid pair separated
    by g residues within a protein sequence, for each gap g in the range
    0..gap. With 400 possible amino acid pairs per gap this creates a
    400 * (gap+1) Dimensional vector, calculated as:

    CKSAAP(XY, g) = N_g(XY) / (L - g - 1) * 100

    where CKSAAP(XY, g) is the composition of the amino acid pair XY with g
    residues between X and Y, N_g(XY) is the number of such pairs in the sequence
    and L is the length of the sequence [20]. A gap of 0 is the same as the
    dipeptide composition. A gap with no pairs (L - g - 1 <= 0) gives values
    of 0.

    With ``alphabet`` the amino acids are first mapped to the groups of a reduced
    alphabet (see ``kmer_composition``), giving G^2 group pairs per gap.
    ``alphabet="ifeature5"`` is the composition of k-spaced amino acid group pairs
    (CKSAAGP) of iFeature [21], also available as ``cksaagp()``.

    Parameters
    ==========
    :sequence: str
        protein sequence.
    :gap: int (default=5)
        maximum number of residues between the two amino acids of a pair.
    :alphabet: int | str | list[str] | None (default=None)
        reduced alphabet; None uses the 20 amino acids.

    Returns
    =======
    :cksaap_df: pd.DataFrame
        pandas dataframe of CKSAAP values for protein sequence. Dataframe will
        be of the shape 1 x 400 * (gap+1) (G^2 * (gap+1) with an alphabet), and 1
        is the input sequence. By default, the shape will be 1 x 2400. Column
        names are of the form CKSAAP_XY_gapg, e.g. CKSAAP_AC_gap2; with an alphabet
        CKSAAP_<alphabet>_G1.G2_gapg, or CKSAAGP_G1.G2_gapg for "ifeature5".
    """
    #validate input protein sequence
    sequence = _validate_sequence(sequence)
    gap = _validate_int("gap", gap)

    if alphabet is None:
        mapping, labels, fmt = None, amino_acids, "CKSAAP_{}{}_gap{}"
    else:
        prefix, groups = _alphabet(alphabet)
        mapping = {aa: label for label, residues in groups for aa in residues}
        labels = [label for label, _ in groups]
        fmt = ("CKSAAGP" if alphabet == "ifeature5" else "CKSAAP_" + prefix.replace("RA_", "")) + "_{}.{}_gap{}"
    pairs = list(itertools.product(labels, repeat=2))
    values, columns = [], []
    for g in range(gap + 1):
        counts = _gapped_pair_counts(sequence, g, mapping)
        total = len(sequence) - g - 1
        values += [counts[p] / total * 100 if total > 0 else 0.0 for p in pairs]
        columns += [fmt.format(x, y, g) for x, y in pairs]
    return _to_df(values, columns)

################# Composition of k-Spaced Amino Acid Group Pairs #######################

def cksaagp(sequence: str, gap: int = 5) -> pd.DataFrame:
    """
    Calculate Composition of k-Spaced Amino Acid Group Pairs (CKSAAGP) of protein
    sequence. CKSAAGP is the same as CKSAAP but the amino acids are first mapped
    to 5 groups: aliphatic (GAVLMI), aromatic (FYW), positive charge (KRH),
    negative charge (DE) and uncharged (STCPNQ). With 25 possible group pairs
    per gap this creates a 25 * (gap+1) Dimensional vector, calculated as:

    CKSAAGP(G1G2, g) = N_g(G1G2) / (L - g - 1) * 100

    where CKSAAGP(G1G2, g) is the composition of the group pair G1G2 with g
    residues between them, N_g(G1G2) is the number of such pairs in the sequence
    and L is the length of the sequence [21]. The 25 values for each gap sum to 100.
    This is a permanent alias of ``cksaap(sequence, gap, alphabet="ifeature5")``.

    Parameters
    ==========
    :sequence: str
        protein sequence.
    :gap: int (default=5)
        maximum number of residues between the two amino acids of a pair.

    Returns
    =======
    :cksaagp_df: pd.DataFrame
        pandas dataframe of CKSAAGP values for protein sequence. Dataframe will
        be of the shape 1 x 25 * (gap+1), where 25 * (gap+1) is the number of
        features calculated from the descriptor (25 group pairs for each gap) and
        1 is the input sequence. By default, the shape will be 1 x 150. Column
        names are of the form CKSAAGP_G1.G2_gapg, e.g. CKSAAGP_aliphatic.aromatic_gap0.
    """
    return cksaap(sequence, gap, alphabet="ifeature5")

############################ Dipeptide Deviation from Expected Mean ######################

def dde(sequence: str) -> pd.DataFrame:
    """
    Calculate Dipeptide Deviation from Expected Mean (DDE) of protein sequence.
    DDE compares the observed dipeptide composition of a protein sequence with
    the composition expected from the number of codons that code for each amino
    acid [22]. This is a permanent alias of
    ``dipeptide_composition(sequence, normalize="dde")``, which documents the
    formula. Unlike the other composition descriptors, DDE values are not
    percentages.

    Parameters
    ==========
    :sequence: str
        protein sequence, at least 2 residues long.

    Returns
    =======
    :dde_df: pd.DataFrame
        pandas dataframe of DDE values for protein sequence. Dataframe will be of
        the shape 1 x 400, where 400 is the number of features calculated from the
        descriptor (20^2 for the 20 canonical amino acids) and 1 is the input
        sequence. Column names are of the form DDE_XY.
    """
    return dipeptide_composition(sequence, normalize="dde")

################################ Adaptive Skip Dipeptide Composition #####################

def asdc(sequence: str) -> pd.DataFrame:
    """
    Calculate Adaptive Skip Dipeptide Composition (ASDC) of protein sequence.
    ASDC is a modified dipeptide composition that counts every ordered pair of
    amino acids (s_i, s_j), where i < j, at all gaps along the sequence, rather
    than just adjacent pairs, calculated as:

    ASDC(XY) = sum_g N_g(XY) / sum_g sum_XY N_g(XY) * 100

    where ASDC(XY) is the composition of the amino acid pair XY, N_g(XY) is the
    number of XY pairs separated by g residues, and the denominator is the total
    number of ordered pairs in the sequence, L(L-1)/2, where L is the length of
    the sequence [23]. The 400 values sum to 100.

    Parameters
    ==========
    :sequence: str
        protein sequence.

    Returns
    =======
    :asdc_df: pd.DataFrame
        pandas dataframe of ASDC values for protein sequence. Dataframe will be of
        the shape 1 x 400, where 400 is the number of features calculated from the
        descriptor (20^2 for the 20 canonical amino acids) and 1 is the input
        sequence. Column names are of the form ASDC_XY.
    """
    #validate input protein sequence
    sequence = _validate_sequence(sequence)

    index = {aa: i for i, aa in enumerate(amino_acids)}
    #integer matmul: exact counts, and avoids spurious BLAS float warnings
    onehot = np.zeros((len(sequence), 20), dtype=np.int64)
    onehot[np.arange(len(sequence)), [index[aa] for aa in sequence]] = 1
    #rows of `before` count each residue type strictly before position j
    before = np.cumsum(onehot, axis=0) - onehot
    pair_counts = before.T @ onehot
    total = len(sequence) * (len(sequence) - 1) / 2
    values = (pair_counts / total * 100).ravel() if total else np.zeros(400)
    return _to_df(list(values), [f"ASDC_{x}{y}" for x, y in itertools.product(amino_acids, repeat=2)])

#################################### Distance Pair #######################################

def distance_pair(sequence: str, rng: int = 3, alphabet: str = "cp14") -> pd.DataFrame:
    """
    Calculate Distance Pair (DP) descriptor of protein sequence. The descriptor
    incorporates amino acid distance pairs and a reduced amino acid alphabet into
    the general pseudo amino acid composition, as used in iDNA-Prot|dis [24, 25].
    The amino acids are first mapped to the clusters of a cluster profile cp(n),
    and the descriptor is the composition of the clusters followed by the
    composition of cluster pairs at each distance d = 1..rng:

    DP(c) = N(c) / L * 100

    DP(c1c2, d) = N_d(c1c2) / (L - d) * 100

    where N(c) is the number of residues in cluster c, N_d(c1c2) is the number of
    residue pairs (s_i, s_i+d) whose residues are in clusters c1 and c2, and L is
    the length of the sequence. A distance of 1 is the nearest-neighbour pair.
    The available cluster profiles are cp13, cp14, cp19 and cp20 (no reduction).
    The published cp(14) lists E twice; the extra E in the first cluster is
    dropped, as it is the only reading that gives 14 clusters.

    Parameters
    ==========
    :sequence: str
        protein sequence.
    :rng: int (default=3)
        maximum distance between the two residues of a pair.
    :alphabet: str (default="cp14")
        cluster profile to use: "cp13", "cp14", "cp19" or "cp20".

    Returns
    =======
    :distance_pair_df: pd.DataFrame
        pandas dataframe of DP values for protein sequence. Dataframe will be of
        the shape 1 x (G + G^2 * rng), where G is the number of clusters,
        G + G^2 * rng is the number of features calculated from the descriptor and
        1 is the input sequence. By default (cp14, rng=3) the shape will be 1 x 602.
        Column names are of the form DP_c and DP_c1.c2_dd, where c are the cluster
        residues, e.g. DP_IMV and DP_IMV.L_d1.
    """
    #validate input protein sequence
    sequence = _validate_sequence(sequence)
    rng = _validate_int("rng", rng)
    if alphabet not in ("cp13", "cp14", "cp19", "cp20"):
        raise ValueError(f"Invalid alphabet '{alphabet}'. Valid alphabets are: ['cp13', 'cp14', 'cp19', 'cp20'].")

    labels = [label for label, _ in _alphabet(alphabet)[1]]
    mapping = {aa: label for label in labels for aa in label}
    comp = Counter(mapping[aa] for aa in sequence)
    values = [comp[c] / len(sequence) * 100 for c in labels]
    columns = [f"DP_{c}" for c in labels]
    pairs = list(itertools.product(labels, repeat=2))
    for d in range(1, rng + 1):
        counts = _gapped_pair_counts(sequence, d - 1, mapping)
        total = len(sequence) - d
        values += [counts[p] / total * 100 if total > 0 else 0.0 for p in pairs]
        columns += [f"DP_{a}.{b}_d{d}" for a, b in pairs]
    return _to_df(values, columns)

#################################### Natural Vector ######################################

def natural_vector(sequence: str, normalize: bool = False) -> pd.DataFrame:
    """
    Calculate Natural Vector (NV) of protein sequence. The natural vector
    describes the number of each amino acid in a protein sequence together with
    the mean and the normalized central moment of their positions, calculated as:

    n_k = number of occurrences of amino acid k

    mu_k = sum of the positions of k / n_k

    D2_k = sum of (position of k - mu_k)^2 / (n_k * L)

    where positions start at 1, L is the length of the sequence and k is each of
    the 20 canonical amino acids [26]. mu_k and D2_k are 0 when amino acid k is
    not in the sequence. The n_k block is the amino acid composition times L/100;
    with ``normalize=True`` n_k and mu_k are divided by L (the composition as a
    fraction, and the relative mean position), so the vector no longer scales
    with sequence length.

    Parameters
    ==========
    :sequence: str
        protein sequence.
    :normalize: bool (default=False)
        divide n_k and mu_k by the length of the sequence.

    Returns
    =======
    :natural_vector_df: pd.DataFrame
        pandas dataframe of NV values for protein sequence. Dataframe will be of
        the shape 1 x 60, where 60 is the number of features calculated from the
        descriptor (3 values for each of the 20 amino acids) and 1 is the input
        sequence. Column names are of the form NV_n_X, NV_mu_X and NV_D2_X.
    """
    #validate input protein sequence
    sequence = _validate_sequence(sequence)

    L = len(sequence)
    positions = {aa: [] for aa in amino_acids}
    for i, aa in enumerate(sequence, start=1):
        positions[aa].append(i)
    n = [len(positions[aa]) for aa in amino_acids]
    mu = [sum(positions[aa]) / len(positions[aa]) if positions[aa] else 0.0 for aa in amino_acids]
    d2 = [sum((p - m) ** 2 for p in positions[aa]) / (len(positions[aa]) * L) if positions[aa] else 0.0
          for aa, m in zip(amino_acids, mu)]
    if normalize:
        n, mu = [v / L for v in n], [v / L for v in mu]
    columns = ([f"NV_n_{aa}" for aa in amino_acids] + [f"NV_mu_{aa}" for aa in amino_acids]
               + [f"NV_D2_{aa}" for aa in amino_acids])
    return _to_df(n + mu + d2, columns)

#################################### Atomic Composition ##################################

def atomic_composition(sequence: str, subtract_water: bool = False) -> pd.DataFrame:
    """
    Calculate Atomic Composition (ATC) of protein sequence. ATC describes the
    percentage of carbon, hydrogen, nitrogen, oxygen and sulphur atoms in a
    protein sequence, using the molecular formulas of the free amino acids,
    calculated as:

    ATC(e) = N_e / (N_C + N_H + N_N + N_O + N_S) * 100

    where ATC(e) is the atomic composition of element e and N_e is the number of
    atoms of element e summed over all residues [27]. By default the formulas of
    the free amino acids are summed without removing water, matching Pfeature
    [27]. With subtract_water=True, one water molecule is removed for each of the
    L - 1 peptide bonds, giving the formula of the peptide (e.g. "GG" -> C4H8N2O3).

    Parameters
    ==========
    :sequence: str
        protein sequence.
    :subtract_water: bool (default=False)
        remove one H2O per peptide bond.

    Returns
    =======
    :atomic_composition_df: pd.DataFrame
        pandas dataframe of ATC values for protein sequence. Dataframe will be of
        the shape 1 x 5, where 5 is the number of features calculated from the
        descriptor (one per element) and 1 is the input sequence. Column names
        are ATC_C, ATC_H, ATC_N, ATC_O and ATC_S.
    """
    #validate input protein sequence
    sequence = _validate_sequence(sequence)

    formulas = _load_json_table("residue_tables.json", "atomic_formulas")["values"]
    atoms = {e: sum(formulas[aa][e] for aa in sequence) for e in "CHNOS"}
    if subtract_water:
        atoms["H"] -= 2 * (len(sequence) - 1)
        atoms["O"] -= len(sequence) - 1
    total = sum(atoms.values())
    return _to_df([atoms[e] / total * 100 for e in "CHNOS"], [f"ATC_{e}" for e in "CHNOS"])
