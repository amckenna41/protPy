################################################################################
#################                    CTD                       #################
################################################################################

from __future__ import annotations

import pandas as pd
import math
from difflib import get_close_matches

#list of amino acids
from ._constants import _validate_sequence
from ._helpers import _load_json_table, _to_df

#str_to_num is an internal encoder and is not exported from the protpy namespace
__all__ = ["ctd", "ctd_", "ctd_composition", "ctd_transition", "ctd_distribution"]

"""
References
==========
[1] Inna Dubchak, Ilya Muchink, Stephen R.Holbrook and Sung-Hou Kim.
    Prediction of protein folding class using global description of amino
    acid sequence. Proc.Natl. Acad.Sci.USA, 1995, 92, 8700-8704.
[2] Inna Dubchak, Ilya Muchink, Christopher Mayor, Igor Dralyuk and Sung-Hou
    Kim. Recognition of a Protein Fold in the Context of the SCOP
    classification. Proteins: Structure, Function and Genetics, 1999, 35, 401-407.
[3] Gromiha, M. M. (2010). Protein Sequence Analysis. In M. M. Gromiha (Ed.), 
    Protein Bioinformatics (pp. 29–62). Elsevier.
[4] Chen, Z. et al. (2018). iFeature: a Python package and web server for features
    extraction and selection from protein and peptide sequences. Bioinformatics,
    34(14), 2499-2502. https://doi.org/10.1093/bioinformatics/bty140
"""

hydrophobicity = {"name": "hydrophobicity", "1": "RKEDQN", "2": "GASTPHY", "3": "CLVIMFW"}
# '1' -> Polar; '2' -> Neutral, '3' -> Hydrophobicity

normalized_vdwv = {"name": "normalized_vdwv", "1": "GASTPDC", "2": "NVEQIL", "3": "MHKFRYW"}
# '1' -> (0-2.78); '2' -> (2.95-4.0), '3' -> (4.03-8.08)

polarity = {"name": "polarity", "1": "LIFWCMVY", "2": "PAGST", "3": "HQRKNED"}
# '1' -> (4.9-6.2); '2' -> (8.0-9.2), '3' -> (10.4-13.0)

charge = {"name": "charge", "1": "KR", "2": "ANCQGHILMFPSTWYV", "3": "DE"}
# '1' -> Positive; '2' -> Neutral, '3' -> Negative

secondary_struct = {"name": "secondary_struct", "1": "EALMQKRH", "2": "VIYCWFT", "3": "GNPSD"}
# '1' -> Helix; '2' -> Strand, '3' -> coil

solvent_accessibility = {"name": "solvent_accessibility", "1": "ALFCGIVW", "2": "RKQEND", "3": "MPSTHY"}
# '1' -> Buried; '2' -> Exposed, '3' -> Intermediate

polarizability = {"name": "polarizability", "1": "GASDT", "2": "CPNVEQIL", "3": "KMHFRYW"}
# '1' -> (0-0.108); '2' -> (0.128-0.186), '3' -> (0.219-0.409)

#object of physicochemical properties to use for calculating CTD descriptors
ctd_properties = {
    "hydrophobicity": hydrophobicity,
    "normalized_vdwv": normalized_vdwv,
    "polarity": polarity,
    "charge": charge,
    "secondary_struct": secondary_struct,
    "solvent_accessibility": solvent_accessibility,
    "polarizability": polarizability
    }

#six extra hydrophobicity partitions used by iFeature's 13-property CTD [4]; kept out of
#ctd_properties so the default 7-property outputs are unchanged
ctd_properties_ifeature_extra = {
    name: {"name": name, **groups}
    for name, groups in _load_json_table("reduced_alphabets.json", "ctd_ifeature_extra")["values"].items()
    }

#property sets accepted by ctd_(property_set=...)
ctd_property_sets = {
    "default": ctd_properties,
    "ifeature13": {**ctd_properties, **ctd_properties_ifeature_extra},
    }

#every property accepted by the single-property CTD functions
_all_ctd_properties = ctd_property_sets["ifeature13"]

def str_to_num(sequence: str, property: dict[str, str]) -> str:
    """
    Convert sequences str to number from input physicochemical property.

    Parameters
    ==========
    :sequence: str
        protein sequence.
    :property: str
        physicochemical property name to use when calculating descriptor.

    Returns
    =======
    :sequence_converted: str
        converted protein sequence into numerical format. 
    """
    #validate input protein sequence
    sequence = _validate_sequence(sequence)

    sequence_converted = sequence

    #convert amino acid into corresponding CTD keys
    for key, value in list(property.items()):
        if (key == "name"):
            continue
        for index in value:
            sequence_converted = sequence_converted.replace(index, key)

    return sequence_converted

def _ctd_property(property: str) -> dict[str, str]:
    """
    Get a CTD property by name, accepting close misspellings (difflib cutoff 0.8).

    Parameters
    ==========
    :property: str
        CTD property name, e.g. "hydrophobicity".

    Returns
    =======
    :prop: dict
        the property's name and its three residue classes.
    """
    matches = get_close_matches(property, _all_ctd_properties.keys(), cutoff=0.8)
    if not matches:
        raise ValueError(f"Invalid property '{property}'. Available properties are: {list(_all_ctd_properties.keys())}.")
    return _all_ctd_properties[matches[0]]

##################################### CTD Composition #####################################

def ctd_composition(sequence: str, property: str = "hydrophobicity") -> pd.DataFrame:
    """
    Calculate composition physicochemical/structural descriptor. Composition 
    is determined as the number of amino acids of a particular property 
    divided by the total number of amino acids. The shape of the output will be 
    1 x 3, with 3 features being generated per sequence.

    Parameters
    ==========
    :sequence: str
        protein sequence.
    :property: str (default="hydrophocity")
        physicochemical property name to use when calculating descriptor.

    Returns
    =======
    :ctd_composition_df: pd.DataFrame
        dataframe of calculated composition values for sequence using
        selected physicochemical property. Output will be of shape 1 x 3,
        with 3 features being generated per sequence.
    """
    #validate input protein sequence
    sequence = _validate_sequence(sequence)

    prop = _ctd_property(property)
    
    #convert sequence to number
    seq = str_to_num(sequence, prop)

    ctd_composition_ = {}
    
    #calculate descriptor values, append to ctd_composition_ dict
    for key in ("1", "2", "3"):
        ctd_composition_[f"CTD_C_0{key}_" + prop["name"]] = seq.count(key) / len(sequence)

    return _to_df(list(ctd_composition_.values()), list(ctd_composition_))

##################################### CTD Transition ######################################

def ctd_transition(sequence: str, property: str = "hydrophobicity") -> pd.DataFrame:
    """
    Calculate transition physicochemical/structural descriptor. Transition 
    is determined as the number of transitions from a particular property 
    to different property divided by (total number of amino acids − 1). 
    The shape of the output will be 1 x 3, with 3 features being generated 
    per sequence.

    Parameters
    ==========
    :sequence: str
        protein sequence.
    :property: str (default="hydrophocity")
        physicochemical property name to use when calculating descriptor.

    Returns
    =======
    :ctd_transition_df: pd.DataFrame
        dataframe of calculated transition values for sequence using
        selected physicochemical property. Output will be of shape 1 x 3,
        with 3 features being generated per sequence.
    """
    #validate input protein sequence
    sequence = _validate_sequence(sequence)

    prop = _ctd_property(property)

    #convert sequence to number
    seq = str_to_num(sequence, prop)
    
    ctd_transition_ = {}

    #a single residue has no transitions
    n_pairs = len(sequence) - 1
    if n_pairs == 0:
        n_pairs = 1

    #calculate descriptor values, append to ctd_transition_ dict
    for a, b in (("1", "2"), ("1", "3"), ("2", "3")):
        ctd_transition_[f"CTD_T_{a}{b}_" + prop["name"]] = (seq.count(a + b) + seq.count(b + a)) / n_pairs

    return _to_df(list(ctd_transition_.values()), list(ctd_transition_))

#################################### CTD Distribution #####################################

def ctd_distribution(sequence: str, property: str = "hydrophobicity") -> pd.DataFrame:
    """
    Calculate distribution physicochemical/structural descriptor. Distribution 
    is the chain length within which the first, 25%, 50%, 75% and 100% of the 
    amino acids of a particular property are located. The shape of the output 
    will be 1 x 15, with 15 features being generated per sequence.

    The position of the q-th quantile residue is the floor(q * n)-th residue of
    the class, where n is the number of residues in the class, clamped to at least
    the first residue as in iFeature [4]. (Before v1.6.0 a class with 1-3 residues
    used the last residue for the 25% quantile.)

    Parameters
    ==========
    :sequence: str
        protein sequence.
    :property: str (default="hydrophocity")
        physicochemical property name to use when calculating descriptor.

    Returns
    =======
    :ctd_distribution_df: pd.DataFrame
        dataframe of calculated distribution values for sequence using
        selected physicochemical property. Output will be of shape 1 x 15,
        with 15 features being generated per sequence.
    """
    #validate input protein sequence
    sequence = _validate_sequence(sequence)

    prop = _ctd_property(property)

    #convert sequence to number
    seq = str_to_num(sequence, prop)
    
    ctd_distribution_ = {}
    
    #iterate through each class, calculating the positions of its first, 25%, 50%, 75% and last residues
    for key in ("1", "2", "3"):
        positions = [i + 1 for i, c in enumerate(seq) if c == key]
        num = len(positions)
        for label, q in (("001", None), ("025", 0.25), ("050", 0.5), ("075", 0.75), ("100", 1.0)):
            if not positions:
                value = 0.0
            else:
                #floor(q * n)-th residue of the class, clamped to the first residue (iFeature CTDD)
                rank = 1 if q is None else max(1, math.floor(num * q))
                value = positions[rank - 1] / len(seq) * 100
            ctd_distribution_["CTD_D_0" + key + "_" + label + "_" + prop["name"]] = value

    return _to_df(list(ctd_distribution_.values()), list(ctd_distribution_))

def ctd(sequence: str, property: str = "hydrophobicity", all_ctd: bool = True,
    property_set: str = "default", parts: str = "CTD") -> pd.DataFrame:
    """
    Calculate all Composition, Transition and Distribution (CTD) features of protein sequences.
    Composition is the number of amino acids of a particular property (e.g., hydrophobicity)
    divided by the total number of amino acids in a protein sequence. Transition
    characterizes the percent frequency with which amino acids of a particular
    property is followed by amino acids of a different property. Distribution
    measures the chain length within which the first, 25%, 50%, 75%, and 100% of
    the amino acids of a particular property are located, respectively [1, 2].

    CTD properties available are: Polarizability, Solvent Accessibility, Secondary 
    Structure, Charge, Polarity, Normalized VDWV, Hydrophobicity. Each property generates an 
    output of shape 1 x 21, 3/21 will be Composition, 3/21 will be Transition, 15/21 will be
    Distribution. When Calculating all available features the generated output will be of shape 
    1 x 147, 21/147 will be composition, 21/147 will be transition and the remaining 105 are 
    distribution.
    
    Parameters
    ==========
    :sequence: str
        protein sequence.
    :property: str (default="hydrophocity")
        physicochemical property name to use when calculating descriptor.
    :all_ctd: bool (default=True)
        calculate all CTD descriptors and concatenate together.
    :property_set: str (default="default")
        property set used when ``all_ctd=True``: "default" (the 7 properties
        above, 1 x 147) or "ifeature13" (adds the six AAindex hydrophobicity
        partitions used by iFeature, 13 properties, 1 x 273) [4].
    :parts: str (default="CTD")
        which parts to calculate, any combination of "C" (composition),
        "T" (transition) and "D" (distribution), e.g. "CT". The single parts are
        also available as ``ctd_composition``, ``ctd_transition`` and
        ``ctd_distribution``, and ``ctd_`` is an alias of this function.

    Returns
    =======
    :ctd_df: pd.DataFrame
        dataframe of CTD descriptor values for all protein sequences. DataFrame will
        be of the shape 1 x 147, where 147 is the total number of features calculated from
        the CTD descriptors per sequence, with each property generating an output of 1 x 21.
    """
    #validate input protein sequence
    sequence = _validate_sequence(sequence)
    if property_set not in ctd_property_sets:
        raise ValueError(f"Invalid property_set '{property_set}'. Available sets are: {list(ctd_property_sets)}.")
    if not parts or set(parts.upper()) - set("CTD"):
        raise ValueError(f"Invalid parts '{parts}'. Use a combination of 'C', 'T' and 'D', e.g. 'CTD'.")

    props = list(ctd_property_sets[property_set]) if all_ctd else [property]
    funcs = {"C": ctd_composition, "T": ctd_transition, "D": ctd_distribution}
    #all composition columns, then all transition columns, then all distribution columns
    frames = [funcs[part](sequence, property=prop) for part in "CTD" if part in parts.upper() for prop in props]
    return pd.concat(frames, axis=1)

def ctd_(sequence: str, property: str = "hydrophobicity", all_ctd: bool = True,
    property_set: str = "default") -> pd.DataFrame:
    """
    Calculate all Composition, Transition and Distribution (CTD) features of
    protein sequences. Alias of ``ctd(sequence, property, all_ctd, property_set)``;
    see ``ctd`` for the description and parameters.

    Returns
    =======
    :ctd_df: pd.DataFrame
        dataframe of CTD descriptor values, by default of the shape 1 x 147.
    """
    return ctd(sequence, property, all_ctd, property_set)