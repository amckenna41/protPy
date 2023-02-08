################################################################################
#################                    CTD                       #################
################################################################################

import pandas as pd
import math
import copy
from varname import nameof
from difflib import get_close_matches

#list of amino acids
amino_acids = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", 
    "Q", "R", "S", "T", "V", "W", "Y"]

"""
References
----------
[1] Inna Dubchak, Ilya Muchink, Stephen R.Holbrook and Sung-Hou Kim.
    Prediction of protein folding class using global description of amino
    acid sequence. Proc.Natl. Acad.Sci.USA, 1995, 92, 8700-8704.

[2] Inna Dubchak, Ilya Muchink, Christopher Mayor, Igor Dralyuk and Sung-Hou
    Kim. Recognition of a Protein Fold in the Context of the SCOP
    classification. Proteins: Structure, Function and
    Genetics, 1999, 35, 401-407.
"""

hydrophobicity = {"name": "hydrophobicity", "1": "RKEDQN", "2": "GASTPHY", "3": "CLVIMFW"}
# '1' -> Polar; '2' -> Neutral, '3' -> Hydrophobicity

normalized_vdwv = {"name": "normalized_vdwv", "1": "GASTPD", "2": "NVEQIL", "3": "MHKFRYW"}
# '1' -> (0-2.78); '2' -> (2.95-4.0), '3' -> (4.03-8.08)

polarity = {"name": "polarity", "1": "LIFWCMVY", "2": "CPNVEQIL", "3": "KMHFRYW"}
# '1' -> (4.9-6.2); '2' -> (8.0-9.2), '3' -> (10.4-13.0)

charge = {"name": "charge", "1": "KR", "2": "ANCQGHILMFPSTWYV", "3": "DE"}
# '1' -> Positive; '2' -> Neutral, '3' -> Negative

secondary_struct = {"name": "secondary_struct", "1": "EALMQKRH", "2": "VIYCWFT", "3": "GNPSD"}
# '1' -> Helix; '2' -> Strand, '3' -> coil

solvent_accessibility = {"name": "solvent_accessibility", "1": "ALFCGIVW", "2": "RKQEND", "3": "MPSTHY"}
# '1' -> Buried; '2' -> Exposed, '3' -> Intermediate

polarizability = {"name": "polarizability", "1": "GASDT", "2": "CPNVEQIL", "3": "KMHFRYW"}
# '1' -> (0-0.108); '2' -> (0.128-0.186), '3' -> (0.219-0.409)

#object of physiochemical properties to use for calculating CTD descriptors 
ctd_properties = {
    nameof(hydrophobicity): hydrophobicity,
    nameof(normalized_vdwv): normalized_vdwv,
    nameof(polarity): polarity,
    nameof(charge): charge,
    nameof(secondary_struct): secondary_struct,
    nameof(solvent_accessibility): solvent_accessibility,
    nameof(polarizability): polarizability
}

def str_to_num(sequence, property):
    """
    Convert sequences str to number from input physiochemical property.

    Parameters
    ----------
    :sequence : str
        protein sequence.
    :property : str
        physiochemical property name to use when calculating descriptor.

    Returns
    -------
    :sequence_converted : str
        converted protein sequence into numerical format. 
    """
    #check input sequence is a string, if not raise type error
    if not isinstance(sequence, str):
        raise TypeError('Input sequence must be a string, got input of type {}'.format(type(sequence)))

    #if invalid amino acids in sequence, raise value error
    for aa in sequence:
        if (aa not in amino_acids):
            raise ValueError("Invalid amino acid in protein sequence: ".format(aa))

    sequence_converted = copy.deepcopy(sequence)

    #convert amino acid into corresponding CTD keys
    for key, value in list(property.items()):
        if (key == "name"):
            continue
        for index in value:
            sequence_converted = sequence_converted.replace(index, key)

    return sequence_converted

def ctd_composition(sequence, property="hydrophobicity"):
    """
    Calculate composition physiochemical/structural descriptor.

    Parameters
    ----------
    :sequence : str
        protein sequence.
    :property : str (default="hydrophocity")
        physiochemical property name to use when calculating descriptor.

    Returns
    -------
    :ctd_composition_df : pd.DataFrame
        dataframe of calculated composition values for sequence using
        selected physiochemical property.
    """
    #check input sequence is a string, if not raise type error
    if not isinstance(sequence, str):
        raise TypeError('Input sequence must be a string, got input of type {}.'.format(type(sequence)))

    #uppercase protein sequence 
    sequence = sequence.upper()

    #if invalid amino acids in sequence, raise value error
    for aa in sequence:
        if (aa not in amino_acids):
            raise ValueError("Invalid amino acid in protein sequence: .".format(aa))

    #get closest matched CTD property 
    property_matches = get_close_matches(property, ctd_properties.keys(), cutoff=0.8)

    #get CTD property values from dict, if no similar match found set to default of hydrophobicity
    if (property_matches != []):
        prop = ctd_properties[property_matches[0]]
    else:
        prop = ctd_properties["hydrophobicity"]
    
    #convert sequence to number
    seq = str_to_num(sequence, prop)

    ctd_composition_ = {}
    
    #calculate descriptor values, append to ctd_composition_ dict
    ctd_composition_[prop["name"] + '_CTD_C_01'] = round(float(seq.count("1"))/len(sequence), 3) 
    ctd_composition_[prop["name"] + '_CTD_C_02'] = round(float(seq.count("2"))/len(sequence), 3)
    ctd_composition_[prop["name"] + '_CTD_C_03'] = round(float(seq.count("3"))/len(sequence), 3)

    #transform values and columns to DataFrame
    ctd_composition_df = pd.DataFrame([list(ctd_composition_.values())], columns=list(ctd_composition_.keys()))

    return ctd_composition_df

def ctd_transition(sequence, property="hydrophobicity"):
    """
    Calculate transition physiochemical/structural descriptor.

    Parameters
    ----------
    :sequence : str
        protein sequence.
    :property : str (default="hydrophocity")
        physiochemical property name to use when calculating descriptor.

    Returns
    -------
    :ctd_transition_df : pd.DataFrame
        dataframe of calculated transition values for sequence using
        selected physiochemical property.
    """
    #check input sequence is a string, if not raise type error
    if not isinstance(sequence, str):
        raise TypeError('Input sequence must be a string, got input of type {}.'.format(type(sequence)))

    #uppercase protein sequence 
    sequence = sequence.upper()

    #if invalid amino acids in sequence, raise value error
    for aa in sequence:
        if (aa not in amino_acids):
            raise ValueError("Invalid amino acid in protein sequence: .".format(aa))

    #get closest matched CTD property 
    property_matches = get_close_matches(property, ctd_properties.keys(), cutoff=0.8)

    #get CTD property values from dict, if no similar match found set to default of hydrophobicity
    if (property_matches != []):
        prop = ctd_properties[property_matches[0]]
    else:
        prop = ctd_properties["hydrophobicity"]

    #convert sequence to number
    seq = str_to_num(sequence, prop)
    
    ctd_transition_ = {}

    #calculate descriptor values, append to ctd_transition_ dict
    ctd_transition_[prop["name"] + "_CTD_T_12"] = round(
        float(seq.count("12") + seq.count("21")) / (len(sequence)-1), 3)
    ctd_transition_[prop["name"] + "_CTD_T_13"] = round(
        float(seq.count("13") + seq.count("31")) / (len(sequence)-1), 3)
    ctd_transition_[prop["name"] + "_CTD_T_23"] = round(
        float(seq.count("23") + seq.count("32")) / (len(sequence)-1), 3)

    #transform values and columns to DataFrame
    ctd_transition_df = pd.DataFrame([list(ctd_transition_.values())], columns=list(ctd_transition_.keys()))

    return ctd_transition_df

def ctd_distribution(sequence, property="hydrophobicity"):
    """
    Calculate distribution physiochemical/structural descriptor.

    Parameters
    ----------
    :sequence : str
        protein sequence.
    :property : str (default="hydrophocity")
        physiochemical property name to use when calculating descriptor.

    Returns
    -------
    :ctd_distribution_df : pd.DataFrame
        dataframe of calculated distribution values for sequence using
        selected physiochemical property.
    """
    #check input sequence is a string, if not raise type error
    if not isinstance(sequence, str):
        raise TypeError('Input sequence must be a string, got input of type {}.'.format(type(sequence)))

    #if invalid amino acids in sequence, raise value error
    for aa in sequence:
        if (aa not in amino_acids):
            raise ValueError("Invalid amino acid in protein sequence: .".format(aa))

    #get closest matched CTD property 
    property_matches = get_close_matches(property, ctd_properties.keys(), cutoff=0.8)

    #get CTD property values from dict, if no similar match found set to default of hydrophobicity
    if (property_matches != []):
        prop = ctd_properties[property_matches[0]]
    else:
        prop = ctd_properties["hydrophobicity"]

    #convert sequence to number
    seq = str_to_num(sequence, prop)
    
    ctd_distribution_ = {}
    
    #iterate through sequence, calculating distribution descriptor values using property
    for key, value in prop.items():
       if (key=="name"):
        continue
       num = seq.count(key)
       ink = 1
       indexk = 0
       cds = []
       while ink <= num:
           indexk = seq.find(key,indexk) + 1
           cds.append(indexk)
           ink = ink + 1
       if cds == []:
           ctd_distribution_[prop["name"] + "_CTD_D_0" + key + "_001"] = 0
           ctd_distribution_[prop["name"] + "_CTD_D_0" + key + "_025"] = 0
           ctd_distribution_[prop["name"] + "_CTD_D_0" + key + "_050"] = 0
           ctd_distribution_[prop["name"] + "_CTD_D_0" + key + "_075"] = 0
           ctd_distribution_[prop["name"] + "_CTD_D_0" + key + "_100"] = 0
       else:
           ctd_distribution_[prop["name"] + "_CTD_D_0" + key + "_001"] = round(
            float(cds[0]) / len(seq) * 100, 3
           )
           ctd_distribution_[prop["name"] + "_CTD_D_0" + key + "_025"] = round(
            float(cds[int(math.floor(num * 0.25)) - 1]) / len(seq) * 100, 3
           )
           ctd_distribution_[prop["name"] + "_CTD_D_0" + key + "_050"] = round(
            float(cds[int(math.floor(num * 0.5)) - 1]) / len(seq) * 100, 3
           )
           ctd_distribution_[prop["name"] + "_CTD_D_0" + key + "_075"] = round(
            float(cds[int(math.floor(num * 0.75)) - 1]) / len(seq) * 100, 3
           )
           ctd_distribution_[prop["name"] + "_CTD_D_0" + key + "_100"] = round(
            float(cds[-1]) / len(seq) * 100, 3)

    #transform values and columns to DataFrame
    ctd_distribution_df = pd.DataFrame([list(ctd_distribution_.values())], columns=list(ctd_distribution_.keys()))

    return ctd_distribution_df

def ctd_(sequence, property="hydrophobicity", all_ctd=True):
    """
    Calculate Composition, transition and distribution (CTD) features of protein sequences.
    Composition is the number of amino acids of a particular property (e.g., hydrophobicity)
    divided by the total number of amino acids in a protein sequence. Transition
    characterizes the percent frequency with which amino acids of a particular
    property is followed by amino acids of a different property. Distribution
    measures the chain length within which the first, 25%, 50%, 75%, and 100% of
    the amino acids of a particular property are located, respectively [6].
    CTD properties used are: Polarizability, Solvent Accessibility, Secondary 
    Structure, Charge, Polarity, Normalized VDWV, Hydrophobicity. The output 
    will be of shape 1 x 147. 21/147 will be
    composition, 21/147 will be transition and the remaining 105 are distribution.
    
    Parameters
    ----------
    :sequence : str
        protein sequence.
    :property : str (default="hydrophocity")
        physiochemical property name to use when calculating descriptor.
    :ctd : bool
        calculate all CTD descriptors and concatenate together.

    Returns
    -------
    :ctd_df : pd.DataFrame
        dataframe of CTD descriptor values for all protein sequences. DataFrame will
        be of the shape 1 x 147, where 147 is the number of features calculated from t
        he descriptors.
    """
    #check input sequence is a string, if not raise type error
    if not isinstance(sequence, str):
        raise TypeError('Input sequence must be a string, got input of type {}.'.format(type(sequence)))

    #if invalid amino acids in sequence, raise value error
    for aa in sequence:
        if (aa not in amino_acids):
            raise ValueError("Invalid amino acid in protein sequence: .".format(aa))
    
    #initialise ctd dataframes
    comp_df = pd.DataFrame()
    trans_df = pd.DataFrame()
    distr_df = pd.DataFrame()

    #if using single property, calculate each of the CTD descriptors individually
    if not all_ctd:
        comp_df = ctd_composition(sequence, property=property)
        trans_df = ctd_transition(sequence, property=property)
        distr_df = ctd_distribution(sequence, property=property)
    else:
        #if using all calculable properties, calculate CTD descriptors for each property
        for prop in ctd_properties:
            comp = ctd_composition(sequence, property=prop)
            comp_df = pd.concat([comp_df, comp], axis=1)
            trans = ctd_transition(sequence, property=prop)
            trans_df = pd.concat([trans_df, trans], axis=1)
            distr = ctd_distribution(sequence, property=prop)
            distr_df = pd.concat([distr_df, distr], axis=1)

    #concatenate all descriptors
    ctd = pd.concat([comp_df, trans_df, distr_df], axis=1)
    
    return ctd