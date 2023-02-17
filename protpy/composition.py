################################################################################
#############                   Composition                       ##############
################################################################################

import re
import pandas as pd
import math
from aaindex import aaindex1

"""
References
----------
[1] Reczko, M. and Bohr, H. (1994) The DEF data base of sequence based protein
    fold class predictions. Nucleic Acids Res, 22, 3616-3619.

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
"""
#list of amino acids
amino_acids = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", 
    "Q", "R", "S", "T", "V", "W", "Y"]

#physiochemical values for all amino acids for hydrophobicity, hydrophilicty and residue mass
#values taken from http://www.csbio.sjtu.edu.cn/bioinf/PseAAC/ParaValue.htm
#property values differ from their equivalent in the aaindex so hard-coding them below
hydrophobicity_ = {
        "A": 0.62, "C": 0.29, "D": -0.90, "E": -0.74, "F": 1.19, "G": 0.48, 
        "H": -0.40, "I": 1.38, "K": -1.50, "L": 1.06, "M": 0.64, "N": -0.78,
        "P": 0.12, "Q": -0.85, "R": -2.53, "S": -0.18, "T": -0.05, "V": 1.08, 
        "W": 0.81, "Y": 0.26
        }

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

def amino_acid_composition(sequence):
    """
    Calculate Amino Acid Composition (AAComp) of protein sequence. AAComp
    describes the fraction of each amino acid type within a protein sequence,
    and is calculated as:

    AA_Comp(s) = AA(t)/N(s)

    where AA_Comp(s) is the AAComp of protein sequence s, AA(t) is the number
    of amino acid types t (where t = 1,2,..,20) and N(s) is the length of the
    sequence s.

    Parameters
    ----------
    :sequence : str
        protein sequence.

    Returns
    -------
    :amino_acid_composition_df : pd.DataFrame
        pandas dataframe of AAComp for protein sequence. Dataframe will
        be of the shape 20 x 1, where 20 is the number of features 
        calculated from the descriptor (for the 20 amino acids).
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

    amino_acid_composition = {}

    #iterate through each amino acid, calculating each composition
    for aa in amino_acids:
        amino_acid_composition[aa] = round(float(sequence.count(aa)) / len(sequence) * 100, 3)

    #transform values and columns to DataFrame
    amino_acid_composition_df = pd.DataFrame([list(amino_acid_composition.values())], 
        columns=list(amino_acid_composition.keys()))

    return amino_acid_composition_df

def dipeptide_composition(sequence):
    """
    Calculate Dipeptide Composition (DPComp) for protein sequence.
    Dipeptide composition is the fraction of each dipeptide type within a
    protein sequence. With dipeptides being of length 2 and there being 20
    canonical amino acids this creates 20^2 different combinations, thus a
    400-Dimensional vector will be produced such that:

    DPComp(s,t) = AA(s,t) / N -1

    where DPComp(s,t) is the dipeptide composition of the protein sequence
    for amino acid type s and t (where s and t = 1,2,..,20), AA(s,t) is the number
    of dipeptides represented by amino acid type s and t and N is the total number
    of dipeptides.

    Parameters
    ----------
    :sequence : str
        protein sequence.

    Returns
    -------
    :dipeptide_composition_df : pd.DataFrame
        pandas dataframe of dipeptide composition for protein sequence. Dataframe will
        be of the shape 400 x 1, where 400 is the number of features calculated 
        from the descriptor (20^2 for the 20 canonical amino acids).
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

    dipeptide_composition = {}

    #iterate through each amino acid, calculating each dipeptide combination's composition
    for i in amino_acids:
        for j in amino_acids:
            dipep = i + j
            dipeptide_composition[dipep] = round(
                float(sequence.count(dipep)) / (len(sequence)-1) * 100, 2)

    #transform values and columns to DataFrame
    dipeptide_composition_df = pd.DataFrame([list(dipeptide_composition.values())], 
        columns=list(dipeptide_composition.keys()))

    return dipeptide_composition_df

def tripeptide_composition(sequence):
    """
    Calculate Tripeptide Composition (TPComp) of protein sequence.
    Tripeptide composition is the fraction of each tripeptide type within a
    protein sequence. With tripeptides being of length 3 and there being 20
    canonical amino acids this creates 20^3 different combinations, thus a
    8000-Dimensional vector will be produced such that:

    TPComp(s,t,u) = AA(s,t,u) / N -1

    where TPComp(s,t,u) is the tripeptide composition of the protein sequence
    for amino acid type s, t and u (where s, t and u = 1,2,..,20), AA(s,t,u) is
    the number of tripeptides represented by amino acid type s and t, and N is
    the total number of tripeptides.

    Parameters
    ----------
    :sequence : str
        protein sequence in str form.

    Returns
    -------
    :tripeptide_composition_df : pd.DataFrame
        pandas DataFrame of tripeptide composition for protein sequence. Dataframe will
        be of the shape 8000 x 1, where 8000 is the number of features calculated 
        from the descriptor (20^3 for the 20 canonical amino acids).
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

    tripeptide_composition = {}
    tripeptides = []    

    #get list of tripeptides
    for i in amino_acids:
        for j in amino_acids:
            for k in amino_acids:
                tripeptides.append(i + j + k)

    #get frequency of each tripeptide in the sequence
    for i in tripeptides:
        tripeptide_composition[i] = len(re.findall(i, sequence))

    #transform values and columns to DataFrame
    tripeptide_composition_df = pd.DataFrame([list(tripeptide_composition.values())], 
        columns=list(tripeptide_composition.keys()))

    return tripeptide_composition_df

############################## Pseudo Amino Acid Composition ##############################

def pseudo_amino_acid_composition(sequence, lamda=30, weight=0.05, properties=[]):
    """
    Pseudo amino acid composition (PAAComp) combines the vanilla amino acid composition descriptor with 
    additional local features, such as correlation between residues of a certain distance, as
    amino acid composition doesn't take into accont sequence order info. The pseudo 
    components of the descriptor are a series rank-different correlation factors [5].
    The first 20 components are a weighted sum of the amino acid composition and 30 are 
    physiochemical square correlations as dictated by the lamda and properties parameters.
    This generates an output of [(20 + lamda), 1] = 50 x 1 when using the default lamda of 30. 
    By default, the physiochemical properties used are hydrophobicity and hydrophillicity, with 
    a lamda of 30 and weight of 0.05.

    Parameters
    ----------
    :lamda: int
        rank of correlation. Number of calculable descriptors depends on lamda.
    :weight: float
        weighting factor.
    :properties : list 
        list of dicts of physiochemical/structural property values for amino acids.

    Returns
    -------
    :pseudo_amino_acid_composition_df: pd.Dataframe
        output dataframe of calculated pseudo amino acid composition descriptors 
        for input sequence.
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

    #set lamda to its default value if <0, or > sequence len or not an int
    if ((lamda < 0) or (lamda > len(sequence)) or not isinstance(lamda, int)):
        lamda = 30  

    #validate weight, set default weight if invalid value input
    if ((weight<0) or not (isinstance(weight, float)) or not (isinstance(weight, int))):
        weight = 0.05

    #cast properties to list 
    if not (isinstance(properties, list)):
        properties = [properties]

    #by default use properties (hydrphocity, hydrophilicty, residue mass),
    #otherwise get the individual property values from aaindex1 using accession number
    if (properties == []):
        properties = [hydrophobicity_, hydrophilicity_, residue_mass_]
    else:
        for prop in properties:
            if (prop not in aaindex1.record_codes()):
                raise ValueError('Accession number not found in aaindex1: {}.'.format(prop))
            properties.append(aaindex1[prop].values)

    right_part = 0.0
    psuedo_aac = {}
    rightpart = []

    #calculate correlation factor for protein sequence using propeties
    for i in range(lamda):
        right_part = right_part + sequence_order_correlation_factor(sequence, i+1, properties)

    #calculate amino acid composition
    aa_comp = amino_acid_composition(sequence)
    
    #compute first 20 pseudo amino acid descriptor components based on properties,
    #append descriptor values to dict
    temp = 1 + weight + right_part
    for index, i in enumerate(amino_acids):
        psuedo_aac["PAAC_" + str(index + 1)] = round(aa_comp[i].values[0] / temp, 3)

    #calculate correlation factor for protein sequence using propeties
    for i in range(lamda):
        rightpart.append(sequence_order_correlation_factor(sequence, i + 1, properties))

    #compute last 20+lambda pseudo amino acid descriptor components based on properties,
    #append descriptor values to dict
    temp = 1 + weight * sum(rightpart)
    for index in range(20, 20 + lamda):
        psuedo_aac["PAAC_" + str(index + 1)] = round(
            weight * rightpart[index - 20] / temp * 100, 3)

    #transform values and columns to DataFrame
    psuedo_amino_acid_composition_df = pd.DataFrame([list(psuedo_aac.values())], columns=list(psuedo_aac.keys()))

    return psuedo_amino_acid_composition_df 

def sequence_order_correlation_factor(sequence, k=1, properties=[]):
    """
    Calculating sequence order correlation factor with gap equal to k based on
    the given input properities for a protein sequence.

    Parameters
    ----------
    :k : int
        gap between amino acids in the sequence.
    :properties : list
        list of dicts of physiochemical/structural property values for amino acids.

    Returns
    -------
    :result : float
        correlation factor value for the sequence using the correlation function
        for the particular sequence.
    """
    correlation_factor = []

    #iterate through sequence using correlation function to calculate sequence order
    #between its amino acids with gap=k 
    for i in range(len(sequence) - k):
        aa1 = sequence[i]
        aa2 = sequence[i + k]
        correlation_factor.append(correlation_function(aa1, aa2, properties))
    
    #round array of correlation factors 
    correlation_factor_ = round(sum(correlation_factor) / (len(sequence) - k), 3)

    return correlation_factor_    

def correlation_function(Ri="S", Rj="D", properties=[]):
    """
    Calculate the correlation between 2 input amino acids based on the 
    physiochemical/structural properties from the protein sequence.

    Parameters
    ----------
    :Ri : str
        1st amino acid.
    :Rj : str
        2nd amino acid.
    :properties : list 
        list of dicts of physiochemical/structural property values for amino acids.

    Returns
    -------
    :correlation : float
        correlation value for the two input amino acids based on input properties.
    """
    theta = 0.0

    #iterate over each property, calculating the correlation from the normalzied values
    for i in range(len(properties)):
        temp = normalize_property(properties[i])
        theta = theta + math.pow(temp[Ri] - temp[Rj], 2)    

    #calculate correlation for all properties
    correlation = round(theta / len(properties), 3)
    
    return correlation

def normalize_property(properties):
    """
    Normalize physiochemical/structural property values using
    their mean and standard deviation.
    
    Parameters
    ----------
    :properties : dict
        dictionary of amino acid values for a physiochemical property. 
        
    Returns
    -------
    :normalized_vals : dict
        dict of normalized property values.
    """
    normalized_vals = {}

    #iterate over all property values, calculate their mean and standard dev
    for i, j in properties.items():
        mean = sum(properties.values()) / len(properties.values())
        normalized_vals[i] = (j - mean) / _std(properties.values(), mean, ddof=0)

    return normalized_vals

def _std(prop, mean, ddof=1):
    """
    Calculate standard deviation from list of physiochemical values using
    the mean of their values.

    Parameters
    ----------
    :prop : dict
        dictionary of property values for each amino acid.
    :mean : float
        mean of physiochemical property values.
    :ddof : int (default=1)
        delta degrees of freedom
        
    Returns
    -------
    :std_ : float
        calculated standard deviation. 
    """
    temp = [math.pow(i - mean, 2) for i in prop]
    std_ = math.sqrt(sum(temp) / (len(prop) - ddof))
    return std_

######################## Amphiphilic Pseudo Amino Acid Composition ########################

def amphiphilic_pseudo_amino_acid_composition(sequence, lamda=30, weight=0.5, 
    properties=[hydrophobicity_, hydrophilicity_]):
    """
    Amphiphillic pseudo amino acid composition has the same form as the amino
    acid composition, but contains much more information that is related to the 
    sequence order of a protein and the distribution of the hydrophobic and 
    hydrophilic amino acids along its chain. The first 20 numbers in the 
    descriptor are the components of the conventional amino acid composition; 
    the next 2*λ numbers are a set of correlation factors that reflect different 
    hydrophobicity and hydrophilicity distribution patterns along a protein chain

    Parameters
    ----------
    :sequence : str
        protein sequence.
    :lamda: int
        rank of correlation. Number of calculable descriptors depends on lamda.
    :weight: float
        weighting factor.
    :properties : list (default=[hydrophobicity_, hydrophilicity_])
        list of dicts of physiochemical/structural property values for amino acids.

    Returns
    -------    
    :amp_pseudo_amino_acid_composition_df: pd.Dataframe
        output dataframe of calculated amphiphilic pseudo amino acid composition descriptors 
        for input sequence.
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

    #set lamda to its default value if <0, or > sequence len or not an int
    if ((lamda < 0) or (lamda > len(sequence)) or not isinstance(lamda, int)):
        lamda = 30  

    #validate weight, set default weight if invalid value input
    if ((weight<0) or not (isinstance(weight, float)) or not (isinstance(weight, int))):
        weight = 0.5

    #cast properties to list 
    if not (isinstance(properties, list)):
        properties = [properties]

    amp_pseudo_amino_acid_composition = {}

    #calculate correlation factor in protein sequence
    rightpart = 0.0
    for i in range(lamda):
        rightpart = rightpart + sum(
            amphiphilic_sequence_order_correllation_factor(sequence, k=i + 1))
    
    #calculate amino acid composition
    aa_composition = amino_acid_composition(sequence)

    #compute first 20 pseudo amino acid descriptor components based on properties
    temp = 1 + weight * rightpart
    for index, i in enumerate(amino_acids):
        amp_pseudo_amino_acid_composition["APAAC_" + str(index + 1)] = round(aa_composition[i].values[0] / temp, 3)


    #calculate correlation factor in protein sequence
    rightpart = []
    for i in range(lamda):
        temp = amphiphilic_sequence_order_correllation_factor(sequence, k=i + 1)
        rightpart.append(temp[0])
        rightpart.append(temp[1])

    #compute remaining (20 + 2 * lamda) amphiphillic composition values
    temp = 1 + weight * sum(rightpart)
    for index in range(20, 20 + 2 * lamda):
        amp_pseudo_amino_acid_composition["APAAC_" + str(index + 1)] = round(
            weight * rightpart[index - 20] / temp * 100, 3)

    #transform values and columns to DataFrame
    amp_pseudo_amino_acid_composition_df = pd.DataFrame([list(amp_pseudo_amino_acid_composition.values())], 
        columns=list(amp_pseudo_amino_acid_composition.keys()))


    return amp_pseudo_amino_acid_composition_df

def amphiphilic_sequence_order_correllation_factor(sequence, k=1):
    """ 
    Calculate Amphipilic sequence order correlation factor for sequence
    with gap=k.

    Parameters
    ----------
    :sequence : str
        protein sequence.
    :k : int (default=1)
        gap between amino acids in the sequence.

    Returns
    -------
    :correlation_factor : list
        list of correlation factors for both hydrophobicity and hydrophillicty.
    """
    correlation_factor_hydrophobicity = []
    correlation_factor_hydrophilicity = []

    #iterate over sequence, calculating correlation factors using both properties
    for i in range(len(sequence) - k):
        aa1 = sequence[i]
        aa2 = sequence[i + k]
        corelation_function = amphiphilic_correllation_function(aa1, aa2)
        correlation_factor_hydrophobicity.append(corelation_function[0])
        correlation_factor_hydrophilicity.append(corelation_function[1])

    #append hydrophobicity and hydrophilicity values to array
    correlation_factor = []
    correlation_factor.append(round(sum(correlation_factor_hydrophobicity) / (len(sequence) - k), 3))
    correlation_factor.append(round(sum(correlation_factor_hydrophilicity) / (len(sequence) - k), 3))

    return correlation_factor

def amphiphilic_correllation_function(Ri="S", Rj="D"):
    """
    Calculate correlation value based on hydrophobicity and hydrophillicity
    property values for input amino acids Ri and Rj in sequence.

    Parameters
    ----------
    :Ri : str
        1st amino acid.
    :Rj : str
        2nd amino acid.

    Returns
    -------
    :theta1, theta2 : float
        correlation values for input property values.
    """
    #normalize properties
    _hydrophobicity = normalize_property(hydrophobicity_)
    _hydrophilicity = normalize_property(hydrophilicity_)

    #calculate correlation using hydrophobicity & hydrophilicty properties
    theta1 = round(_hydrophobicity[Ri] * _hydrophobicity[Rj], 3)
    theta2 = round(_hydrophilicity[Ri] * _hydrophilicity[Rj], 3)

    return theta1, theta2