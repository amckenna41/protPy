################################################################################
#############                   Autocorrelation                   ##############
################################################################################

from __future__ import annotations

import pandas as pd
import math
from aaindex import aaindex1

"""
References
==========
[1] B. Hollas, “An analysis of the autocorrelation descriptor for molecules,”
    J. Math. Chem., vol. 33, no. 2, pp. 91–101, 2003.

[2] S. A. K. Ong, H. H. Lin, Y. Z. Chen, Z. R. Li, and Z. Cao, “Efficacy
    of different protein descriptors in predicting protein functional families,”
    BMC Bioinformatics, vol. 8, p. 300, 2007.

[3] D. Raimondi, G. Orlando, W. F. Vranken, and Y. Moreau,
    “Exploring the limitations of biophysical propensity scales coupled with
    machine learning for protein sequence analysis,” Sci. Rep., vol. 9, no. 1, p. 16932, 2019.  
"""

#list of amino acids
from ._constants import amino_acids, _validate_sequence, _validate_int
from ._helpers import _minmax
    
############################### MoreauBroto Autocorrelation ###############################

def moreaubroto_autocorrelation(sequence: str, lag: int = 30,
    properties: list[str] = ["CIDH920105", "BHAR880101", "CHAM820101", "CHAM820102",
    "CHOC760101", "BIGC670101", "CHAM810101", "DAYM780201"], normalize: bool = True,
    normalized: bool = False) -> pd.DataFrame:
    """
    Calculate MoreauBrotoAuto Autocorrelation (MBAuto) descriptor for sequence.
    Autocorrelation descriptors are a class of topological descriptors,
    also known as molecular connectivity indices, that describe the level of
    correlation between two objects (protein or peptide sequences) in terms of
    their specific structural or physicochemical properties, which are
    defined based on the distribution of amino acid properties along the sequence.
    By default, 8 amino acid properties are used for deriving the descriptors. 
    The derivations and detailed explanations of this type of descriptor is 
    outlined in [1]. 
    
    The MBAuto descriptor is a type of Autocorrelation descriptor that uses the
    property values as the basis for measurement [2]. Each autocorrelation will 
    generate the number of features depending on the lag value and number of 
    properties input with total features = lag * number of properties. The output 
    autocorrelation can also be normalized by setting the normalize parameter to 
    true, this occurs by default. Using the default 8 properties with default lag 
    value of 30, 240 features are generated, the default 8 properties are:

    AccNo. CIDH920105 - Normalized Average Hydrophobicity Scales
    AccNo. BHAR880101 - Average Flexibility Indices
    AccNo. CHAM820101 - Polarizability Parameter
    AccNo. CHAM820102 - Free Energy of Solution in Water, kcal/mole
    AccNo. CHOC760101 - Residue Accessible Surface Area in Tripeptide
    AccNo. BIGC670101 - Residue Volume
    AccNo. CHAM810101 - Steric Parameter
    AccNo. DAYM780201 - Relative Mutability

    Parameters
    ==========
    :sequence: str
        protein sequence.
    :lag: int (default=30)
        A value for a lag, the max value is equal to the length of the shortest peptide minus one.
    :properties: list
        list of AAI index record codes/accession numbers for the physicochemical properties to use in the calculation. Defaults to the 8 standard AAIndex properties listed above.
    :normalize: bool (default=True)
        standardise each property (zero mean, unit standard deviation over the 20
        amino acids) before the calculation. This does not bound the output; use
        ``normalized`` for [0, 1] output.
    :normalized: bool (default=False)
        rescale every value to [0, 1] with its theoretical bounds (see below);
        columns get a ``_norm`` suffix. Independent of ``normalize``, which
        standardises the property values before the calculation.
        Bounds: the minimum and maximum product P(a) x P(b) over all amino-acid
        pairs (extended to include 0, the value of lags with no residue pairs),
        since each value is a mean of such products.

    Returns
    =======
    :moreaubroto_autocorrelation_df: pd.Dataframe
        pandas Dataframe of MBAuto values for protein sequence. Output will
        be of the shape 1 x N, where N is the number of features calculated from
        the descriptor and 1 is the input sequence. By default, the shape will
        be 1 x 240 (30 features per property - using 8 properties and lag=30).
    """
    #validate input protein sequence
    sequence = _validate_sequence(sequence)

    #validate lag; fall back to the default of 30 when lag exceeds the sequence length
    lag = _validate_int("lag", lag)
    if lag >= len(sequence):
        lag = 30

    #get (optionally normalized) property values for each AAIndex accession number
    aai_properties = _get_properties(properties, normalize)

    #store resultant autocorrelation values
    moreaubroto_autocorrelation = {}

    #iterate through list of properties, calculating autocorrelation values for the sequence, append to results dict
    for key in aai_properties:
        for i in range(1, lag+1):
            temp = 0  # reset per lag
            for j in range(len(sequence)-i):
                temp = temp + aai_properties[key][sequence[j]] * aai_properties[key][sequence[j+i]]
            #lags with no residue pairs (lag >= sequence length) are 0
            moreaubroto_autocorrelation["MBAuto_" + key + "_" + str(i)] = (
                round(temp / (len(sequence) - i), 3) if len(sequence) > i else 0.0)

    #transform values and columns to DataFrame
    moreaubroto_autocorrelation_df = pd.DataFrame([list(moreaubroto_autocorrelation.values())], 
        columns=list(moreaubroto_autocorrelation.keys()))

    if normalized:
        lo, hi = [], []
        for key in aai_properties:
            products = [aai_properties[key][a] * aai_properties[key][b] for a in amino_acids for b in amino_acids]
            #include 0: lags with no residue pairs are reported as 0
            lo += [min(min(products), 0.0)] * lag
            hi += [max(max(products), 0.0)] * lag
        return _minmax(moreaubroto_autocorrelation_df, lo, hi)
    return moreaubroto_autocorrelation_df

################################## Moran Autocorrelation ##################################

def moran_autocorrelation(sequence: str, lag: int = 30,
    properties: list[str] = ["CIDH920105", "BHAR880101", "CHAM820101", "CHAM820102",
    "CHOC760101", "BIGC670101", "CHAM810101", "DAYM780201"], normalize: bool = True,
    normalized: bool = False) -> pd.DataFrame:
    """
    Refer to MBAuto docstring for autocorrelation explanation.
    Moran autocorrelation (MAuto) utilizes property deviations from the
    average values [2].

    Parameters
    ==========
    :sequence: str
        protein sequence.
    :lag: int (default=30)
        A value for a lag, the max value is equal to the length of the shortest peptide minus one.
    :properties: list
        list of AAI index record codes/accession numbers for the physicochemical properties to use in the calculation. Defaults to the same 8 standard AAIndex properties as MBAuto.
    :normalize: bool (default=True)
        standardise each property (zero mean, unit standard deviation over the 20
        amino acids) before the calculation. This does not bound the output; use
        ``normalized`` for [0, 1] output.
    :normalized: bool (default=False)
        rescale every value to [0, 1] with its theoretical bounds (see below);
        columns get a ``_norm`` suffix. Independent of ``normalize``, which
        standardises the property values before the calculation.
        Bounds: abs(I(d)) <= N / (N - d) by the Cauchy-Schwarz inequality, so I is
        mapped from [-N/(N-d), N/(N-d)]; 0.5 means no autocorrelation.

    Returns
    =======
    :moran_autocorrelation_df: pd.DataFrame
        pandas Dataframe of MAuto values for protein sequence. Output will
        be of the shape 1 x N, where N is the number of features calculated from
        the descriptor and 1 is the input sequence. By default, the shape will 
        be 1 x 240 (30 features per property - using 8 properties and lag=30).  
    """
    #validate input protein sequence
    sequence = _validate_sequence(sequence)

    #validate lag; fall back to the default of 30 when lag exceeds the sequence length
    lag = _validate_int("lag", lag)
    if lag >= len(sequence):
        lag = 30

    #get (optionally normalized) property values for each AAIndex accession number
    aai_properties = _get_properties(properties, normalize)

    #store resultant autocorrelation values
    moran_autocorrelation = {}

    #iterate through list of properties, calculating autocorrelation values for the sequence, append to results dict
    for key in aai_properties:
        cds = 0
        for aa in amino_acids:
            cds = cds + sequence.count(aa) * aai_properties[key][aa]
        prop_mean = cds / len(sequence)
        cc = []
        for aa in sequence:
            cc.append(aai_properties[key][aa])
        k = (_std(cc, ddof=0)) ** 2
        for i in range(1, lag+1):
            temp = 0
            for j in range(len(sequence) - i):
                temp = temp + (aai_properties[key][sequence[j]] - prop_mean) * (
                    aai_properties[key][sequence[j+i]] - prop_mean)
            #Moran: (1/(N-d) * sum) / variance; 0 when there are no pairs or no variance
            moran_autocorrelation["MAuto_" + key + "_" + str(i)] = (
                round(temp / (len(sequence) - i) / k, 5) if (len(sequence) > i and k > 0) else 0.0)

    #transform values and columns to DataFrame
    moran_autocorrelation_df = pd.DataFrame([list(moran_autocorrelation.values())], 
        columns=list(moran_autocorrelation.keys()))

    if normalized:
        N = len(sequence)
        bound = [N / (N - i) if N > i else 1.0 for i in range(1, lag + 1)] * len(aai_properties)
        return _minmax(moran_autocorrelation_df, [-b for b in bound], bound)
    return moran_autocorrelation_df

################################## Geary Autocorrelation ##################################

def geary_autocorrelation(sequence: str, lag: int = 30,
    properties: list[str] = ["CIDH920105", "BHAR880101", "CHAM820101", "CHAM820102",
    "CHOC760101", "BIGC670101", "CHAM810101", "DAYM780201"], normalize: bool = True,
    normalized: bool = False) -> pd.DataFrame:
    """
    Refer to MBAuto docstring for autocorrelation description.
    Geary Autocorrelation (GAuto) utilizes the square-difference of property
    values instead of vector-products (of property values or deviations) [2].
    
    Parameters
    ==========
    :sequence: str
        protein sequence.
    :lag: int (default=30)
        A value for a lag, the max value is equal to the length of the shortest peptide minus one.
    :properties: list
        list of AAI index record codes/accession numbers for the physicochemical properties to use in the calculation. Defaults to the same 8 standard AAIndex properties as MBAuto.
    :normalize: bool (default=True)
        standardise each property (zero mean, unit standard deviation over the 20
        amino acids) before the calculation. This does not bound the output; use
        ``normalized`` for [0, 1] output.
    :normalized: bool (default=False)
        rescale every value to [0, 1] with its theoretical bounds (see below);
        columns get a ``_norm`` suffix. Independent of ``normalize``, which
        standardises the property values before the calculation.
        Bounds: 0 <= C(d) <= 2(N - 1) / (N - d), since each squared difference is at
        most twice the sum of the two squared deviations.

    Returns
    =======
    :geary_autocorrelation_df: pd.DataFrame
        pandas Dataframe of MAuto values for protein sequence. Output will
        be of the shape 1 x N, where N is the number of features calculated from
        the descriptor and 1 is the input sequence. By default, the shape will 
        be 1 x 240 (30 features per property - using 8 properties and lag=30).  
    """
    #validate input protein sequence
    sequence = _validate_sequence(sequence)
            
    #validate lag; fall back to the default of 30 when lag exceeds the sequence length
    lag = _validate_int("lag", lag)
    if lag >= len(sequence):
        lag = 30

    #get (optionally normalized) property values for each AAIndex accession number
    aai_properties = _get_properties(properties, normalize)

    #store resultant autocorrelation values
    geary_autocorrelation = {}

    #iterate through list of properties, calculating autocorrelation values for the sequence, append to results dict
    for key in aai_properties:
        cc = []
        for aa in sequence:
            cc.append(aai_properties[key][aa])
        k = ((_std(cc, ddof=0)) ** 2) * len(sequence) / (len(sequence) - 1) if len(sequence) > 1 else 0.0
        for i in range(1, lag+1):
            temp = 0
            for j in range(len(sequence) - i):
                temp = (temp + (
                    aai_properties[key][sequence[j]] - aai_properties[key][sequence[j+i]]) ** 2)
            #0 when there are no pairs or no variance
            geary_autocorrelation["GAuto_" + key + "_" + str(i)] = (
                round(temp / (2* (len(sequence) -i)) / k, 3) if (len(sequence) > i and k > 0) else 0.0)

    #transform values and columns to DataFrame
    geary_autocorrelation_df = pd.DataFrame([list(geary_autocorrelation.values())], 
        columns=list(geary_autocorrelation.keys()))

    if normalized:
        N = len(sequence)
        bound = [2 * (N - 1) / (N - i) if N > i else 1.0 for i in range(1, lag + 1)] * len(aai_properties)
        return _minmax(geary_autocorrelation_df, 0.0, bound)
    return geary_autocorrelation_df

def _get_properties(properties: list[str] | str, normalize: bool) -> dict[str, dict[str, float]]:
    """
    Get the values of the 20 canonical amino acids for each AAIndex1 accession number,
    optionally normalized to zero mean and unit standard deviation. Values are copied
    so the aaindex records themselves are never modified.

    Parameters
    ==========
    :properties: list[str] | str
        AAIndex1 accession number(s).
    :normalize: bool
        normalize the property values.

    Returns
    =======
    :aai_properties: dict
        {accession_number: {amino_acid: value}}.
    """
    #validate that at least 1 property input to function
    if (properties == "" or properties == [] or properties is None):
        raise ValueError('At least one property must be passed into function.')

    #if 1 property passed into function, wrap it to an iterable list
    if (isinstance(properties, str)):
        properties = [properties]

    record_codes = aaindex1.record_codes()
    aai_properties = {}
    for prop in properties:
        #raise value error if accession number not found in aaindex1
        if prop not in record_codes:
            raise ValueError(f"Property {prop} not found in list of available properties in the aaindex1:\n{properties}.")
        #keep only the 20 canonical amino acids (drops the '-' gap entry)
        vals = {aa: v for aa, v in aaindex1[prop].values.items() if aa in amino_acids}
        if normalize:
            mean = sum(vals.values()) / len(vals)
            std = _std(list(vals.values()), ddof=0)
            vals = {aa: (v - mean) / std if std else 0.0 for aa, v in vals.items()}
        aai_properties[prop] = vals

    return aai_properties

def _std(array: list[float], ddof: int = 1) -> float:
    """
    Calculate the standard deviation of the array data, with associated means Delta Degrees
    of Freedom (ddof).

    Parameters
    ==========
    :array: np.array
        numpy array of float values.
    :ddof: int (default=1)
        means delta degrees of freedom.

    Returns
    =======
    :result: np.array
        input array after standard deviation transformation.
    """
    mean = sum(array) / len(array)
    return math.sqrt(sum(math.pow(v - mean, 2) for v in array) / (len(array) - ddof))