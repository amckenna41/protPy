################################################################################
#############                   Autocorrelation                   ##############
################################################################################

from __future__ import annotations

import pandas as pd
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

import numpy as np

#list of amino acids
from ._constants import amino_acids, _validate_sequence
from ._helpers import _check_lag, _property_matrix, _to_df

__all__ = ["autocorrelation", "moreaubroto_autocorrelation", "moran_autocorrelation", "geary_autocorrelation"]

#default AAindex1 properties of the autocorrelation descriptors:
#CIDH920105 normalized average hydrophobicity, BHAR880101 average flexibility,
#CHAM820101 polarizability, CHAM820102 free energy of solution in water,
#CHOC760101 residue accessible surface area in tripeptide, BIGC670101 residue volume,
#CHAM810101 steric parameter, DAYM780201 relative mutability
default_properties_ = ("CIDH920105", "BHAR880101", "CHAM820101", "CHAM820102",
    "CHOC760101", "BIGC670101", "CHAM810101", "DAYM780201")

#column prefix of each method
_methods = {"moreaubroto": "MBAuto", "moran": "MAuto", "geary": "GAuto"}

################################### Autocorrelation #######################################

def autocorrelation(sequence: str, method: str = "geary", lag: int = 30,
    properties: list[str] | tuple[str, ...] = default_properties_, normalize: bool = True) -> pd.DataFrame:
    """
    Calculate an Autocorrelation descriptor of protein sequence. Autocorrelation
    descriptors are topological descriptors that describe the level of correlation
    between residues separated by d positions in terms of a physicochemical
    property [1, 2]. For a property P (by default standardised over the 20 amino
    acids) and each lag d = 1..lag they are calculated as:

    Moreau-Broto: MB(d) = sum_{i=1}^{L-d} P_i * P_{i+d} / (L - d)

    Moran:        I(d)  = [sum_{i=1}^{L-d} (P_i - P') * (P_{i+d} - P') / (L - d)]
                          / [sum_{i=1}^{L} (P_i - P')^2 / L]

    Geary:        C(d)  = [sum_{i=1}^{L-d} (P_i - P_{i+d})^2 / (2 (L - d))]
                          / [sum_{i=1}^{L} (P_i - P')^2 / (L - 1)]

    where P_i is the property value of the residue at position i, P' is the mean
    of P over the sequence and L is the length of the sequence. Moran and Geary
    are 0 for a sequence whose property values do not vary. With the default 8
    AAindex properties and lag=30, 240 features are generated. Ong et al. found no
    clearly preferred method among the three [2]; Moran and Geary are strongly
    (negatively) correlated at the same property and lag, and Moreau-Broto equals
    ``auto_cross_covariance(mode="ac")`` on the same properties. The method names
    are also available as ``moreaubroto_autocorrelation``,
    ``moran_autocorrelation`` and ``geary_autocorrelation``.

    Parameters
    ==========
    :sequence: str
        protein sequence.
    :method: str (default="geary")
        "moreaubroto", "moran" or "geary".
    :lag: int (default=30)
        maximum lag d, less than the length of the sequence.
    :properties: list[str] (default=the 8 properties above)
        AAindex1 accession numbers (or protpy scale names) of the properties.
    :normalize: bool (default=True)
        standardise each property (zero mean, unit standard deviation over the 20
        amino acids) before the calculation. This does not bound the output to [0, 1].

    Returns
    =======
    :autocorrelation_df: pd.DataFrame
        pandas Dataframe of autocorrelation values for protein sequence, of the
        shape 1 x (number of properties * lag); by default 1 x 240. Column names
        are of the form <prefix>_<property>_<d>, with the prefix MBAuto, MAuto or
        GAuto.
    """
    #validate input protein sequence
    sequence = _validate_sequence(sequence)
    if method not in _methods:
        raise ValueError(f"Invalid method '{method}'. Valid methods are: {list(_methods)}.")
    lag = _check_lag("lag", lag, len(sequence))
    properties = [properties] if isinstance(properties, str) else list(properties or [])
    if not all(properties) or not properties:
        raise ValueError("At least one property must be passed into function.")

    p = _property_matrix(sequence, properties, standardize=normalize)
    L = len(sequence)
    centred = p - p.mean(axis=0)
    var = (centred ** 2).mean(axis=0)
    values = np.zeros((len(properties), lag))
    for d in range(1, lag + 1):
        if method == "moreaubroto":
            values[:, d - 1] = np.mean(p[:-d] * p[d:], axis=0)
        elif method == "moran":
            values[:, d - 1] = np.divide(np.mean(centred[:-d] * centred[d:], axis=0), var,
                                         out=np.zeros(len(properties)), where=var > 0)
        else:
            values[:, d - 1] = np.divide(np.mean((p[:-d] - p[d:]) ** 2, axis=0) / 2, var * L / (L - 1),
                                         out=np.zeros(len(properties)), where=var > 0)
    prefix = _methods[method]
    columns = [f"{prefix}_{prop}_{d}" for prop in properties for d in range(1, lag + 1)]
    return _to_df(list(values.ravel()), columns)

def moreaubroto_autocorrelation(sequence: str, lag: int = 30,
    properties: list[str] | tuple[str, ...] = default_properties_, normalize: bool = True) -> pd.DataFrame:
    """
    Calculate the Moreau-Broto Autocorrelation (MBAuto) of protein sequence, which
    uses the products of the property values [2]. Permanent alias of
    ``autocorrelation(sequence, method="moreaubroto", ...)``; see ``autocorrelation``
    for the formula and parameters.

    Returns
    =======
    :moreaubroto_autocorrelation_df: pd.Dataframe
        pandas Dataframe of shape 1 x (number of properties * lag), by default
        1 x 240, with columns MBAuto_<property>_<d>.
    """
    return autocorrelation(sequence, "moreaubroto", lag, properties, normalize)

def moran_autocorrelation(sequence: str, lag: int = 30,
    properties: list[str] | tuple[str, ...] = default_properties_, normalize: bool = True) -> pd.DataFrame:
    """
    Calculate the Moran Autocorrelation (MAuto) of protein sequence, which uses the
    deviations of the property values from their mean [2]. Permanent alias of
    ``autocorrelation(sequence, method="moran", ...)``; see ``autocorrelation`` for
    the formula and parameters.

    Returns
    =======
    :moran_autocorrelation_df: pd.DataFrame
        pandas Dataframe of shape 1 x (number of properties * lag), by default
        1 x 240, with columns MAuto_<property>_<d>.
    """
    return autocorrelation(sequence, "moran", lag, properties, normalize)

def geary_autocorrelation(sequence: str, lag: int = 30,
    properties: list[str] | tuple[str, ...] = default_properties_, normalize: bool = True) -> pd.DataFrame:
    """
    Calculate the Geary Autocorrelation (GAuto) of protein sequence, which uses the
    squared differences of the property values [2]. Permanent alias of
    ``autocorrelation(sequence, method="geary", ...)``; see ``autocorrelation`` for
    the formula and parameters. Geary is the autocorrelation in the classic
    ``all_descriptors`` set.

    Returns
    =======
    :geary_autocorrelation_df: pd.DataFrame
        pandas Dataframe of shape 1 x (number of properties * lag), by default
        1 x 240, with columns GAuto_<property>_<d>.
    """
    return autocorrelation(sequence, "geary", lag, properties, normalize)

def _get_properties(properties: list[str] | str) -> dict[str, dict[str, float]]:
    """
    Get the values of the 20 canonical amino acids for each AAIndex1 accession
    number. Values are copied so the aaindex records themselves are never
    modified. Used by protpy's scale lookup (_helpers._scale_values); the
    autocorrelation descriptors standardise the values there.

    Parameters
    ==========
    :properties: list[str] | str
        AAIndex1 accession number(s).

    Returns
    =======
    :aai_properties: dict
        {accession_number: {amino_acid: value}}.
    """
    if isinstance(properties, str):
        properties = [properties]
    record_codes = aaindex1.record_codes()
    aai_properties = {}
    for prop in properties:
        if prop not in record_codes:
            raise ValueError(f"Property {prop} not found in list of available properties in the aaindex1:\n{properties}.")
        #keep only the 20 canonical amino acids (drops the '-' gap entry)
        aai_properties[prop] = {aa: v for aa, v in aaindex1[prop].values.items() if aa in amino_acids}
    return aai_properties
