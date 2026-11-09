################################################################################
#############                Internal Helpers                     ##############
################################################################################

from __future__ import annotations

"""
Internal engines shared by the descriptor modules added in v1.5.0. Nothing in
this module is exported from the protpy namespace.
"""

import copy
import json
import os
import warnings
from collections import Counter
from functools import lru_cache

import numpy as np
import pandas as pd

from ._constants import amino_acids

# data directory for protpy JSON data files (aa_scales.json, aaindex.json)
_DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")


@lru_cache(maxsize=None)
def _read_json(filename: str) -> dict:
    """
    Read a JSON data file from the protpy/data directory. Files are cached
    after the first read.

    Parameters
    ==========
    :filename: str
        name of the JSON file in protpy/data.

    Returns
    =======
    :tables: dict
        contents of the JSON file.
    """
    with open(os.path.join(_DATA_DIR, filename)) as f:
        return json.load(f)


def _load_json_table(filename: str, key: str) -> dict:
    """
    Get a copy of a data table from a JSON file in the protpy/data directory.
    Each table has the keys name, values, reference, doi and source_note. A copy
    is returned so the cached table can't be modified.

    Parameters
    ==========
    :filename: str
        name of the JSON file in protpy/data.
    :key: str
        name of the table in the file.

    Returns
    =======
    :table: dict
        copy of the data table.

    Raises
    ======
    KeyError
        if the table is not in the file; the message lists the available tables.
    """
    tables = _read_json(filename)
    if key not in tables:
        raise KeyError(f"Table '{key}' not found in {filename}. Valid tables are: {sorted(tables)}.")
    #copy so callers cannot mutate the cached table
    return copy.deepcopy(tables[key])


def _scale_values(scale: str | dict) -> dict:
    """
    Get the values of an amino acid scale for each amino acid. The scale can be a
    table name in the aa_scales.json data file, an AAindex accession number or a
    dict of amino acid values.

    Parameters
    ==========
    :scale: str | dict
        scale name, AAindex accession number or dict of amino acid values.

    Returns
    =======
    :values: dict
        dict of scale values (a float or a list of floats) for each amino acid;
        a copy for protpy tables, so callers cannot change the cached data.

    Raises
    ======
    ValueError
        if the scale is not a protpy scale or an AAindex accession number.
    """
    if isinstance(scale, dict):
        return scale
    tables = _read_json("aa_scales.json")
    if scale == "eisenberg":
        _warn_deprecated("scale='eisenberg'", "scale='eisenberg_arb1984' (the same table, AAindex EISD840101)")
        scale = "eisenberg_arb1984"
    if scale in tables:
        if tables[scale].get("experimental"):
            warnings.warn(f"Scale '{scale}' is experimental: some of its values could not be checked against "
                          "the primary paper (see 'source_note' in protpy/data/aa_scales.json).",
                          UserWarning, stacklevel=3)
        return copy.deepcopy(tables[scale]["values"])
    from .autocorrelation import _get_properties
    try:
        return _get_properties(scale)[scale]
    except ValueError:
        raise ValueError(f"Unknown scale '{scale}'. Use a protpy scale name {sorted(tables)}, "
                         "an AAindex1 accession code, or a {amino_acid: value} dict.") from None


def _property_matrix(sequence: str, scales: str | list[str] | dict, standardize: bool = False) -> np.ndarray:
    """
    Convert a protein sequence into a matrix of amino acid property values, with
    one row per residue and one column per property dimension. The scales can be a
    factor-scale family or single-value scale from the aa_scales.json data file,
    an AAindex accession number, a list of AAindex accession numbers or a dict of
    amino acid values (a float or a list of floats per amino acid).

    Parameters
    ==========
    :sequence: str
        protein sequence.
    :scales: str | list[str] | dict
        scale name, AAindex accession number, list of accession numbers or dict of
        amino acid values.
    :standardize: bool (default=False)
        standardise each dimension to zero mean and unit standard deviation over
        the 20 amino acids.

    Returns
    =======
    :matrix: np.ndarray
        property values of the sequence, of the shape L x d, where L is the length
        of the sequence and d is the number of property dimensions.

    Raises
    ======
    ValueError
        if the scale is missing values for any of the 20 amino acids.
    """
    if isinstance(scales, (list, tuple)):
        tables = [_scale_values(s) for s in scales]
        table = {aa: [t[aa] for t in tables] for aa in amino_acids}
    else:
        table = _scale_values(scales)
    missing = [aa for aa in amino_acids if aa not in table]
    if missing:
        raise ValueError(f"Scale is missing values for amino acids: {missing}.")
    #rows ordered by amino_acids so standardisation is over exactly the 20 residues
    aa_matrix = np.array([np.atleast_1d(np.asarray(table[aa], dtype=float)) for aa in amino_acids])
    if standardize:
        std = aa_matrix.std(axis=0)
        aa_matrix = np.divide(aa_matrix - aa_matrix.mean(axis=0), std,
                              out=np.zeros_like(aa_matrix), where=std != 0)
    index = {aa: i for i, aa in enumerate(amino_acids)}
    return aa_matrix[[index[aa] for aa in sequence]]


def _sliding_windows(values: np.ndarray, window: int, edge: float = 1.0) -> np.ndarray:
    """
    Calculate the weighted mean of a profile over a sliding window. With edge=1.0
    the weights are all 1 (a plain mean); with edge < 1 the weights fall linearly
    from 1 at the centre of the window to edge at its ends, as in ExPASy ProtScale
    and Biopython's protein_scale. A window longer than the profile uses the whole
    profile.

    Parameters
    ==========
    :values: np.ndarray
        profile values along the sequence.
    :window: int
        length of the sliding window.
    :edge: float (default=1.0)
        relative weight of the ends of the window, between 0 and 1.

    Returns
    =======
    :windows: np.ndarray
        weighted mean of each window, of length L - window + 1.
    """
    values = np.asarray(values, dtype=float)
    window = min(window, len(values))
    half = (window - 1) / 2
    dist = np.minimum(np.arange(window), np.arange(window)[::-1])
    weights = edge + (1 - edge) * dist / half if half else np.ones(window)
    weights = np.minimum(weights, 1.0)
    return np.convolve(values, weights[::-1], mode="valid") / weights.sum()


def _profile_stats(profile: np.ndarray, prefix: str, stats=("mean", "max", "min", "std"),
                   threshold: float | None = None) -> dict:
    """
    Summarise a profile by its mean, maximum, minimum and standard deviation, and
    optionally the fraction of values above a threshold.

    Parameters
    ==========
    :profile: np.ndarray
        profile values.
    :prefix: str
        prefix for the names of the statistics.
    :stats: tuple (default=("mean", "max", "min", "std"))
        statistics to calculate.
    :threshold: float | None (default=None)
        if set, also calculate the fraction of values above it (frac_above).

    Returns
    =======
    :stats: dict
        dict of statistic names (prefix_stat) and values.
    """
    profile = np.asarray(profile, dtype=float)
    funcs = {"mean": np.mean, "max": np.max, "min": np.min, "std": np.std}
    out = {f"{prefix}_{s}": float(funcs[s](profile)) for s in stats}
    if threshold is not None:
        out[f"{prefix}_frac_above"] = float(np.mean(profile > threshold))
    return out


def _gapped_pair_counts(sequence: str, gap: int, alphabet: dict[str, str] | None = None) -> Counter:
    """
    Count the pairs of residues (s_i, s_i+gap+1) in a protein sequence, with gap
    residues between them. The residues can first be mapped to groups using a
    reduced alphabet.

    Parameters
    ==========
    :sequence: str
        protein sequence.
    :gap: int
        number of residues between the two residues of a pair.
    :alphabet: dict | None (default=None)
        dict mapping each amino acid to its group; None counts the amino acids.

    Returns
    =======
    :counts: Counter
        count of each pair, keyed by a tuple of the two residues or groups.
    """
    seq = [alphabet[aa] for aa in sequence] if alphabet else list(sequence)
    step = gap + 1
    return Counter((seq[i], seq[i + step]) for i in range(len(seq) - step))


def _kmer_counts(sequence: str, k: int, alphabet: dict[str, str] | None = None) -> Counter:
    """
    Count the overlapping k-mers of a protein sequence. The residues can first be
    mapped to groups using a reduced alphabet.

    Parameters
    ==========
    :sequence: str
        protein sequence.
    :k: int
        length of the k-mers.
    :alphabet: dict | None (default=None)
        dict mapping each amino acid to its group; None counts the amino acids.

    Returns
    =======
    :counts: Counter
        count of each k-mer, keyed by a tuple of its residues or groups.
    """
    seq = [alphabet[aa] for aa in sequence] if alphabet else list(sequence)
    return Counter(tuple(seq[i:i + k]) for i in range(len(seq) - k + 1))


def _charge_vector(sequence: str, his: bool = False) -> np.ndarray:
    """
    Get the charge of each residue of a protein sequence: +1 for K and R, -1 for
    D and E and 0 otherwise. His is neutral unless his=True.

    Parameters
    ==========
    :sequence: str
        protein sequence.
    :his: bool (default=False)
        treat His as positive.

    Returns
    =======
    :charges: np.ndarray
        charge of each residue.
    """
    positive = "KRH" if his else "KR"
    return np.array([1 if aa in positive else -1 if aa in "DE" else 0 for aa in sequence], dtype=float)


def _runs(mask: np.ndarray, min_len: int = 1) -> list[tuple[int, int]]:
    """
    Find the runs of consecutive True values in a boolean array, e.g. the
    segments of a profile above a threshold.

    Parameters
    ==========
    :mask: np.ndarray
        boolean values along the sequence.
    :min_len: int (default=1)
        minimum length of a run.

    Returns
    =======
    :runs: list
        list of (start, end) indices of each run, where end is exclusive.
    """
    mask = np.concatenate(([False], np.asarray(mask, dtype=bool), [False]))
    edges = np.flatnonzero(np.diff(mask.astype(int)))
    return [(int(s), int(e)) for s, e in zip(edges[::2], edges[1::2]) if e - s >= min_len]


def _warn_deprecated(old: str, new: str | None = None) -> None:
    """
    Emit a DeprecationWarning for a descriptor that will be removed in protpy 2.0.

    Parameters
    ==========
    :old: str
        name of the deprecated function or option.
    :new: str | None (default=None)
        what to use instead, if anything.
    """
    msg = f"{old} is deprecated and will be removed in protpy 2.0"
    warnings.warn(msg + (f"; use {new} instead." if new else "."), DeprecationWarning, stacklevel=3)


def _check_lag(name: str, lag, length: int, minimum: int = 1) -> int:
    """
    Validate a lag (or lamda/gap) against the sequence length. A lag needs at
    least one residue pair, so it must be less than the length of the sequence.
    This is the single lag policy of protpy: an invalid lag raises a ValueError
    rather than being replaced by a default.

    Parameters
    ==========
    :name: str
        name of the parameter, used in the error message.
    :lag: int
        lag value.
    :length: int
        length of the sequence.
    :minimum: int (default=1)
        smallest allowed lag.

    Returns
    =======
    :lag: int
        validated lag.
    """
    from ._constants import _validate_int
    lag = _validate_int(name, lag, minimum)
    if lag >= length:
        raise ValueError(f"{name} ({lag}) must be less than the sequence length ({length}); "
                         f"pass a smaller {name} for short sequences.")
    return lag


def _to_df(values: list, columns: list[str]) -> pd.DataFrame:
    """
    Convert a list of descriptor values into a single-row dataframe. Values are
    not rounded (rounding is opt-in through the ``decimals`` argument of
    all_descriptors and encode_dataset); booleans and numpy integers are converted
    to ints and numpy floats to Python floats.

    Parameters
    ==========
    :values: list
        descriptor values.
    :columns: list[str]
        column names.

    Returns
    =======
    :descriptor_df: pd.DataFrame
        pandas dataframe of the descriptor values, of the shape 1 x N, where N is
        the number of values.
    """
    values = [float(v) if isinstance(v, (float, np.floating)) else
              int(v) if isinstance(v, (bool, np.bool_, np.integer)) else v for v in values]
    return pd.DataFrame([values], columns=columns)
