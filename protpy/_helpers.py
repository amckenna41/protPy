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
from collections import Counter
from functools import lru_cache

import numpy as np
import pandas as pd

from ._constants import amino_acids

_DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")


@lru_cache(maxsize=None)
def _read_json(filename: str) -> dict:
    with open(os.path.join(_DATA_DIR, filename)) as f:
        return json.load(f)


def _load_json_table(filename: str, key: str) -> dict:
    """
    Return a copy of table ``key`` from ``protpy/data/<filename>``. Each table
    has the keys name, values, reference, doi and source_note.

    Raises
    ======
    KeyError
        if ``key`` is not a table in the file; the message lists the valid keys.
    """
    tables = _read_json(filename)
    if key not in tables:
        raise KeyError(f"Table '{key}' not found in {filename}. Valid tables are: {sorted(tables)}.")
    #copy so callers cannot mutate the cached table
    return copy.deepcopy(tables[key])


def _scale_values(scale: str | dict) -> dict:
    """
    Resolve a scale to {aa: float | list[float]}. ``scale`` may be a table name
    in aa_scales.json, an AAindex1 accession code, or a user dict.
    """
    if isinstance(scale, dict):
        return scale
    tables = _read_json("aa_scales.json")
    if scale in tables:
        return tables[scale]["values"]
    from .autocorrelation import _get_properties
    try:
        return _get_properties(scale, normalize=False)[scale]
    except ValueError:
        raise ValueError(f"Unknown scale '{scale}'. Use a protpy scale name {sorted(tables)}, "
                         "an AAindex1 accession code, or a {amino_acid: value} dict.") from None


def _property_matrix(sequence: str, scales: str | list[str] | dict, standardize: bool = False) -> np.ndarray:
    """
    Return the L x d matrix of property values for ``sequence``.

    ``scales`` may be a factor-scale family / single-value scale name from
    aa_scales.json, an AAindex1 accession code, a list of AAindex1 codes, or a
    user dict {aa: float} / {aa: list[float]}. With ``standardize=True`` each
    dimension is z-scored over the 20 canonical amino acids (population std).
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
    Weighted sliding-window means of a 1-D profile.

    Weights fall linearly from ``edge`` at the window ends to 1 at the centre
    (ExPASy ProtScale / Biopython ``protein_scale`` semantics); ``edge=1`` is a
    plain mean. If ``window`` exceeds the profile length a single window of the
    full length is used.
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
    """Reduce a profile to {prefix_stat: value}; adds ``frac_above`` when ``threshold`` is set."""
    profile = np.asarray(profile, dtype=float)
    funcs = {"mean": np.mean, "max": np.max, "min": np.min, "std": np.std}
    out = {f"{prefix}_{s}": float(funcs[s](profile)) for s in stats}
    if threshold is not None:
        out[f"{prefix}_frac_above"] = float(np.mean(profile > threshold))
    return out


def _gapped_pair_counts(sequence: str, gap: int, alphabet: dict[str, str] | None = None) -> Counter:
    """Count pairs (x_i, x_{i+gap+1}), optionally mapped through a residue -> group alphabet."""
    seq = [alphabet[aa] for aa in sequence] if alphabet else list(sequence)
    step = gap + 1
    return Counter((seq[i], seq[i + step]) for i in range(len(seq) - step))


def _kmer_counts(sequence: str, k: int, alphabet: dict[str, str] | None = None) -> Counter:
    """Overlapping k-mer counts (as tuples) over the raw or reduced alphabet."""
    seq = [alphabet[aa] for aa in sequence] if alphabet else list(sequence)
    return Counter(tuple(seq[i:i + k]) for i in range(len(seq) - k + 1))


def _charge_vector(sequence: str, his: bool = False) -> np.ndarray:
    """+1 for K/R (and H when ``his``), -1 for D/E, else 0."""
    positive = "KRH" if his else "KR"
    return np.array([1 if aa in positive else -1 if aa in "DE" else 0 for aa in sequence], dtype=float)


def _runs(mask: np.ndarray, min_len: int = 1) -> list[tuple[int, int]]:
    """Return (start, end) half-open index pairs of consecutive True runs of length >= min_len."""
    mask = np.concatenate(([False], np.asarray(mask, dtype=bool), [False]))
    edges = np.flatnonzero(np.diff(mask.astype(int)))
    return [(int(s), int(e)) for s, e in zip(edges[::2], edges[1::2]) if e - s >= min_len]


def _to_df(values: list, columns: list[str]) -> pd.DataFrame:
    """Round floats to 3 dp and return a 1 x N DataFrame."""
    values = [round(float(v), 3) if isinstance(v, (float, np.floating)) else
              int(v) if isinstance(v, (bool, np.bool_, np.integer)) else v for v in values]
    return pd.DataFrame([values], columns=columns)
