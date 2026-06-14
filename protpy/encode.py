################################################################################
#############                      Encode                         ##############
################################################################################

from __future__ import annotations

import pandas as pd

from .composition import (
    amino_acid_composition,
    dipeptide_composition,
    tripeptide_composition,
    gravy,
    aromaticity,
    instability_index,
    isoelectric_point,
    molecular_weight,
    charge_distribution,
    hydrophobic_polar_charged_composition,
    secondary_structure_propensity,
    kmer_composition,
    reduced_alphabet_composition,
    motif_composition,
    amino_acid_pair_composition,
    aliphatic_index,
    extinction_coefficient,
    boman_index,
    aggregation_propensity,
    hydrophobic_moment,
    shannon_entropy,
    pseudo_amino_acid_composition,
    amphiphilic_pseudo_amino_acid_composition,
)
from .autocorrelation import (
    moreaubroto_autocorrelation,
    moran_autocorrelation,
    geary_autocorrelation,
)
from .conjoint_triad import conjoint_triad as _conjoint_triad
from .ctd import ctd_
from .sequence_order import (
    sequence_order_coupling_number_all,
    quasi_sequence_order_all,
)
from ._constants import _validate_sequence

#descriptors that accept extra keyword arguments beyond sequence
_PARAMETERISED_DESCRIPTORS: frozenset[str] = frozenset({
    "charge_distribution",
    "kmer_composition",
    "reduced_alphabet_composition",
    "motif_composition",
    "aggregation_propensity",
    "hydrophobic_moment",
    "pseudo_amino_acid_composition",
    "amphiphilic_pseudo_amino_acid_composition",
    "moreaubroto_autocorrelation",
    "moran_autocorrelation",
    "geary_autocorrelation",
    "ctd_",
    "sequence_order_coupling_number_all",
    "quasi_sequence_order_all",
})


def all_descriptors(sequence: str, params: dict[str, dict] | None = None) -> pd.DataFrame:
    """
    Calculate all available protpy descriptors for a protein sequence and
    return them concatenated into a single DataFrame (1 x N_total).

    The combined variants of multi-property descriptors are used to avoid
    redundancy: ``ctd_`` is called with ``all_ctd=True`` (covers
    ctd_composition, ctd_transition and ctd_distribution across all 7
    properties), and ``sequence_order_coupling_number_all`` /
    ``quasi_sequence_order_all`` cover both the schneider-wrede and grantham
    distance matrices.

    Parameters
    ==========
    :sequence: str
        Protein sequence string (case-insensitive, canonical amino acids only).
        Must be longer than 30 residues to satisfy the default lag requirements
        of autocorrelation and sequence-order descriptors (unless overridden
        via ``params``).

    :params: dict[str, dict] | None (default=None)
        Optional per-descriptor keyword arguments. Keys must be descriptor
        function names; values are dicts of keyword arguments passed directly
        to that descriptor. Unrecognised descriptor names raise a ``ValueError``.
        Descriptors not present in ``params`` use their default parameter values.

        Only the following descriptors accept extra parameters::

            charge_distribution           → ph
            kmer_composition              → k
            reduced_alphabet_composition  → alphabet_size
            motif_composition             → motifs
            aggregation_propensity        → window, hydrophobicity_threshold,
                                            charge_threshold
            hydrophobic_moment            → window, angle
            pseudo_amino_acid_composition → lamda, weight, properties
            amphiphilic_pseudo_amino_acid_composition
                                          → lamda, weight, properties
            moreaubroto_autocorrelation   → lag, properties, normalize
            moran_autocorrelation         → lag, properties, normalize
            geary_autocorrelation         → lag, properties, normalize
            ctd_                          → all_ctd
            sequence_order_coupling_number_all → lag
            quasi_sequence_order_all      → lag, weight

        Example::

            df = protpy.all_descriptors(seq, params={
                'charge_distribution':            {'ph': 5.5},
                'kmer_composition':               {'k': 3},
                'pseudo_amino_acid_composition':  {'lamda': 10, 'weight': 0.1},
                'moreaubroto_autocorrelation':    {'lag': 15, 'normalize': False},
                'quasi_sequence_order_all':       {'lag': 20, 'weight': 0.2},
            })

    Returns
    =======
    :descriptors_df: pd.DataFrame
        DataFrame of shape 1 x N containing all descriptor values, where N is
        the total number of features across all descriptors.

    Raises
    ======
    TypeError
        If ``sequence`` is not a string, or if ``params`` is not a dict.
    ValueError
        If ``sequence`` is empty or contains non-canonical amino acid characters,
        or if a key in ``params`` is not a recognised descriptor name.

    Examples
    ========
    >>> import protpy
    >>> seq = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGD"
    >>> df = protpy.all_descriptors(seq)
    >>> df.shape[0]
    1
    """
    sequence = _validate_sequence(sequence)

    if params is None:
        params = {}

    if not isinstance(params, dict):
        raise TypeError(f"params must be a dict, got {type(params)}.")

    unknown = set(params) - _PARAMETERISED_DESCRIPTORS
    if unknown:
        raise ValueError(
            f"Unknown or non-parameterisable descriptor(s) in params: {sorted(unknown)}. "
            f"Valid keys are: {sorted(_PARAMETERISED_DESCRIPTORS)}."
        )

    def _p(name: str) -> dict:
        return params.get(name, {})

    frames: list[pd.DataFrame] = []

    # ── Composition ────────────────────────────────────────────────────────────
    frames.append(amino_acid_composition(sequence))
    frames.append(dipeptide_composition(sequence))
    frames.append(tripeptide_composition(sequence))
    frames.append(gravy(sequence))
    frames.append(aromaticity(sequence))
    frames.append(instability_index(sequence))
    frames.append(isoelectric_point(sequence))
    frames.append(molecular_weight(sequence))
    frames.append(charge_distribution(sequence, **_p("charge_distribution")))
    frames.append(hydrophobic_polar_charged_composition(sequence))
    frames.append(secondary_structure_propensity(sequence))
    frames.append(kmer_composition(sequence, **_p("kmer_composition")))
    frames.append(reduced_alphabet_composition(sequence, **_p("reduced_alphabet_composition")))
    frames.append(motif_composition(sequence, **_p("motif_composition")))
    frames.append(amino_acid_pair_composition(sequence))
    frames.append(aliphatic_index(sequence))
    frames.append(extinction_coefficient(sequence))
    frames.append(boman_index(sequence))
    frames.append(aggregation_propensity(sequence, **_p("aggregation_propensity")))
    frames.append(hydrophobic_moment(sequence, **_p("hydrophobic_moment")))
    frames.append(shannon_entropy(sequence))
    frames.append(pseudo_amino_acid_composition(sequence, **_p("pseudo_amino_acid_composition")))
    frames.append(amphiphilic_pseudo_amino_acid_composition(sequence, **_p("amphiphilic_pseudo_amino_acid_composition")))

    # ── Autocorrelation ────────────────────────────────────────────────────────
    frames.append(moreaubroto_autocorrelation(sequence, **_p("moreaubroto_autocorrelation")))
    frames.append(moran_autocorrelation(sequence, **_p("moran_autocorrelation")))
    frames.append(geary_autocorrelation(sequence, **_p("geary_autocorrelation")))

    # ── Conjoint Triad ─────────────────────────────────────────────────────────
    frames.append(_conjoint_triad(sequence))

    # ── CTD (Composition + Transition + Distribution across all 7 properties) ──
    frames.append(ctd_(sequence, all_ctd=True, **_p("ctd_")))

    # ── Sequence Order ─────────────────────────────────────────────────────────
    frames.append(sequence_order_coupling_number_all(sequence, **_p("sequence_order_coupling_number_all")))
    frames.append(quasi_sequence_order_all(sequence, **_p("quasi_sequence_order_all")))

    return pd.concat(frames, axis=1)