################################################################################
#############                      Encode                         ##############
################################################################################

from __future__ import annotations

import os
import csv
import pandas as pd
from tqdm import tqdm as _tqdm

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
from . import (composition as _composition, autocorrelation as _autocorrelation, ctd as _ctd,
               conjoint_triad as _conjoint_triad_module, sequence_order as _sequence_order, scales as _scales,
               physicochemical as _physicochemical, complexity as _complexity, signal as _signal,
               structural as _structural)
from .utils import _DESCRIPTOR_REGISTRY

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
}) | frozenset(
    #v1.5.0 descriptors outside the classic set that take parameters
    name for name, meta in _DESCRIPTOR_REGISTRY.items()
    if not meta.get("in_default_set", True) and set(meta["parameters"]) - {"sequence"}
)

#descriptor sets accepted by all_descriptors / encode_dataset
_DESCRIPTOR_SETS = ("classic", "extended", "all")

#modules holding each registry "module" value
_MODULES = {
    "composition": _composition, "autocorrelation": _autocorrelation, "ctd": _ctd,
    "conjoint_triad": _conjoint_triad_module, "sequence_order": _sequence_order, "scales": _scales,
    "physicochemical": _physicochemical, "complexity": _complexity, "signal": _signal,
    "structural": _structural,
}


def _set_members(descriptor_set: str) -> list[str]:
    """Registry descriptors (outside the classic set) that belong to ``descriptor_set``, in registry order."""
    if descriptor_set == "classic":
        return []
    return [name for name, meta in _DESCRIPTOR_REGISTRY.items()
            if not meta.get("in_default_set", True) and descriptor_set in meta.get("descriptor_sets", [])]


def all_descriptors(
    sequence: str,
    params: dict[str, dict] | None = None,
    show_progress: bool = False,
    output_file: str | None = None,
    save_each_sequence: bool = True,
    descriptor_set: str = "classic",
    include: list[str] | None = None,
    exclude: list[str] | None = None,
    _show_descriptor_list: bool = True,
) -> pd.DataFrame:
    """
    Calculate all available protpy descriptors for a protein sequence and
    return them concatenated into a single DataFrame (1 x N_total).

    ``descriptor_set`` chooses which descriptors are computed:

    - ``"classic"`` (default): the v1.4 descriptor set described below, with an
      unchanged output.
    - ``"extended"``: classic plus every v1.5 descriptor whose default output
      has at most 150 columns (factor scales, charge patterning,
      hydrophobicity/stability, complexity, signal and structural descriptors).
    - ``"all"``: extended plus the high-dimensional v1.5 descriptors (CKSAAP,
      DDE, ASDC, distance pair and auto/cross covariance).

    Set membership comes from the descriptor registry (``descriptor_sets``).
    ``include`` adds further registered descriptors by name and ``exclude``
    removes descriptors by name (classic ones included, e.g. "ctd_" or
    "tripeptide_composition"). Some extended descriptors are NaN where their
    method leaves them undefined (``plaac`` LLR/COREscore for short or
    non-prion sequences, ``half_life`` E. coli for an N-terminal Pro).

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

        plus any v1.5 descriptor with parameters (see
        ``protpy.get_descriptor_info(name)["parameters"]``), e.g.
        ``{"kappa": {"blob_sizes": (5,)}, "aggrescan": {"window": 7}}``.

        Example::

            df = protpy.all_descriptors(seq, params={
                'charge_distribution':            {'ph': 5.5},
                'kmer_composition':               {'k': 3},
                'pseudo_amino_acid_composition':  {'lamda': 10, 'weight': 0.1},
                'moreaubroto_autocorrelation':    {'lag': 15, 'normalize': False},
                'quasi_sequence_order_all':       {'lag': 20, 'weight': 0.2},
            })

    :show_progress: bool (default=False)
        When ``True``, display a ``tqdm`` progress bar that names the
        descriptor currently being computed and tracks overall progress.

    :output_file: str | None (default=None)
        Path to the CSV file where results are written.  ``None`` (default)
        disables all file output.  When calling this function in a loop over
        multiple sequences the rows are accumulated in the file across calls;
        each call appends one row to any previously saved sequences.  A
        ``ValueError`` is raised, and the file left untouched, if an existing
        file has different columns.

    :save_each_sequence: bool (default=True)
        When ``True``, the output CSV is updated once after all descriptors
        for the current sequence have been computed, appending a complete
        (no-NaN) row.  When calling this function in a loop over multiple
        sequences, each completed sequence is written immediately so that
        progress is preserved if the process is interrupted.  When ``False``,
        nothing is written to the CSV within this call.  Has no effect when
        ``output_file`` is ``None``.

    :descriptor_set: str (default="classic")
        "classic", "extended" or "all" (see above).

    :include: list[str] | None (default=None)
        extra registered descriptor names to compute.

    :exclude: list[str] | None (default=None)
        descriptor names to skip.

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
    >>> df = protpy.all_descriptors(seq, show_progress=True)  # with progress bar
    >>> df = protpy.all_descriptors(seq, output_file="out.csv")  # save to file
    >>> df = protpy.all_descriptors(seq, output_file=None)  # no file output
    """
    sequence = _validate_sequence(sequence)

    if params is None:
        params = {}

    if not isinstance(params, dict):
        raise TypeError(f"params must be a dict, got {type(params)}.")

    if descriptor_set not in _DESCRIPTOR_SETS:
        raise ValueError(f"Invalid descriptor_set '{descriptor_set}'. Valid sets are: {list(_DESCRIPTOR_SETS)}.")
    include, exclude = list(include or []), set(exclude or [])
    unregistered = [n for n in include + sorted(exclude) if n not in _DESCRIPTOR_REGISTRY]
    if unregistered:
        raise ValueError(f"Unknown descriptor name(s) in include/exclude: {unregistered}. "
                         "See protpy.list_descriptors().")

    unknown = set(params) - _PARAMETERISED_DESCRIPTORS
    if unknown:
        raise ValueError(
            f"Unknown or non-parameterisable descriptor(s) in params: {sorted(unknown)}. "
            f"Valid keys are: {sorted(_PARAMETERISED_DESCRIPTORS)}."
        )

    def _p(name: str) -> dict:
        return params.get(name, {})

    # k=2 duplicates dipeptide_composition; k=3 duplicates tripeptide_composition.
    # Only include kmer_composition when the user requests a k value that is unique.
    _kmer_k = _p("kmer_composition").get("k", 2)

    steps: list[tuple[str, str, object]] = [
        # ── Composition ────────────────────────────────────────────────────────
        ("[Composition]    Amino acid composition",                    "amino_acid_composition",                     lambda: amino_acid_composition(sequence)),
        ("[Composition]    Dipeptide composition",                     "dipeptide_composition",                      lambda: dipeptide_composition(sequence)),
        ("[Composition]    Tripeptide composition",                    "tripeptide_composition",                     lambda: tripeptide_composition(sequence)),
        ("[Composition]    GRAVY index",                               "gravy",                                      lambda: gravy(sequence)),
        ("[Composition]    Aromaticity",                               "aromaticity",                                lambda: aromaticity(sequence)),
        ("[Composition]    Instability index",                         "instability_index",                          lambda: instability_index(sequence)),
        ("[Composition]    Isoelectric point",                         "isoelectric_point",                          lambda: isoelectric_point(sequence)),
        ("[Composition]    Molecular weight",                          "molecular_weight",                           lambda: molecular_weight(sequence)),
        ("[Composition]    Charge distribution",                       "charge_distribution",                        lambda: charge_distribution(sequence, **_p("charge_distribution"))),
        ("[Composition]    Hydrophobic / polar / charged comp.",       "hydrophobic_polar_charged_composition",      lambda: hydrophobic_polar_charged_composition(sequence)),
        ("[Composition]    Secondary structure propensity",            "secondary_structure_propensity",             lambda: secondary_structure_propensity(sequence)),
        *([("[Composition]    K-mer composition",                         "kmer_composition",                           lambda: kmer_composition(sequence, **_p("kmer_composition")))]
          if _kmer_k not in (2, 3) else []),
        ("[Composition]    Reduced alphabet composition",              "reduced_alphabet_composition",               lambda: reduced_alphabet_composition(sequence, **_p("reduced_alphabet_composition"))),
        ("[Composition]    Motif composition",                         "motif_composition",                          lambda: motif_composition(sequence, **_p("motif_composition"))),
        ("[Composition]    Aliphatic index",                           "aliphatic_index",                            lambda: aliphatic_index(sequence)),
        ("[Composition]    Extinction coefficient",                    "extinction_coefficient",                     lambda: extinction_coefficient(sequence)),
        ("[Composition]    Boman index",                               "boman_index",                                lambda: boman_index(sequence)),
        ("[Composition]    Aggregation propensity",                    "aggregation_propensity",                     lambda: aggregation_propensity(sequence, **_p("aggregation_propensity"))),
        ("[Composition]    Hydrophobic moment",                        "hydrophobic_moment",                         lambda: hydrophobic_moment(sequence, **_p("hydrophobic_moment"))),
        ("[Composition]    Shannon entropy",                           "shannon_entropy",                            lambda: shannon_entropy(sequence)),
        ("[Composition]    Pseudo amino acid composition",             "pseudo_amino_acid_composition",              lambda: pseudo_amino_acid_composition(sequence, **_p("pseudo_amino_acid_composition"))),
        ("[Composition]    Amphiphilic pseudo amino acid comp.",       "amphiphilic_pseudo_amino_acid_composition",  lambda: amphiphilic_pseudo_amino_acid_composition(sequence, **_p("amphiphilic_pseudo_amino_acid_composition"))),
        # ── Autocorrelation ────────────────────────────────────────────────────
        ("[Autocorrelation] Moreaubroto autocorrelation",              "moreaubroto_autocorrelation",                lambda: moreaubroto_autocorrelation(sequence, **_p("moreaubroto_autocorrelation"))),
        ("[Autocorrelation] Moran autocorrelation",                    "moran_autocorrelation",                      lambda: moran_autocorrelation(sequence, **_p("moran_autocorrelation"))),
        ("[Autocorrelation] Geary autocorrelation",                    "geary_autocorrelation",                      lambda: geary_autocorrelation(sequence, **_p("geary_autocorrelation"))),
        # ── Conjoint Triad ─────────────────────────────────────────────────────
        ("[Conjoint Triad]  Conjoint triad",                           "conjoint_triad",                             lambda: _conjoint_triad(sequence)),
        # ── CTD ────────────────────────────────────────────────────────────────
        ("[CTD]             Composition, transition & distribution",   "ctd_",                                       lambda: ctd_(sequence, **{"all_ctd": True, **_p("ctd_")})),
        # ── Sequence Order ─────────────────────────────────────────────────────
        ("[Sequence Order]  Sequence order coupling number",           "sequence_order_coupling_number_all",         lambda: sequence_order_coupling_number_all(sequence, **_p("sequence_order_coupling_number_all"))),
        ("[Sequence Order]  Quasi sequence order",                     "quasi_sequence_order_all",                   lambda: quasi_sequence_order_all(sequence, **_p("quasi_sequence_order_all"))),
    ]

    #v1.5 descriptors from the chosen set and include, then apply exclude
    extra = _set_members(descriptor_set) + [n for n in include if n not in _set_members(descriptor_set)]
    classic_names = {name for _, name, _ in steps}
    for name in extra:
        if name in classic_names:
            continue
        meta = _DESCRIPTOR_REGISTRY[name]
        func = getattr(_MODULES[meta["module"]], name)
        steps.append((f"[{meta['category']}] {name}", name,
                      lambda func=func, name=name: func(sequence, **_p(name))))
    steps = [step for step in steps if step[1] not in exclude]

    descriptor_names = [label.split("]", 1)[1].strip() for label, _, _ in steps]
    if _show_descriptor_list:
        print("Calculating descriptors: " + ", ".join(descriptor_names))

    frames: list[pd.DataFrame] = []

    if show_progress:
        bar = _tqdm(
            total=len(steps),
            unit="descriptor",
            dynamic_ncols=True,
            bar_format="{desc} {percentage:3.0f}%|{bar:28}| {n:02d}/{total} [{elapsed}<{remaining}]",
        )
        for label, _, fn in steps:
            bar.set_description(label)
            frames.append(fn())
            bar.update(1)
        bar.set_description(f"All {len(steps)} descriptors computed")
        bar.close()
    else:
        for _, _, fn in steps:
            frames.append(fn())

    result = pd.concat(frames, axis=1)

    if output_file is not None and save_each_sequence:
        _prepare_output_file(output_file, list(result.columns))
        result.to_csv(output_file, mode="a", header=False, index=False)

    return result


def encode_dataset(
    sequences: "list[str] | pd.DataFrame",
    seq_col: str = "sequence",
    params: dict[str, dict] | None = None,
    show_progress: bool = False,
    output_file: str | None = "descriptors.csv",
    save_each_sequence: bool = True,
    strip_gaps: bool = False,
    verbose: bool = True,
    descriptor_set: str = "classic",
    include: list[str] | None = None,
    exclude: list[str] | None = None,
) -> pd.DataFrame:
    """
    Compute all available descriptors for every sequence in a dataset and
    return a DataFrame whose **first column is a unique sequence identifier**.

    Parameters
    ==========
    :sequences: list[str] | pd.DataFrame
        Input sequences in one of two forms:

        - ``list[str]`` — a plain list of protein sequences.  Identifiers are
          generated automatically as ``seq1, seq2, seq3, …`` and the
          identifier column is named ``"sequence_id"``.
        - ``pd.DataFrame`` — typically the output of
          ``read_sequences_from_file(..., id_col=…)``, which has two columns
          ``[id_col, seq_col]``.  The first column is used as the identifier;
          ``seq_col`` selects the sequence data.

    :seq_col: str (default="sequence")
        Name of the column containing sequences when *sequences* is a
        ``pd.DataFrame``.

    :params: dict[str, dict] | None (default=None)
        Per-descriptor keyword arguments forwarded to each ``all_descriptors``
        call.  See ``all_descriptors`` for the full list of accepted keys.

    :show_progress: bool (default=False)
        When ``True``, a ``tqdm`` progress bar is displayed while computing
        descriptors for each individual sequence.

    :output_file: str | None (default="descriptors.csv")
        Path to the CSV file where results are written.  Each completed row
        (identifier prepended) is appended immediately after that sequence
        finishes, so progress is preserved if the process is interrupted.
        Rows already present in the file from a previous run are retained,
        except rows whose identifier is encoded again in this run, which are
        replaced.  A ``ValueError`` is raised, and the file left untouched, if
        an existing file has different columns.  Pass ``None`` to disable all
        file output.

    :save_each_sequence: bool (default=True)
        When ``True``, the CSV is updated after each sequence completes.
        When ``False``, the file is never written inside this call.  Has no
        effect when ``output_file`` is ``None``.

    :strip_gaps: bool (default=False)
        When ``True``, gap characters (``-``) are removed from each sequence
        before encoding.  Useful for sequences derived from multiple sequence
        alignments (MSAs) where ``-`` denotes an alignment gap rather than an
        amino acid.  When ``False`` (default), a ``-`` in any sequence raises
        a ``ValueError``.

    :verbose: bool (default=True)
        When ``True``, print a status line for each sequence as it is encoded
        (``"Encoding sequence N/total: id"``).  Set to ``False`` to suppress
        all stdout output from this function.

    :descriptor_set: str (default="classic")
        "classic", "extended" or "all"; see ``all_descriptors``.

    :include: list[str] | None (default=None)
        extra descriptor names; see ``all_descriptors``.

    :exclude: list[str] | None (default=None)
        descriptor names to skip; see ``all_descriptors``.

    Returns
    =======
    :result: pd.DataFrame
        Shape ``(n_sequences, 1 + n_features)`` where the first column holds
        the sequence identifier and the remaining columns are descriptor
        values.  Returns an empty DataFrame when *sequences* is empty.

    Raises
    ======
    TypeError
        If *sequences* is neither a ``list`` nor a ``pd.DataFrame``.
    ValueError
        If *sequences* is a ``pd.DataFrame`` and *seq_col* is not present.

    Examples
    ========
    >>> import protpy
    >>> # Plain list — identifiers auto-generated as seq1, seq2, …
    >>> seqs = ["MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGD",
    ...         "ACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWY"]
    >>> df = protpy.encode_dataset(seqs)
    >>> df.columns[0]
    'sequence_id'
    >>> df["sequence_id"].tolist()
    ['seq1', 'seq2']

    >>> # DataFrame with named identifiers (e.g. from thermostability.txt)
    >>> data = protpy.read_sequences_from_file(
    ...     "thermostability.txt", seq_col="sequence", id_col="name"
    ... )
    >>> df = protpy.encode_dataset(data, seq_col="sequence")
    >>> df.columns[0]
    'name'
    """
    if isinstance(sequences, list):
        id_col_name = "sequence_id"
        ids = [f"seq{i + 1}" for i in range(len(sequences))]
        seqs = [str(s) for s in sequences]
    elif isinstance(sequences, pd.DataFrame):
        if sequences.empty:
            return pd.DataFrame()
        id_col_name = sequences.columns[0]
        if seq_col not in sequences.columns:
            raise ValueError(
                f"seq_col '{seq_col}' not found in DataFrame. "
                f"Available columns: {list(sequences.columns)}."
            )
        ids = sequences[id_col_name].tolist()
        seqs = sequences[seq_col].tolist()
    else:
        raise TypeError(
            f"sequences must be a list[str] or pd.DataFrame, got {type(sequences)}."
        )

    completed: list[pd.DataFrame] = []
    total = len(seqs)

    for idx, (seq_id, seq) in enumerate(zip(ids, seqs), start=1):
        if verbose:
            print(f"Encoding sequence {idx}/{total}: {seq_id}")
        if strip_gaps:
            seq = seq.replace("-", "")
        elif "-" in seq:
            raise ValueError(
                f"Sequence '{seq_id}' contains gap characters ('-'). "
                "If the sequences are from a multiple sequence alignment (MSA), "
                "set strip_gaps=True to remove them before encoding."
            )
        row_df = all_descriptors(
            seq,
            params=params,
            show_progress=show_progress,
            output_file=None,
            save_each_sequence=False,
            descriptor_set=descriptor_set,
            include=include,
            exclude=exclude,
            _show_descriptor_list=(idx == 1 and verbose),
        )
        row_df.insert(0, id_col_name, seq_id)
        completed.append(row_df)

        if output_file is not None and save_each_sequence:
            if idx == 1:
                _prepare_output_file(output_file, list(row_df.columns), replace_ids=ids)
            row_df.to_csv(output_file, mode="a", header=False, index=False)

    if not completed:
        return pd.DataFrame()
    return pd.concat(completed, ignore_index=True)


def _prepare_output_file(output_file: str, columns: list, replace_ids: list | None = None) -> None:
    """
    Get ``output_file`` ready for rows to be appended. A missing or empty file
    gets a header row. An existing file must have exactly ``columns`` as its
    header, otherwise a ValueError is raised and the file is not modified.
    When ``replace_ids`` is given, existing rows whose first column matches
    one of them are removed so a re-run does not duplicate them.
    """
    if not os.path.isfile(output_file) or os.path.getsize(output_file) == 0:
        with open(output_file, "w", newline="") as fh:
            csv.writer(fh).writerow(columns)
        return

    with open(output_file, newline="") as fh:
        rows = list(csv.reader(fh))
    if rows[0] != [str(c) for c in columns]:
        raise ValueError(
            f"Existing output file '{output_file}' has different columns from this run; "
            "refusing to overwrite it. Pass a different output_file or remove the file."
        )
    if replace_ids:
        replace = {str(i) for i in replace_ids}
        kept = [rows[0]] + [r for r in rows[1:] if r and r[0] not in replace]
        if len(kept) != len(rows):
            with open(output_file, "w", newline="") as fh:
                csv.writer(fh).writerows(kept)


def _parse_fasta(filepath: str) -> list[str]:
    sequences: list[str] = []
    current: list[str] = []
    with open(filepath) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current:
                    sequences.append("".join(current).upper())
                    current = []
            else:
                current.append(line)
    if current:
        sequences.append("".join(current).upper())
    return sequences


def _parse_txt(filepath: str) -> list[str]:
    sequences: list[str] = []
    with open(filepath) as fh:
        for line in fh:
            seq = line.strip().upper()
            if seq:
                sequences.append(seq)
    return sequences


def _detect_sep(filepath: str) -> str:
    with open(filepath) as fh:
        first_line = fh.readline()
    return "\t" if "\t" in first_line else ","


def _resolve_col(
    df: pd.DataFrame,
    col: str | int | None,
    param_name: str,
) -> str:
    if col is None:
        return df.columns[0]
    if isinstance(col, int):
        try:
            return df.columns[col]
        except IndexError:
            raise ValueError(
                f"{param_name} index {col} is out of range for a file with "
                f"{df.shape[1]} column(s)."
            )
    if isinstance(col, str):
        if col not in df.columns:
            raise ValueError(
                f"{param_name} '{col}' not found. Available columns: {list(df.columns)}."
            )
        return col
    raise TypeError(f"{param_name} must be a str, int, or None, got {type(col)}.")


def _parse_delimited(
    filepath: str,
    sep: str,
    seq_col: str | int | None,
    id_col: str | int | None,
    has_header: bool,
) -> "list[str] | pd.DataFrame":
    df = pd.read_csv(filepath, sep=sep, header=0 if has_header else None)

    seq_label = _resolve_col(df, seq_col, "seq_col")

    if id_col is None:
        raw = df[seq_label]
        return [str(s).strip().upper() for s in raw if pd.notna(s) and str(s).strip()]

    id_label = _resolve_col(df, id_col, "id_col")
    result = df[[id_label, seq_label]].copy()
    result = result[result[seq_label].notna()].copy()
    result[seq_label] = result[seq_label].astype(str).str.strip().str.upper()
    result = result[result[seq_label].str.len() > 0].reset_index(drop=True)
    return result


def read_sequences_from_file(
    filepath: str | os.PathLike,
    seq_col: str | int | None = None,
    id_col: str | int | None = None,
    has_header: bool = True,
) -> "list[str] | pd.DataFrame":
    """
    Parse protein sequences from a FASTA, TXT, TSV, or CSV file.

    The file format is inferred from the extension:

    - ``.fasta`` / ``.fa`` — standard FASTA format; multi-line sequences are
      joined automatically; every record in the file is returned as a
      ``list[str]``.
    - ``.txt`` without ``seq_col`` — one raw sequence per non-blank line,
      returned as a ``list[str]``.
    - ``.txt`` with ``seq_col`` — treated as a delimited file; the separator
      (``,`` or ``\\t``) is auto-detected from the first line.
    - ``.csv`` — comma-separated values.
    - ``.tsv`` — tab-separated values.

    When ``id_col`` is provided the function returns a ``pd.DataFrame``
    containing only the identifier and sequence columns; all other columns
    are dropped and rows with blank sequences are removed.  When ``id_col``
    is omitted the function returns a plain ``list[str]`` of sequences.

    Parameters
    ==========
    :filepath: str | os.PathLike
        Path to the input file.

    :seq_col: str | int | None (default=None)
        The column that contains sequences.  Pass a column *name* (str) to
        select by header, or a zero-based *index* (int) to select by
        position.  ``None`` selects the first column.  Ignored for plain
        FASTA and plain TXT files.

    :id_col: str | int | None (default=None)
        The column that identifies each sequence row (e.g. ``"name"``).
        When provided, all other columns are dropped and a ``pd.DataFrame``
        with columns ``[id_col, seq_col]`` is returned.  ``None`` returns a
        ``list[str]`` of sequences only.  Ignored for FASTA files.

    :has_header: bool (default=True)
        Whether the first row of a delimited file is a header.  Set to
        ``False`` for headerless files (column names are then integers).
        Ignored for FASTA and plain TXT files.

    Returns
    =======
    :result: list[str] | pd.DataFrame
        - ``list[str]`` when ``id_col`` is ``None``: uppercase sequences,
          blank entries excluded.
        - ``pd.DataFrame`` when ``id_col`` is given: two-column frame with
          the identifier and (uppercased) sequence columns; rows with blank
          sequences are excluded and the index is reset.

    Raises
    ======
    FileNotFoundError
        If the file at ``filepath`` does not exist.
    ValueError
        If the extension is unsupported, or if a requested column is not
        found in a delimited file.
    TypeError
        If ``seq_col`` or ``id_col`` is not a str, int, or None.

    Examples
    ========
    >>> import protpy
    >>> # FASTA / plain TXT
    >>> seqs = protpy.read_sequences_from_file("proteins.fasta")
    >>> seqs = protpy.read_sequences_from_file("seqs.txt")
    >>> # CSV / TSV — sequences only
    >>> seqs = protpy.read_sequences_from_file("data.csv", seq_col="sequence")
    >>> seqs = protpy.read_sequences_from_file("data.tsv", seq_col=1)
    >>> # Multi-column TXT/CSV/TSV — keep id + sequence, drop everything else
    >>> df = protpy.read_sequences_from_file(
    ...     "thermostability.txt", seq_col="sequence", id_col="name"
    ... )
    """
    filepath = str(filepath)

    if not os.path.isfile(filepath):
        raise FileNotFoundError(f"No such file: '{filepath}'")

    ext = os.path.splitext(filepath)[1].lower()

    if ext in (".fasta", ".fa"):
        return _parse_fasta(filepath)
    elif ext == ".txt":
        if seq_col is not None:
            sep = _detect_sep(filepath)
            return _parse_delimited(filepath, sep=sep, seq_col=seq_col, id_col=id_col, has_header=has_header)
        return _parse_txt(filepath)
    elif ext == ".csv":
        return _parse_delimited(filepath, sep=",", seq_col=seq_col, id_col=id_col, has_header=has_header)
    elif ext == ".tsv":
        return _parse_delimited(filepath, sep="\t", seq_col=seq_col, id_col=id_col, has_header=has_header)
    else:
        raise ValueError(
            f"Unsupported file extension '{ext}'. "
            "Expected one of: .fasta, .fa, .txt, .csv, .tsv."
        )