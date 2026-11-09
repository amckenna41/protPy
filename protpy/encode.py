################################################################################
#############                      Encode                         ##############
################################################################################

from __future__ import annotations

import os
import csv
import inspect
import importlib
import warnings
import pandas as pd
from tqdm import tqdm as _tqdm

from ._constants import _validate_sequence
from .utils import _DESCRIPTOR_REGISTRY

#descriptor sets accepted by all_descriptors / encode_dataset
_DESCRIPTOR_SETS = ("classic", "extended", "all")

#arguments used when a descriptor is computed as part of a set (user params override them)
_BUNDLE_KWARGS: dict[str, dict] = {
    "quasi_sequence_order": {"distance_matrix": "both"},
}

#columns kept when a descriptor is computed as part of a set without user params
#(ExtCoeff_Reduced is ExtCoeff_Oxidized minus 125 per cystine, so one column is kept)
_BUNDLE_COLUMNS: dict[str, list[str]] = {
    "extinction_coefficient": ["ExtCoeff_Oxidized"],
}


#long sequence used to get the full column layout of a descriptor that a short
#sequence cannot fill (the 20 amino acids repeated, 1,000 residues)
_REFERENCE_SEQUENCE = "ACDEFGHIKLMNPQRSTVWY" * 50


def _run_in_bundle(func, sequence: str, kwargs: dict) -> tuple[pd.DataFrame, bool]:
    """
    Run a descriptor as part of a descriptor set. When the sequence is too short
    for the descriptor's lag (lag/lamda must be less than the sequence length), the
    descriptor is computed with the largest lag the sequence allows and the
    columns of the higher lags are NaN, so every sequence gets the same columns.
    If the descriptor cannot be computed at all, every column is NaN. Errors that
    a long sequence also raises (invalid parameters) are re-raised.

    Returns
    =======
    :df, padded: tuple[pd.DataFrame, bool]
        descriptor values, and whether any column was filled with NaN.
    """
    try:
        return func(sequence, **kwargs), False
    except ValueError as error:
        try:
            columns = func(_REFERENCE_SEQUENCE, **kwargs).columns
        except ValueError:
            raise error from None
        #the reference sequence works, so the error came from the short sequence
        df = pd.DataFrame(index=[0])
        lag_param = next((p for p in ("lag", "lamda") if p in inspect.signature(func).parameters), None)
        if lag_param and len(sequence) > 1:
            try:
                df = func(sequence, **{**kwargs, lag_param: len(sequence) - 1})
            except ValueError:
                pass
        return df.reindex(columns=columns), True


def _csv_safe(value) -> str:
    """
    Make a value safe to write as a CSV cell: text starting with =, +, -, @, tab
    or carriage return is prefixed with a single quote so spreadsheet programs do
    not run it as a formula (CSV injection). Numbers are written unchanged.
    """
    text = str(value)
    if isinstance(value, str) and text[:1] in ("=", "+", "-", "@", "\t", "\r"):
        return "'" + text
    return text


def _descriptor_function(name: str):
    """
    Get the function of a registered descriptor from its module. The module is
    imported by name because ``protpy.conjoint_triad`` and ``protpy.ctd`` are
    rebound to the functions of the same name by the package namespace.
    """
    return getattr(importlib.import_module(f"protpy.{_DESCRIPTOR_REGISTRY[name]['module']}"), name)


def _set_members(descriptor_set: str) -> list[str]:
    """
    Get the descriptors of a descriptor set, in the order of the descriptor
    registry. Membership is taken from the descriptor_sets field of each registry
    entry; the sets are nested (classic within extended within all).

    Parameters
    ==========
    :descriptor_set: str
        descriptor set: "classic", "extended" or "all".

    Returns
    =======
    :members: list[str]
        names of the descriptors in the set.
    """
    return [name for name, meta in _DESCRIPTOR_REGISTRY.items() if descriptor_set in meta["descriptor_sets"]]


def all_descriptors(
    sequence: str,
    params: dict[str, dict] | None = None,
    show_progress: bool = False,
    output_file: str | None = None,
    save_each_sequence: bool = True,
    descriptor_set: str = "classic",
    include: list[str] | None = None,
    exclude: list[str] | None = None,
    decimals: int | None = None,
    _show_descriptor_list: bool = True,
) -> pd.DataFrame:
    """
    Calculate a set of protpy descriptors for a protein sequence and return them
    concatenated into a single DataFrame (1 x N_total).

    ``descriptor_set`` chooses which descriptors are computed; the sets come from
    the descriptor registry (``descriptor_sets``) and are nested:

    - ``"classic"`` (default, about 1,400 columns): amino acid and dipeptide
      composition, the ProtParam scalars (GRAVY, aromaticity, instability index,
      pI, MW, aliphatic index, extinction coefficient), Boman index, charge
      distribution, hydrophobic moment, Shannon entropy, motifs, PseAAC, APseAAC,
      Geary autocorrelation, conjoint triad, CTD (7 properties) and
      quasi-sequence-order with both distance matrices.
    - ``"extended"``: classic plus the low-dimensional descriptors (natural vector,
      atomic composition, CKSAAGP, charge fractions, kappa, SCD, hydrophobicity
      profile, thermostability indices, Wootton-Federhen and Lempel-Ziv
      complexity, spectral descriptors, AGGRESCAN, charge-hydropathy, PLAAC).
    - ``"all"``: extended plus the high-dimensional descriptors (CKSAAP, ASDC,
      distance pair, auto/cross covariance, Moreau-Broto and Moran
      autocorrelation, sequence order coupling numbers).

    ``include`` adds further registered descriptors by name, e.g.
    ``include=["tripeptide_composition"]`` (8,000 columns, in no set) or ``"dde"``,
    and ``exclude`` removes descriptors by name. Some extended descriptors are NaN
    where their method leaves them undefined (``plaac`` LLR/COREscore for short or
    non-prion sequences, ``kappa`` for a sequence without charged residues).

    Since v1.6.0 the classic set no longer contains tripeptide composition, the
    Moreau-Broto and Moran autocorrelations (which duplicate Geary), the sequence
    order coupling numbers (contained in QSO), or the deprecated descriptors.

    Parameters
    ==========
    :sequence: str
        Protein sequence string (case-insensitive, canonical amino acids only).
        A sequence that is not longer than a descriptor's lag (30 by default for
        the autocorrelation, PseAAC and sequence-order descriptors) is computed
        with the largest lag it allows (length - 1) and the columns of the higher
        lags are NaN, with a ``UserWarning``; the columns are always the same, so
        short peptides can be encoded with long proteins. Lower the lags via
        ``params`` to avoid the NaN columns.

    :params: dict[str, dict] | None (default=None)
        Optional per-descriptor keyword arguments. Keys must be registered
        descriptor names that take parameters; values are dicts of keyword
        arguments passed to that descriptor, e.g.
        ``{"geary_autocorrelation": {"lag": 10}, "kappa": {"blob_sizes": (5,)}}``.
        See ``protpy.get_descriptor_info(name)["parameters"]``. Parameters for a
        descriptor that is not computed (not in the set or ``include``, or in
        ``exclude``) raise a ``UserWarning``.

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
        for the current sequence have been computed.  When ``False``, nothing is
        written to the CSV within this call.  Has no effect when
        ``output_file`` is ``None``.

    :descriptor_set: str (default="classic")
        "classic", "extended" or "all" (see above).

    :include: list[str] | None (default=None)
        extra registered descriptor names to compute.

    :exclude: list[str] | None (default=None)
        descriptor names to skip.

    :decimals: int | None (default=None)
        round every value to this many decimal places. Descriptors are not
        rounded otherwise.

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
        or if a key in ``params``, ``include`` or ``exclude`` is not a
        recognised descriptor name.

    Examples
    ========
    >>> import protpy
    >>> seq = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGD"
    >>> df = protpy.all_descriptors(seq)
    >>> df.shape[0]
    1
    >>> df = protpy.all_descriptors(seq, show_progress=True)  # with progress bar
    >>> df = protpy.all_descriptors(seq, output_file="out.csv")  # save to file
    >>> df = protpy.all_descriptors(seq, descriptor_set="extended", decimals=3)
    """
    sequence = _validate_sequence(sequence)

    if params is None:
        params = {}

    if not isinstance(params, dict):
        raise TypeError(f"params must be a dict, got {type(params)}.")

    if descriptor_set not in _DESCRIPTOR_SETS:
        raise ValueError(f"Invalid descriptor_set '{descriptor_set}'. Valid sets are: {list(_DESCRIPTOR_SETS)}.")
    #a single name may be passed as a str
    include = [include] if isinstance(include, str) else list(include or [])
    exclude = {exclude} if isinstance(exclude, str) else set(exclude or [])
    unregistered = [n for n in include + sorted(exclude) if n not in _DESCRIPTOR_REGISTRY]
    if unregistered:
        raise ValueError(f"Unknown descriptor name(s) in include/exclude: {unregistered}. "
                         "See protpy.list_descriptors().")

    parameterised = {name for name, meta in _DESCRIPTOR_REGISTRY.items() if set(meta["parameters"]) - {"sequence"}}
    unknown = set(params) - parameterised
    if unknown:
        raise ValueError(
            f"Unknown or non-parameterisable descriptor(s) in params: {sorted(unknown)}. "
            f"Valid keys are: {sorted(parameterised)}."
        )

    members = _set_members(descriptor_set)
    names = [n for n in members + [n for n in include if n not in members] if n not in exclude]
    not_run = sorted(set(params) - set(names))
    if not_run:
        warnings.warn(f"params given for {not_run}, which are not computed with descriptor_set="
                      f"'{descriptor_set}' and these include/exclude; add them with include=[...].",
                      UserWarning, stacklevel=2)
    for name, kwargs in params.items():
        if not isinstance(kwargs, dict):
            raise TypeError(f"params['{name}'] must be a dict of keyword arguments, got {type(kwargs)}.")

    def _step(name: str):
        func = _descriptor_function(name)

        def run() -> pd.DataFrame:
            df, was_padded = _run_in_bundle(func, sequence, {**_BUNDLE_KWARGS.get(name, {}), **params.get(name, {})})
            if was_padded:
                padded.append(name)
            if name in _BUNDLE_COLUMNS and name not in params:
                df = df[_BUNDLE_COLUMNS[name]]
            return df
        return f"[{_DESCRIPTOR_REGISTRY[name]['category']}] {name}", name, run

    padded: list[str] = []
    steps = [_step(name) for name in names]

    if _show_descriptor_list:
        print("Calculating descriptors: " + ", ".join(names))

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

    result = pd.concat(frames, axis=1) if frames else pd.DataFrame(index=[0])
    if padded:
        warnings.warn(f"The sequence ({len(sequence)} residues) is too short for the lags of {padded}: they "
                      "were computed with lag = length - 1 and the higher-lag columns are NaN. Lower the lags "
                      "via params to avoid this.", UserWarning, stacklevel=2)
    if decimals is not None:
        result = result.round(decimals)

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
    decimals: int | None = None,
) -> pd.DataFrame:
    """
    Compute a set of descriptors for every sequence in a dataset and
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

    :decimals: int | None (default=None)
        round every value to this many decimal places; see ``all_descriptors``.

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
        seqs = list(sequences)
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
        if not isinstance(seq, str):
            raise TypeError(f"Sequence '{seq_id}' must be a string, got {type(seq).__name__}: {seq!r}.")
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
            decimals=decimals,
            _show_descriptor_list=(idx == 1 and verbose),
        )
        row_df.insert(0, id_col_name, seq_id)
        completed.append(row_df)

        if output_file is not None and save_each_sequence:
            if idx == 1:
                _prepare_output_file(output_file, list(row_df.columns), replace_ids=ids)
            #escape the identifier so a spreadsheet does not run it as a formula
            row_df.assign(**{id_col_name: _csv_safe(seq_id)}).to_csv(output_file, mode="a", header=False, index=False)

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
    columns = [_csv_safe(c) for c in columns]
    if not os.path.isfile(output_file) or os.path.getsize(output_file) == 0:
        with open(output_file, "w", newline="") as fh:
            csv.writer(fh).writerow(columns)
        return

    with open(output_file, newline="") as fh:
        rows = list(csv.reader(fh))
    if rows[0] != columns:
        raise ValueError(
            f"Existing output file '{output_file}' has different columns from this run; "
            "refusing to overwrite it. Pass a different output_file or remove the file."
        )
    if replace_ids:
        replace = {_csv_safe(i) for i in replace_ids}
        kept = [rows[0]] + [r for r in rows[1:] if r and r[0] not in replace]
        if len(kept) != len(rows):
            with open(output_file, "w", newline="") as fh:
                csv.writer(fh).writerows(kept)


def _parse_fasta(filepath: str) -> tuple[list[str], list[str]]:
    """
    Parse a FASTA file into record identifiers (the first word of each header)
    and sequences. Multi-line sequences are joined, a trailing stop codon ``*``
    is removed and records without a sequence are skipped.
    """
    ids: list[str] = []
    sequences: list[str] = []
    header, current = None, []

    def flush():
        seq = "".join(current).upper().rstrip("*")
        if seq:
            ids.append(header.split()[0] if header and header.split() else f"seq{len(ids) + 1}")
            sequences.append(seq)

    with open(filepath) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                flush()
                header, current = line[1:].strip(), []
            else:
                current.append(line)
    flush()
    return ids, sequences


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
      joined automatically and a trailing stop codon ``*`` is removed; every
      record in the file is returned as a ``list[str]``, or with ``id_col`` as
      a DataFrame of record identifiers (the first word of each header, e.g.
      ``sp|P59594|SPIKE_SARS``) and sequences.
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
        ``list[str]`` of sequences only.  For FASTA files any ``id_col`` returns
        the header identifiers, in a column named ``id_col`` if it is a str (else
        ``"id"``); the sequence column is named ``seq_col`` if it is a str (else
        ``"sequence"``).

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
    >>> # FASTA with the header identifiers, ready for encode_dataset
    >>> df = protpy.read_sequences_from_file("proteins.fasta", id_col="id")
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
        ids, sequences = _parse_fasta(filepath)
        if id_col is None:
            return sequences
        return pd.DataFrame({id_col if isinstance(id_col, str) else "id": ids,
                             seq_col if isinstance(seq_col, str) else "sequence": sequences})
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