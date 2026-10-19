from __future__ import annotations

#canonical 20 amino acids used throughout the protpy package
amino_acids: list[str] = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P",
    "Q", "R", "S", "T", "V", "W", "Y"]


def _validate_sequence(sequence: str) -> str:
    """
    Validate a protein sequence string.

    Checks that the input is a str, uppercases it, and confirms every
    character is one of the 20 canonical amino acids.

    Parameters
    ==========
    :sequence: str
        raw protein sequence (case-insensitive).

    Returns
    =======
    :sequence: str
        uppercased, validated protein sequence.

    Raises
    ======
    TypeError
        if ``sequence`` is not a str.
    ValueError
        if ``sequence`` is empty or contains a non-canonical amino acid character.
    """
    if not isinstance(sequence, str):
        raise TypeError(
            f"Input sequence must be a string, got input of type {type(sequence)}."
        )
    sequence = sequence.upper()
    #raise error on empty sequence before the aa loop
    if not sequence:
        raise ValueError("Input sequence cannot be empty.")
    for aa in sequence:
        if aa not in amino_acids:
            raise ValueError(
                f"Invalid amino acid found in protein sequence: {aa}."
            )
    return sequence


def _validate_int(name: str, value, minimum: int = 0) -> int:
    """
    Parse ``value`` as an int and check it is >= ``minimum``.

    Raises
    ======
    ValueError
        if ``value`` cannot be parsed as an int or is below ``minimum``.
    """
    if isinstance(value, bool):
        raise ValueError(f"Invalid {name} value input: {value}.")
    try:
        value = int(value)
    except (ValueError, TypeError):
        raise ValueError(f"Invalid {name} value input: {value}.")
    if value < minimum:
        raise ValueError(f"{name} must be >= {minimum}, got {value}.")
    return value


def _validate_weight(weight) -> float:
    """
    Check that ``weight`` is a non-negative number.

    Raises
    ======
    ValueError
        if ``weight`` is not a number or is negative.
    """
    if isinstance(weight, bool) or not isinstance(weight, (int, float)) or weight < 0:
        raise ValueError(f"weight must be a non-negative number, got {weight}.")
    return weight
