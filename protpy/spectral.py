################################################################################
#############                 Signal Processing                   ##############
################################################################################

from __future__ import annotations

"""
Spectral descriptors of amino-acid property signals: binned FFT power spectrum
and DCT coefficients. numpy only. (This module was protpy/signal.py before
v1.6.0; it was renamed so it no longer shadows the standard library module.)

References
==========
[1] Cosic, I. (1994). Macromolecular bioactivity: is it resonant interaction
    between macromolecules? Theory and applications. IEEE Trans Biomed Eng,
    41(12), 1101-1114. https://doi.org/10.1109/10.335859
[2] Wang, Y., Ding, Y., Guo, F., Wei, L., & Tang, J. (2017). Improved detection
    of DNA-binding proteins via compression technology on PSSM information.
    PLoS ONE, 12(9), e0185587. https://doi.org/10.1371/journal.pone.0185587
    (DCT used as a sequence-length-invariant compression of a residue signal)
"""

import numpy as np
import pandas as pd

from ._constants import _validate_sequence, _validate_int
from ._helpers import _property_matrix, _to_df, _warn_deprecated

__all__ = ["spectral_descriptors",
    #deprecated, removed in protpy 2.0
    "fft_descriptors", "dct_descriptors"]


def _signal(sequence: str, scale: str | dict) -> tuple[np.ndarray, str]:
    """
    Convert a protein sequence into a numerical signal using a single-value amino
    acid scale, with the mean of the signal subtracted.

    Parameters
    ==========
    :sequence: str
        protein sequence.
    :scale: str | dict
        single-value scale name, AAindex accession number or dict of amino acid
        values.

    Returns
    =======
    :signal: np.ndarray
        mean-centred property values along the sequence.
    :name: str
        name of the scale used in column names ("custom" for a dict).
    """
    values = _property_matrix(sequence, scale)
    if values.shape[1] != 1:
        raise ValueError(f"Scale '{scale}' is multi-dimensional; use a single-value scale.")
    x = values[:, 0]
    return x - x.mean(), scale if isinstance(scale, str) else "custom"

################################# Spectral Descriptors ####################################

def spectral_descriptors(sequence: str, transform: str = "fft", scale: str | dict | None = None,
    n: int = 10) -> pd.DataFrame:
    """
    Calculate Spectral descriptors of protein sequence. Each amino acid in the
    sequence is replaced by its value on a property scale and the mean-centred
    signal x is transformed into a fixed number of features.

    transform="fft": the power spectrum of the signal [1],

    P(f) = abs(FFT(x)(f))^2 / L

    where f is the frequency in cycles per residue (0 to 0.5) and L is the length
    of the sequence. The frequency range 0-0.5 is split into n equal-width bins and
    the mean power in each bin is returned; bins containing no frequencies are 0.
    The default scale is the electron-ion interaction potential (EIIP) of the
    resonant recognition model [1].

    transform="dct": coefficients 1..n of the orthonormal type-II discrete cosine
    transform of the signal, which compresses it into a fixed number of values [2]:

    X_k = sqrt(2/L) * sum_{n=0}^{L-1} x_n * cos(pi * k * (2n + 1) / (2L))

    X_0 is always 0 for a mean-centred signal, so it is not returned. The default
    scale is the Kyte-Doolittle hydropathy.

    Parameters
    ==========
    :sequence: str
        protein sequence.
    :transform: str (default="fft")
        "fft" or "dct".
    :scale: str | dict | None (default=None)
        single-value scale name, AAindex accession number or dict of amino acid
        values; None uses "eiip" for FFT and "kyte_doolittle" for DCT.
    :n: int (default=10)
        number of FFT bins, or of DCT coefficients (less than the sequence length).

    Returns
    =======
    :spectral_df: pd.DataFrame
        pandas dataframe of shape 1 x n. Column names are of the form
        FFT_<scale>_bin<i> (i = 1..n) or DCT_<scale>_<k> (k = 1..n).
    """
    #validate input protein sequence
    sequence = _validate_sequence(sequence)
    if transform not in ("fft", "dct"):
        raise ValueError(f"Invalid transform '{transform}'. Valid transforms are: ['fft', 'dct'].")
    n = _validate_int("n", n, 1)
    if scale is None:
        scale = "eiip" if transform == "fft" else "kyte_doolittle"
    x, name = _signal(sequence, scale)
    L = len(x)

    if transform == "fft":
        power = np.abs(np.fft.rfft(x)) ** 2 / L
        freqs = np.fft.rfftfreq(L)
        #bin index of each frequency; f = 0.5 belongs to the last bin
        bins = np.minimum((freqs / 0.5 * n).astype(int), n - 1)
        sums = np.bincount(bins, weights=power, minlength=n)
        counts = np.bincount(bins, minlength=n)
        values = np.divide(sums, counts, out=np.zeros(n), where=counts > 0)
        return _to_df(list(values), [f"FFT_{name}_bin{i}" for i in range(1, n + 1)])

    if n >= L:
        raise ValueError(f"n ({n}) must be less than the sequence length ({L}) for the DCT.")
    k = np.arange(1, n + 1)[:, None]
    basis = np.cos(np.pi * k * (2 * np.arange(L) + 1) / (2 * L))
    #elementwise sum instead of matmul: macOS Accelerate BLAS raises spurious float warnings
    coeffs = np.sqrt(2 / L) * (basis * x).sum(axis=1)
    return _to_df(list(coeffs), [f"DCT_{name}_{i}" for i in range(1, n + 1)])

def fft_descriptors(sequence: str, scale: str | dict = "eiip", n_bins: int = 10) -> pd.DataFrame:
    """
    Calculate Fast Fourier Transform (FFT) descriptors of protein sequence: the
    binned power spectrum of a property signal [1].

    Deprecated since v1.6.0, removed in 2.0: use
    ``spectral_descriptors(sequence, transform="fft", scale=scale, n=n_bins)``.

    Parameters
    ==========
    :sequence: str
        protein sequence.
    :scale: str | dict (default="eiip")
        single-value scale name, AAindex accession number or dict of amino acid values.
    :n_bins: int (default=10)
        number of frequency bins.

    Returns
    =======
    :fft_df: pd.DataFrame
        pandas dataframe of shape 1 x n_bins, with columns FFT_<scale>_bin<n>.
    """
    _warn_deprecated("fft_descriptors", 'spectral_descriptors(sequence, transform="fft")')
    return spectral_descriptors(sequence, "fft", scale, n_bins)

def dct_descriptors(sequence: str, scale: str | dict = "kyte_doolittle", n_coeffs: int = 10) -> pd.DataFrame:
    """
    Calculate Discrete Cosine Transform (DCT) descriptors of protein sequence:
    coefficients of the type-II DCT of a property signal [2].

    Deprecated since v1.6.0, removed in 2.0: use
    ``spectral_descriptors(sequence, transform="dct", scale=scale, n=n_coeffs)``.
    Since v1.6.0 the always-zero coefficient X_0 is dropped, so the coefficients
    are 1..n_coeffs.

    Parameters
    ==========
    :sequence: str
        protein sequence.
    :scale: str | dict (default="kyte_doolittle")
        single-value scale name, AAindex accession number or dict of amino acid values.
    :n_coeffs: int (default=10)
        number of coefficients, less than the length of the sequence.

    Returns
    =======
    :dct_df: pd.DataFrame
        pandas dataframe of shape 1 x n_coeffs, with columns DCT_<scale>_<k>.
    """
    _warn_deprecated("dct_descriptors", 'spectral_descriptors(sequence, transform="dct")')
    return spectral_descriptors(sequence, "dct", scale, n_coeffs)
