################################################################################
#############                 Signal Processing                   ##############
################################################################################

from __future__ import annotations

"""
Signal-processing descriptors of amino-acid property signals: binned FFT power
spectrum and DCT coefficients. numpy only.

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
from ._helpers import _property_matrix, _to_df

__all__ = ["fft_descriptors", "dct_descriptors"]


def _signal(sequence: str, scale: str | dict) -> tuple[np.ndarray, str]:
    """Mean-centred single-value property signal and the name used in column labels."""
    values = _property_matrix(sequence, scale)
    if values.shape[1] != 1:
        raise ValueError(f"Scale '{scale}' is multi-dimensional; use a single-value scale.")
    x = values[:, 0]
    return x - x.mean(), scale if isinstance(scale, str) else "custom"

#################################### FFT ##################################################

def fft_descriptors(sequence: str, scale: str | dict = "eiip", n_bins: int = 10) -> pd.DataFrame:
    """
    Calculate a length-invariant binned power spectrum of a property signal, as
    in the resonant recognition model of Cosic [1]. Each residue is replaced by
    its scale value (default: electron-ion interaction potential, EIIP), the
    signal is mean-centred, and its power spectrum is

    P(f) = |FFT(x)(f)|^2 / L,   f = 0 .. 0.5 cycles/residue

    The normalised frequency range 0-0.5 is divided into ``n_bins`` equal-width
    bins (the last bin includes 0.5) and the mean power in each bin is returned;
    bins containing no frequency are 0.

    Parameters
    ==========
    :sequence: str
        protein sequence.
    :scale: str | dict (default="eiip")
        single-value scale name, AAindex1 accession code or {aa: value} dict.
    :n_bins: int (default=10)
        number of frequency bins.

    Returns
    =======
    :fft_df: pd.DataFrame
        pandas DataFrame of shape 1 x n_bins, columns ``FFT_{scale}_bin{1..n}``.
    """
    #validate input protein sequence
    sequence = _validate_sequence(sequence)
    n_bins = _validate_int("n_bins", n_bins, 1)

    x, name = _signal(sequence, scale)
    power = np.abs(np.fft.rfft(x)) ** 2 / len(x)
    freqs = np.fft.rfftfreq(len(x))
    #bin index of each frequency; f = 0.5 belongs to the last bin
    bins = np.minimum((freqs / 0.5 * n_bins).astype(int), n_bins - 1)
    sums = np.bincount(bins, weights=power, minlength=n_bins)
    counts = np.bincount(bins, minlength=n_bins)
    values = np.divide(sums, counts, out=np.zeros(n_bins), where=counts > 0)
    return _to_df(list(values), [f"FFT_{name}_bin{i}" for i in range(1, n_bins + 1)])

#################################### DCT ##################################################

def dct_descriptors(sequence: str, scale: str | dict = "kyte_doolittle", n_coeffs: int = 10) -> pd.DataFrame:
    """
    Calculate the first ``n_coeffs`` coefficients of the orthonormal type-II
    discrete cosine transform (DCT-II) of a mean-centred property signal [2]:

    X_k = s_k sum_{n=0}^{L-1} x_n cos(pi k (2n + 1) / (2L)),
    s_0 = sqrt(1/L), s_k = sqrt(2/L) for k > 0

    X_0 is 0 because the signal is mean-centred.

    Parameters
    ==========
    :sequence: str
        protein sequence.
    :scale: str | dict (default="kyte_doolittle")
        single-value scale name, AAindex1 accession code or {aa: value} dict.
    :n_coeffs: int (default=10)
        number of coefficients, at most the sequence length.

    Returns
    =======
    :dct_df: pd.DataFrame
        pandas DataFrame of shape 1 x n_coeffs, columns ``DCT_{scale}_{0..n-1}``.
    """
    #validate input protein sequence
    sequence = _validate_sequence(sequence)
    n_coeffs = _validate_int("n_coeffs", n_coeffs, 1)
    if n_coeffs > len(sequence):
        raise ValueError(f"n_coeffs ({n_coeffs}) cannot exceed the sequence length ({len(sequence)}).")

    x, name = _signal(sequence, scale)
    L = len(x)
    k = np.arange(n_coeffs)[:, None]
    basis = np.cos(np.pi * k * (2 * np.arange(L) + 1) / (2 * L))
    scale_factors = np.where(k[:, 0] == 0, np.sqrt(1 / L), np.sqrt(2 / L))
    #elementwise sum instead of matmul: macOS Accelerate BLAS raises spurious float warnings
    coeffs = scale_factors * (basis * x).sum(axis=1)
    return _to_df(list(coeffs), [f"DCT_{name}_{i}" for i in range(n_coeffs)])
