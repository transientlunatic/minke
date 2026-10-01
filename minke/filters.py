"""
This file defines some simple filters for use in minke
"""

import numpy as np

def inner_product(a, b, s=None):
    """
    Detetermine the inner product of two vectors with optional weighting

    Parameters
    ----------
    a, b : array-like
       The vectors over which to determine the inner product
    s : array-like
       The weighting vector. Default is None.
    """

    if s is not None:
        return np.sum(np.real(2 * ( a * np.conj(b)) / s)[2:-2])
    else:
        return np.sum(np.real(2 * (a * np.conj(b))))


def optimal_snr_squared(strain, psd_model, sample_rate, f_min=20.0, f_max=None):
    """
    Optimal SNR squared of a time-domain signal, rho^2 = 4 int |h(f)|^2 / S(f) df.

    This is the same definition used by bilby and simple-pe.

    Parameters
    ----------
    strain : array-like
       The noise-free time-domain signal, uniformly sampled.
    psd_model : minke PSD model
       Must provide ``frequency_domain(frequencies=..., lower_frequency=...)``.
    sample_rate : float
       Sample rate of ``strain`` in Hz.
    f_min, f_max : float
       Band over which to integrate. ``f_max`` defaults to Nyquist.
    """
    x = np.asarray(strain, dtype=float)
    sample_rate = float(sample_rate)
    N = len(x)
    df = sample_rate / N
    frequencies = np.fft.rfftfreq(N, 1.0 / sample_rate)
    h = np.fft.rfft(x) / sample_rate

    upper = sample_rate / 2 if f_max is None else f_max
    band = (frequencies >= f_min) & (frequencies <= upper)
    frequencies, h = frequencies[band], h[band]

    # The LAL-backed PSD models tabulate from lower_frequency, so the axis
    # passed in must start at that frequency.
    psd = np.asarray(
        psd_model.frequency_domain(
            frequencies=frequencies, lower_frequency=frequencies[0]
        ).data,
        dtype=float,
    )
    ok = np.isfinite(psd) & (psd > 0)
    return 4.0 * df * np.sum(np.abs(h[ok]) ** 2 / psd[ok])
