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
