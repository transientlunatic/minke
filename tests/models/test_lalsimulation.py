"""
These tests are designed to test the GPyTorch-based models in heron.
"""

import unittest

import astropy.units as u
import torch
import numpy as np
import gpytorch

from heron.training.makedata import make_manifold, make_optimal_manifold
from heron.models.gpytorch import HeronNonSpinningApproximant
from heron.models.lalsimulation import SEOBNRv3, IMRPhenomPv2
from heron.models.lalnoise import AdvancedLIGO
from heron.filters import Overlap
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use("agg")

from . import _GenericWaveform

class test_IMRPhenomPv2(_GenericWaveform):

    @classmethod
    def setUpClass(cls):
        # initialize likelihood and model
        cls.model = IMRPhenomPv2()


class test_SEOBNRv3(_GenericWaveform):

    @classmethod
    def setUpClass(cls):
        # initialize likelihood and model
        cls.model = SEOBNRv3()
