import logging

import numpy as array_library
import numpy as np
import torch
from scipy.interpolate import CubicSpline
import scipy

# Astropy to handle units sanely
from astropy import units as u

# Import LAL tools
import lalsimulation
import lal
from lal import cached_detector_by_prefix, TimeDelayFromEarthCenter, LIGOTimeGPS

# Import heron types
from ..types import Waveform, WaveformDict, PSD
from . import WaveformApproximant, PSDApproximant


class LALBurstApproximant(WaveformApproximant):
    """
    This is the base class for LALSimulationBurst models
    """

    def __init__(self):
        self._cache_key = {}
        self._args = {
            "m1": None,
            "m2": None,
            "S1x": 0.0,
            "S1y": 0.0,
            "S1z": 0.0,
            "S2x": 0.0,
            "S2y": 0.0,
            "S2z": 0.0,
            "distance": 10 * u.Mpc,
            "inclination": 0,
            "phi ref": 0.0 * u.Hertz,
            "longAscNodes": 0.0,
            "eccentricity": 0.0,
            "meanPerAno": 0.0,
            "delta T": 1 / (4096.0 * u.Hertz),
            "f_min": 20.0 * u.Hertz,
            "f_ref": 20.0 * u.Hertz,
            "params": lal.CreateDict(),
            "approximant": None,
        }

class SineGaussian(LALBurstApproximant):

    def time_domain(self, parameters, times=None, sample_rate=None, **kwargs):
        """
        Create the time domain representation of a sine gaussian.
        """

        if sample_rate is not None:
            deltaT = (1./sample_rate)
        elif times is not None:
            deltaT = times[1] - times[0]
        else:
            deltaT = (1./16384) * u.second
        epoch = 0 #FIXME
        hp, hx = lalsimulation.SimBurstSineGaussian(parameters['q'],
                                                  parameters['centre_frequency'],
                                                  parameters['hrss'],
                                                  parameters['eccentricity'],
                                                  parameters['phase'],
                                                  deltaT.value,
                                                  )

        hp_ts = Waveform(data=hp.data.data, dt=hp.deltaT, t0=hp.epoch + epoch)
        hx_ts = Waveform(data=hx.data.data, dt=hx.deltaT, t0=hx.epoch + epoch)

        waveform = WaveformDict(parameters=parameters, plus=hp_ts, cross=hx_ts)
        
        return waveform
