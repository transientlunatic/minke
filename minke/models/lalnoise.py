"""
Noise models from lalsimulation.
"""

try:
    import torch
    TORCH_AVAILABLE = True
except ImportError:
    TORCH_AVAILABLE = False
import scipy
import numpy as np

import lal
import lalsimulation

from ..types import PSD, TimeSeries
from . import PSDApproximant


class LALSimulationPSD(PSDApproximant):
    """A class to generate power spectral densities using LALSimulation."""

    def __init__(self):

        super().__init__()

    def frequency_domain(
        self,
        df=1,
        frequencies=None,
        lower_frequency=20,
        upper_frequency=1024,
        mask_below=20,
    ):
        if frequencies is None:
            if TORCH_AVAILABLE:
                frequencies = torch.arange(lower_frequency, upper_frequency + df, df)
            else:
                frequencies = np.arange(lower_frequency, upper_frequency + df, df)

        N = int(len(frequencies))
        df = float(frequencies[1] - frequencies[0])
        psd_data = lal.CreateREAL8FrequencySeries(
            None, lal.LIGOTimeGPS(0), lower_frequency, df, lal.HertzUnit, N
        )
        self.psd_function(psd_data, flow=lower_frequency)
        psd_data = psd_data.data.data
        #psd_data[frequencies < mask_below] = psd_data[frequencies > mask_below][0]
        psd = PSD(psd_data, frequencies=frequencies)
        return psd

    def twocolumn(self, *args, **kwargs):
        """
        Produce the PSD in two-column format.
        """
        psd = self.frequency_domain(*args, **kwargs)
        frequencies = psd.frequencies.value
        data = np.array(psd.data)

        return np.vstack([frequencies, data]).T

    
    def covariance_matrix(self, times):
        """
        Return a time-domain representation of this power spectral density.
        """
                
        dt = times[1] - times[0]
        N = len(times)
        T = times[-1] - times[0]
        df = 1 / T
        if TORCH_AVAILABLE:
            frequencies = torch.arange(len(times) // 2 + 1) * df.value
        else:
            frequencies = np.arange(len(times) // 2 + 1) * df
        psd = np.array(self.frequency_domain(df=df, frequencies=frequencies).data)
        psd[-1] = psd[-2]
        # import matplotlib.pyplot as plt
        # f, ax = plt.subplots(1,1)
        # ax.plot(frequencies, psd)
        # f.savefig("psd.png")
        # Calculate the autocovariance from a one-sided PSD
        acf = 0.5*np.real(np.fft.irfft(psd*df, n=(N)))*T
        # The covariance is then the Toeplitz matrix formed from the acf
        # f, ax = plt.subplots(1,1)
        # ax.plot(acf)
        # f.savefig("acf.png")
        return scipy.linalg.toeplitz(acf)

    def time_domain(self, times):
        return self.covariance_matrix(times)

    def time_series(self, times=None, duration=None, sample_rate=None, **kwargs):
        """Create a timeseries filled with noise from a specific PSD.
        Parameters
        ----------
        psd : `heron:PSD`
           The power spectral density of the noise to be generated in the time
           series.
        times : array-like
           The time axis of the timeseries.
        device : str or torch.device
           The device which the noise series should be stored in.

        Notes
        -----
        In order to follow with the conventions in lalsuite this has been
        adapted from the code in lalnoise.
        """

        if times is None and sample_rate is not None and duration is not None:
            dt = 1./sample_rate
            N = int(duration * sample_rate)
            df = 1/duration
            times = np.linspace(0, duration-dt, N)
            T = duration
        else:
            dt = times[1] - times[0]
            N = len(times)
            T = times[-1] - times[0]
            df = 1 / T
            
        if TORCH_AVAILABLE:
            frequencies = torch.arange(0, N // 2 + 1) * df
        else:
            frequencies = np.arange(0, N // 2 + 1) * df
        reals = np.random.randn(len(frequencies))
        imags = np.random.randn(len(frequencies))
        psd = np.array(self.frequency_domain(df=df, frequencies=frequencies).data)
        psd[-1] = psd[-2]

        S = 0.5 * np.sqrt(psd)# / df) #* T inside sqrt # np.sqrt(N * N / 4 / (T) * psd.value)

        noise_r = S * (reals)
        noise_i = S * (imags)

        noise_f = noise_r + 1j * noise_i

        times += (kwargs.get("epoch", 0))

        return TimeSeries(data=2*np.fft.irfft(noise_f, n=(N))*df*N, times=times, dt=dt)


class AdvancedLIGODesignSensitivity2018(LALSimulationPSD):
    psd_function = lalsimulation.SimNoisePSDaLIGODesignSensitivityT1800044


class AdvancedLIGO(AdvancedLIGODesignSensitivity2018):
    pass


KNOWN_PSDS = {"AdvancedLIGO": AdvancedLIGO}
