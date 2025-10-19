"""
Noise models from lalsimulation.
"""

import scipy
import numpy as np
try:
    import torch
    array_lib = torch
    array_lib_s = "torch"
except ImportError:
    array_lib = np
    array_lib_s = "numpy"

import lal
import lalsimulation

from .types import PSD, TimeSeries
from .models import PSDApproximant


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
            frequencies = array_lib.arange(lower_frequency, upper_frequency + df, df)
            f0 = float(lower_frequency)
        else:
            f0 = float(frequencies[0])

        N = int(len(frequencies))
        df = float(frequencies[1] - frequencies[0])
        psd_data = lal.CreateREAL8FrequencySeries(
            None, lal.LIGOTimeGPS(0), f0, df, lal.HertzUnit, N
        )
        self.psd_function(psd_data, flow=f0)
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
        if array_lib_s == "torch":
            df = df.value
        frequencies = array_lib.arange(len(times) // 2 + 1) * df
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
        """Create a time series of zero-mean Gaussian noise with this one-sided PSD.

        Normalization follows standard one-sided PSD conventions:
        - Interior frequency bins (k=1..N/2-1): X_k = sqrt(S1(f_k) * fs / 2) * (a_k + i b_k)
        - DC bin (k=0): X_0 = 0 (no offset)
        - Nyquist bin (k=N/2, when N even): X_{N/2} = sqrt(S1(f_{N/2}) * fs) * a_{N/2} (real)
        - Before applying the inverse real FFT, the frequency-domain noise is scaled by sqrt(N/2) to compensate for numpy's normalization.
        """

        # Establish time grid and frequency resolution
        if times is None and sample_rate is not None and duration is not None:
            N = int(duration * sample_rate)
            times = np.linspace(0, duration, N, endpoint=False)
            dt = 1.0 / sample_rate
        elif times is not None:
            times = np.asarray(times)
            N = len(times)
            dt = float(times[1] - times[0])
            sample_rate = 1.0 / dt
        else:
            raise ValueError("Provide either times, or duration and sample_rate")

        df = sample_rate / N
        freqs = array_lib.arange(0, N // 2 + 1) * df

        # One-sided PSD sampled on our frequency grid
        psd = np.asarray(self.frequency_domain(df=df, frequencies=freqs).data)
        if len(psd) > 1:
            psd[-1] = psd[-2]  # avoid edge artifacts at Nyquist

        # Random Gaussian draws
        reals = np.random.randn(len(freqs))
        imags = np.random.randn(len(freqs))

        # Construct one-sided frequency-domain noise
        noise_f = np.zeros(len(freqs), dtype=np.complex128)
        if len(freqs) > 2:
            amp_mid = np.sqrt(psd[1:-1] * sample_rate / 2.0)
            noise_f[1:-1] = amp_mid * (reals[1:-1] + 1j * imags[1:-1])

        # DC = 0 to avoid mean offset
        noise_f[0] = 0.0 + 0.0j

        # Nyquist bin is purely real with full factor (no 1/2)
        if N % 2 == 0:
            noise_f[-1] = np.sqrt(psd[-1] * sample_rate) * reals[-1]

        # Compensate for numpy's 1/N inverse FFT normalization
        noise_f *= np.sqrt(N / 2.0)

        # Epoch shift if requested
        times = times + (kwargs.get("epoch", 0))

        # Inverse real FFT; numpy uses backward normalization so this is consistent
        data = np.fft.irfft(noise_f, n=N)
        return TimeSeries(data=data, times=times)


class AdvancedLIGODesignSensitivity2018(LALSimulationPSD):
    psd_function = lalsimulation.SimNoisePSDaLIGODesignSensitivityT1800044


class AdvancedLIGO(AdvancedLIGODesignSensitivity2018):
    pass


KNOWN_PSDS = {"AdvancedLIGO": AdvancedLIGO}
