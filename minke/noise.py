import numpy as np
import lal, lalsimulation


class PSD():
    """
    Create a power spectral distribution from a file which specifies
    the PSD.

    """

    def __init__(self, filename, fmin = 20, fmax = 3000, rate = 16384.0):
        """
        Load in the PSD file, assuming that it is in units of Hertz and strain.
        """

        df = 1.0 / rate
        N = int((fmax - fmin) / df)
        self.df = df
        
        self.psd = lal.CreateREAL8FrequencySeries('psd', 
                                                  lal.LIGOTimeGPS(0,0),
                                                  fmin,
                                                  df,
                                                  lal.HertzUnit,
                                                  N)
        lalsimulation.SimNoisePSDFromFile(self.psd, fmin, filename)

        
    def SNR(self, waveform, detectors, fmin=0, fmax=-1):
        """
        Calculate the SNR for a given waveform compared to this SNR.

        Parameters
        ==========
        waveform : minke source object
           The waveform which should have its snr measured
        detectors : list
           A list of detector names in the format "L1" for the Livingston 4km detector, etc.
        fmin : float
           The minimum frequency at which to evaluate the SNR. If 0 the lower frequency is unbounded.
        fmax : float
           The maximum frequency at which the SNR is evaluated. If negative the upper frequency is unbounded.


        Returns
        =======
        network snr : float
           The SNR over the whole network described by the detector list
        SNRs : list
           A list of SNRs for each detector
        """
        snrs = []
        for detector in detectors:
            snrs.append(lalsimulation.MeasureSNR(waveform._generate_for_detector([detector]), self.psd, fmin, fmax))

        snrs_r = np.array(snrs)
        snrs_r = snrs_r**2
        return np.sqrt(snrs_r.sum()), snrs
        

    def plot(self):
        """
        Plot the PSD
        """
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fmin = self.psd.f0
        length = len(self.psd.data.data)
        df = self.psd.deltaF
        frequencies = np.linspace(fmin, fmin + length * df, length)

        f = plt.semilogy(frequencies, self.psd.data.data)

        return f

class NoiseTimeseries(object):
    """
    Create timeseries of noise coloured by a PSD.
    """
    def __init__(self, psd, length, epoch, rate = 16384.0):
        """
        Create a noisy timeseries.
        
        Parameters
        ----------
        psd : `minke.noise.PSD`
           The power spectral density of the noise.
        """
        self.psd = psd
        randomness = lal.gsl_rng("ranlux", seed)
        seg = lal.CreateREAL8TimeSeries("STRAIN", epoch, 0.0, 1.0/rate, lal.StrainUnit, length)

        #lal.simNoise
