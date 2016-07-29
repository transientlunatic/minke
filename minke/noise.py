import numpy as np
import lal, lalsimulation


class PSD():
    """
    Create a power spectral distribution from a file which specifies
    the PSD.

    """

    def __init__(self, filename, fmin = 20, fmax = 3000, df = 0.1):
        """
        Load in the PSD file, assuming that it is in units of Hertz and strain.
        """

        N = int((fmax - fmin) / df)
        self.psd = lal.CreateREAL8FrequencySeries('psd', 
                                                  lal.LIGOTimeGPS(0,0),
                                                  fmin,
                                                  df,
                                                  lal.HertzUnit,
                                                  N)
        lalsimulation.SimNoisePSDFromFile(self.psd, fmin, filename)

        
    def SNR(self, waveform, detectors):
        """
        Calculate the SNR for a given waveform compared to this SNR.

        Parameters
        ==========
        waveform : minke source object
           The waveform which should have its snr measured
        detectors : list
           A list of detector names in the format "L1" for the Livingston 4km detector, etc.

        Returns
        =======
        network snr : float
           The SNR over the whole network described by the detector list
        SNRs : list
           A list of SNRs for each detector
        """
        snrs = []
        for detector in detectors:
            snrs.append(lalsimulation.MeasureSNR(waveform._generate_for_detector([detector]), self.psd))

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
