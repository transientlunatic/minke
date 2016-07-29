import numpy as np
import lal, lalsimulation
import matplotlib.pyplot as plt

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

        

    def plot(self):
        """
        Plot the PSD
        """

        fmin = self.psd.f0
        length = len(self.psd.data.data)
        df = self.psd.deltaF
        frequencies = np.linspace(fmin, fmin + length * df, length)

        f = plt.semilogy(frequencies, self.psd.data.data)

        return f
