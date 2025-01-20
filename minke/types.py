from itertools import cycle

# GWPy to help with timeseries
from gwpy.timeseries import TimeSeriesBase, TimeSeries
from gwpy.frequencyseries import FrequencySeries

from lal import cached_detector_by_prefix, TimeDelayFromEarthCenter, LIGOTimeGPS
from lalinference import DetFrameToEquatorial

import numpy as array_library
import matplotlib.pyplot as plt


class TimeSeries(TimeSeries):
    """
    Overload the GWPy timeseries so that additional methods can be defined upon it.
    """

    pass


class PSD(FrequencySeries):
    def __init__(self, data, frequencies, *args, **kwargs):
        super(PSD).__init__(*args, **kwargs)


class WaveformBase(TimeSeries):
    def __init__(self, *args, **kwargs):
        super(WaveformBase).__init__()


class Waveform(WaveformBase):
    def __init__(self, variance=None, covariance=None, *args, **kwargs):
        # if "covariance" in kwargs:
        #     self.covariance = kwargs.pop("covariance")
        self.covariance = covariance
        self.variance = variance
        super(Waveform, self).__init__(*args, **kwargs)

    def __new__(self, variance=None, covariance=None, *args, **kwargs):
        # if "covariance" in kwargs:
        #     self.covariance = kwargs.pop("covariance")
        waveform = super(Waveform, self).__new__(TimeSeriesBase, *args, **kwargs)
        waveform.covariance = covariance
        waveform.variance = variance

        return waveform

    # @property
    # def dt(self):
    #     return self.waveform.times[1] - self.waveform.times[0]

    def align(self, waveform_b):
        """
        Align this waveform with another one by altering the phase.
        """
        pass


class WaveformDict:
    def __init__(self, parameters=None, **kwargs):
        self.waveforms = kwargs
        self._parameters = parameters

    def __getitem__(self, item):
        return self.waveforms[item]

    @property
    def parameters(self):
        return self._parameters

    @property
    def hrss(self):
        if "plus" in self.waveforms and "cross" in self.waveforms:
            return array_library.sqrt(
                self.waveforms["plus"] ** 2 + self.waveforms["cross"] ** 2
            )
        else:
            raise NotImplementedError

    def project(
        self,
        detector,
        ra=None,
        dec=None,
        psi=None,
        time=None,
        iota=None,
        phi_0=None,
        **kwargs
    ):
        """
        Project this waveform onto a detector.

        Parameters
        ----------
        detector : `heron.detectors.Detector`
          The detector onto which the waveform should be projected.
        ra : float, optional
          The roght ascension of the signal source
        dec : float, optional
          The declination of the signal source.
        """

        
        if not time:
            time = self.waveforms["plus"].epoch.value
        
        if ((ra is None) and (dec is None)) and (
            ("ra" in self._parameters) and ("dec" in self._parameters)
        ):
            ra = self._parameters["ra"]
            dec = self._parameters["dec"]

            dt = detector.geocentre_delay(ra=ra, dec=dec, times=time)

        elif ("azimuth" in self._parameters.keys()) and ("zenith" in self._parameters.keys()) and ("reference_frame" in self._parameters.keys()):
            # Use horizontal coordinates.
            det1 = cached_detector_by_prefix[self._parameters["reference_frame"][0]]
            det2 = cached_detector_by_prefix[self._parameters["reference_frame"][1]]
            tg, ra, dec = DetFrameToEquatorial(
                det1, det2, time, self._parameters["azimuth"], self._parameters["zenith"]
            )
            dt = time - tg
        elif (ra is None) and (dec is None):
            raise ValueError("Right ascension and declination must both be specified.")

        else:
            dt = detector.geocentre_delay(ra=ra, dec=dec, times=time)
        if "plus" in self.waveforms and "cross" in self.waveforms:
            
            if not iota and "theta_jn" in self._parameters:
                iota = self._parameters["theta_jn"]
            elif isinstance(iota, type(None)):
                raise ValueError("Theta_jn must be specified!")

            if not phi_0 and "phase" in self._parameters:
                phi_0 = self._parameters["phase"]
            elif isinstance(phi_0, type(None)):
                raise ValueError("Initial phase must be specified!")

            if not psi and "psi" in self._parameters:
                psi = self._parameters["psi"]
            elif isinstance(psi, type(None)):
                raise ValueError("Polarisation must be specified!")

            print(psi)
            response = detector.antenna_response(ra, dec, psi, time=time)

            plus_prefactor = (
                array_library.cos(phi_0)
                * (1 + array_library.cos(iota) ** 2)
                * response.plus
                + array_library.sin(phi_0) * array_library.cos(iota) * response.cross
            )
            cross_prefactor = (
                array_library.cos(phi_0) * array_library.cos(iota) * response.cross
                - array_library.sin(phi_0)
                * (1 + array_library.cos(iota) ** 2)
                * response.plus
            )

            projected_data = (
                self.waveforms["plus"].data * plus_prefactor
                + self.waveforms["cross"].data * cross_prefactor
            )

            if self.waveforms["plus"].variance is not None:
                projected_variance = (
                    self.waveforms["plus"].variance * plus_prefactor**2
                    + self.waveforms["cross"].variance * cross_prefactor**2
                )
            else:
                projected_variance = None

            if self.waveforms["plus"].covariance is not None:
                projected_covariance = (
                    self.waveforms["plus"].covariance * plus_prefactor**2
                    + self.waveforms["cross"].covariance * cross_prefactor**2
                )
            else:
                projected_covariance = None

            bins = dt / (self.waveforms["plus"].dt)
                
            projected_waveform = Waveform(
                data=array_library.roll(array_library.pad(projected_data, 5000), int(bins.value))[5000:-5000],
                variance=projected_variance,
                covariance=projected_covariance,
                times=self.waveforms["plus"].times,
            )

            return projected_waveform

        else:
            raise NotImplementedError


class WaveformManifold:
    """
    Store a manifold of different waveform points.
    """

    def __init__(self):
        self.locations = []
        self.data = []

    def add_waveform(self, waveforms: WaveformDict):
        self.locations.append(waveforms.parameters)
        self.data.append(waveforms)

    def array(self, component="plus", parameter="m1"):
        all_data = []
        for wn in range(len(self.locations)):
            data = array_library.array(
                list(
                    zip(
                        cycle([self.locations[wn][parameter]]),
                        self.data[wn][component].times.value,
                        self.data[wn][component].value,
                    )
                )
            )
            all_data.append(data)
        return array_library.vstack(all_data)

    def plot(self, component="plus", parameter="m1"):
        f, ax = plt.subplots(1, 1)
        for wn in range(len(self.locations)):
            data = array_library.array(
                list(
                    zip(
                        cycle([self.locations[wn][parameter]]),
                        self.data[wn][component].times.value,
                    )
                )
            )
            plt.scatter(data[:, 1], data[:, 0], c=self.data[wn][component], marker=".")
        return f
