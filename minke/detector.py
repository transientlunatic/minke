"""
Logic for handling detectors in heron, largely wrapping functions from lalsimulation.
"""

from collections import namedtuple

from lal import (
    cached_detector_by_prefix,
    TimeDelayFromEarthCenter,
    LIGOTimeGPS,
    antenna,
)
from lalinference import DetFrameToEquatorial

AntennaResponse = namedtuple("AntennaResponse", "plus cross")


class DetectorBase:
    pass


class Detector(DetectorBase):

    def antenna_response(self, ra: float, dec: float, psi: float, time):
        """
        Get the antenna responses for a given detector.

        Parameters
        ----------
        ra : float
           The right-ascension of the source, in radians.
        dec : float
           The declination of the source, in radians.
        psi : float
           The polarisation angle, in radians.
        time : float, or array of floats
           The GPS time, or an array of GPS times, at which
           the response should be evaluated.

        Returns
        -------
        plus : float, or array
           The 'plus' component of the response function.
        cross : float, or array
           The 'cross' component of the response function.
        """
        responses = antenna.AntennaResponse(
            self.abbreviation, ra, dec, psi=psi, times=time
        )
        return AntennaResponse(responses.plus, responses.cross)

    def geocentre_delay(self, ra: float, dec: float, times):
        """
        Calculate the delay in arrival time for a signal at
        this detector relative to the centre of the Earth.

        Parameters
        ----------
        ra : float
           The right-ascension of the source, in radians.
        dec : float
           The declination of the source, in radians.
        times : float, or array of floats
           The GPS time, or an array of GPS times, at which
           the response should be evaluated.

        Returns
        -------
        float
           The time in seconds by which the response is delayed.
        """
        dt = TimeDelayFromEarthCenter(
            self._lal_detector.location, ra, dec, LIGOTimeGPS(times)
        )
        return dt


class AdvancedLIGO(Detector):
    pass


class AdvancedVirgo(Detector):
    _lal_detector = cached_detector_by_prefix["V1"]
    abbreviation = "V1"


class AdvancedLIGOLivingston(Detector):
    _lal_detector = cached_detector_by_prefix["L1"]
    abbreviation = "L1"


class AdvancedLIGOHanford(Detector):
    _lal_detector = cached_detector_by_prefix["H1"]
    abbreviation = "H1"


KNOWN_IFOS = {
    "AdvancedLIGOHanford": AdvancedLIGOHanford,
    "AdvancedLIGOLivingston": AdvancedLIGOLivingston,
    "AdvancedVirgo": AdvancedVirgo,
}
