"""
This module provides functions to calculate antenna factors for a given time, a given sky location and a given detector.

Adapted from the implementation in pylal
"""
import sys
from math import *
import lal
import lalsimulation

__author__ = "Alexander Dietz <Alexander.Dietz@astro.cf.ac.uk>; Daniel Williams <daniel.williams@ligo.org>"

def response( gpsTime, rightAscension, declination, inclination,
              polarization, unit, det ):

  """
  response( gpsTime, rightAscension, declination, inclination,
              polarization, unit, detector )
  
  Calculates the antenna factors for a detector 'detector' (e.g. 'H1')
  at a given gps time (as integer) for a given sky location
  (rightAscension, declination) in some unit (degree/radians).
  This computation also takes into account a specific inclination
  and polarization.
  
  The returned values are: (f-plus, f-cross, f-average, q-value).
  
  Example: antenna.response( 854378604.780, 11.089, 42.308, 0, 0, 'radians', 'H1' )
  """
  
  # check the input arguments
  if unit =='radians':
    ra_rad = rightAscension
    de_rad = declination
    psi_rad = polarization
    iota_rad = inclination
  elif unit =='degree':
    ra_rad = rightAscension/180.0*pi
    de_rad = declination/180.0*pi
    psi_rad = polarization/180.0*pi
    iota_rad = inclination/180.0*pi
  else:
    raise ValueError("Unknown unit %s" % unit)

  # calculate GMST if the GPS time
  gps = lal.LIGOTimeGPS( gpsTime )
  gmst_rad = lal.GreenwichMeanSiderealTime(gps)

  # Get the detector from its prefix
  try:
    detector = lalsimulation.DetectorPrefixToLALDetector(det)
  except KeyError:
    raise ValueError("ERROR. Key %s is not a valid detector prefix." % (det))

  # get the correct response data
  response = detector.response

  # actual computation of antenna factors
  f_plus, f_cross = lal.ComputeDetAMResponse(response, ra_rad, de_rad,
                                                    psi_rad, gmst_rad)

  f_ave=sqrt( (f_plus*f_plus + f_cross*f_cross)/2.0 );
  ci=cos( iota_rad );
  cc=ci*ci;

  # calculate q-value, e.g. ratio of effective to real distance
  # ref: Duncans PhD, eq. (4.3) on page 57 
  f_q=sqrt( f_plus*f_plus*(1+cc)*(1+cc)/4.0 + f_cross*f_cross*cc ); 

  # output
  return f_plus, f_cross, f_ave, f_q



def timeDelay( gpsTime, rightAscension, declination, unit, det1, det2 ):
  """
  timeDelay( gpsTime, rightAscension, declination, unit, det1, det2 )
  
  Calculates the time delay in seconds between the detectors
  'det1' and 'det2' (e.g. 'H1') for a sky location at (rightAscension
  and declination) which must be given in certain units
  ('radians' or 'degree'). The time is passes as GPS time.
  A positive time delay means the GW arrives first at 'det2', then at 'det1'.
    
  Example:
  antenna.timeDelay( 877320548.000, 355.084,31.757, 'degree','H1','L1')
  0.0011604683260994519
  Given these values, the signal arrives first at detector L1,
  and 1.16 ms later at H2
  """

  # check the input arguments
  if unit =='radians':
    ra_rad = rightAscension
    de_rad = declination
  elif unit =='degree':
    ra_rad = rightAscension/180.0*pi
    de_rad = declination/180.0*pi
  else:
    raise ValueError("Unknown unit %s" % unit)

  # check input values
  if ra_rad<0.0 or ra_rad> 2*pi:
    raise ValueError( "ERROR. right ascension=%f "\
          "not within reasonable range."\
          % (rightAscension))

  if de_rad<-pi or de_rad> pi:
    raise ValueError( "ERROR. declination=%f not within reasonable range."\
          % (declination))

  if det1 == det2:
    return 0.0
  
  gps = lal.LIGOTimeGPS( gpsTime )

  x1 = lalsimulation.DetectorPrefixToLALDetector(det1).location
  x2 = lalsimulation.DetectorPrefixToLALDetector(det2).location
  timedelay = lal.ArrivalTimeDiff(list(x1), list(x2), ra_rad, de_rad, gps)

  return timedelay
