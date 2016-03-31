import sys
import math
from optparse import OptionParser, Option, OptionGroup

from scipy import random
import numpy

from glue.ligolw import lsctables
from glue.ligolw import utils
from glue.ligolw import ligolw
from glue.ligolw import ilwd
from glue.segments import segment
from glue.lal import LIGOTimeGPS as GPS
from glue.ligolw.utils import process
from pylal.antenna import response

import lal
import lalburst
import lalsimulation

from minke.distribution import *

import matplotlib.pyplot as plt

class Waveform(object):
    """
    Generic container for different source types. 
    Currently, it checks for the waveform type and initializes itself appropriately. 
    In the future, different sources should subclass this and override the generation routines.
    """

    sim = lsctables.New(lsctables.SimBurstTable)

    numrel_data = []
    waveform = "Generic"
    expnum = 1

    def _clear_params(self):
        self.params = {}
        for a in lsctables.SimBurstTable.validcolumns.keys():
            self.params[a] = None
        

    def parse_polarisation(self, polarisation):
        """
        Convert a string description of a polarisation to an ellipse eccentricity and an ellipse angle.

        Parameters
        ----------
        polarisation : str, {'linear', 'circular', 'elliptical', 'inclination'}
           The description of the polarisation, in words.

        Outputs
        -------
        e : float
           The ellipse's eccentricity.
        angle : float
           The ellipse angle.
        """
        if polarisation == "linear":
            pol_ellipse_e = 1.0
            pol_ellipse_angle = 0
        elif polarisation == "circular":
            pol_ellipse_e = 0.0
            pol_ellipse_angle = 0
        elif polarisation == "elliptical":
            pol_ellipse_e = uniform_interval((0,1),1)[0]
            pol_ellipse_angle = uniform_interval((0,2*numpy.pi),1)[0]
        elif polarisation == "inclination":
            cosincl = uniform_interval((-1, 1), 1)[0]**2
            pol_ellipse_e = (1 - cosincl) / (1 + cosincl)
            pol_ellipse_angle = -numpy.pi/2 if uniform_interval((0, 1), 1)[0] < 0.5 else numpy.pi/2

        return pol_ellipse_e, pol_ellipse_angle

    def plot(self):
        """
        Produce a plot of the injection.
        """
        self._generate()
        plt.plot(self.hp.data.data, label="+ polarisation")
        plt.plot(self.hx.data.data, label="x polarisation")

    def _generate(self, rate=16384.0):
        """
        Generate the burst described in a given row, so that it can be 
        measured.
        
        Parameters
        ----------
        rate : float
            The sampling rate of the signal, in Hz. Defaults to 16384.0Hz
            
        Returns
        -------
        hp : 
            The strain in the + polarisation
        hx : 
            The strain in the x polarisation
        hp0 : 
            A copy of the strain in the + polarisation
        hx0 : 
            A copy of the strain in the x polarisation
        """
        row = self._row()
        self.swig_row = lalburst.CreateSimBurst()
        for a in lsctables.SimBurstTable.validcolumns.keys():
            try:
                setattr(self.swig_row, a, getattr( row, a ))
                print getattr(row,a)
            except AttributeError: continue
            except TypeError: 
                continue
        try:
            self.swig_row.numrel_data = row.numrel_data
        except:
            pass
        hp, hx = lalburst.GenerateSimBurst(self.swig_row, 1.0/rate)
        # FIXME: Totally inefficent --- but can we deep copy a SWIG SimBurst?
        # DW: I tried that, and it doesn't seem to work :/
        hp0, hx0 = lalburst.GenerateSimBurst(self.swig_row, 1.0/rate)
        self.hp, self.hx, self.hp0, self.hx0 = hp, hx, hp0, hx0

    def _row(self, sim=None):
        """
        Produce a simburst table row for this waveform.

        Todo
        ----
        We also need to add the process_id back in.
        """
        if not sim: sim = self.sim
        self.row = sim.RowType()

        for a in lsctables.SimBurstTable.validcolumns.keys():
            setattr(self.row, a, self.params[a])
        
        self.row.waveform = self.waveform
        # Fill in the time
        self.row.set_time_geocent(GPS(float(self.time)))
        # Get the sky locations
        self.row.ra, self.row.dec, self.row.psi = self.sky_dist()
        self.row.simulation_id = sim.get_next_id()
        self.row.waveform_number = random.randint(0,int(2**32)-1)
        ### !! This needs to be updated.
        self.row.process_id = "process:process_id:0" #procrow.process_id
        self.row.time_slide_id = ilwd.ilwdchar("time_slide:time_slide_id:%d" % 1)#options.time_slide_id)

        return self.row



class SineGaussian(Waveform):
    """
    A class to represent a SineGaussian injection.
    """
    waveform = "SineGaussian"
    
    def __init__(self, q, frequency, hrss, polarisation, time, sky_dist=uniform_sky, seed=0):
        """A class to represent a SineGaussian ad-hoc waveform.

        Parameters
        ----------
        q : float
           The quality factor.

        frequency : float
           The frequency of the injection.

        hrss : float
           The strain magnitude of the injection.

        polarisation : str {'linear', 'elliptical', 'circular'}
           The type of polarisation of the waveform.

        time : float
           The central time of the injection.

        sky_dist : func
           The function describing the sky distribution which the injections
           should be made over. Defaults to a uniform sky.

        seed : int
           The random seed used to make the injection time of the waveform.
           The default seed is 0.

        """
        self._clear_params()
        self.sky_dist = sky_dist
        self.params['hrss'] = hrss
        self.params['frequency'] = frequency
        self.params['q'] = q
        self.time = time
        self.polarisation = polarisation
        self.params['pol_ellipse_e'], self.params['pol_ellipse_angle'] = self.parse_polarisation(self.polarisation)    


class Gaussian(Waveform):
    """
    A class to represent a Gaussian injection.
    """

    waveform = "Gaussian"

    def __init__(self, duration, hrss, time, sky_dist=uniform_sky, seed=0):
        """
        A class to represent a Gaussian ad-hoc waveform.

        Parameters
        ----------
        duration : float or list
           The duration, in milliseconds, of the Gaussian waveform.

        hrss : float or list
           The strain magnitude of the injection.
           If a float is provided then the hrss will be fixed, 
           if a list is provided then this will be the 
           minimum and maximum hrss.

        polarisation : str {'linear', 'elliptical', 'circular'}
           The type of polarisation of the waveform.

        time : float or list 
           The time period over which the injection should be made. If
           a list is given they should be the start and end times, and
           the waveform will be produced at some random point in that
           time range. If a float is given then the injection will be
           made at that specific time.

        sky_dist : func
           The function describing the sky distribution which the injections
           should be made over. Defaults to a uniform sky.

        seed : float 
           The random seed used to make the injection time of the waveform.
           The default seed is 0.

        """
        self._clear_params()
        self.sky_dist = sky_dist
        self.params['duration'] = duration
        self.params['hrss'] = hrss
        self.time = time
        self.params['pol_ellipse_e'] = 1.0
        self.params['pol_ellipse_angle'] = 0
