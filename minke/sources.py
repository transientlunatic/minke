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


class Waveform(object):
    """
    Generic container for different source types. 
    Currently, it checks for the waveform type and initializes itself appropriately. 
    In the future, different sources should subclass this and override the generation routines.
    """

    numrel_data = []
    waveform = "Generic"
    expnum = 1
        
    def log_hrss(self):
        """
        Draw uniformly in the log of a predefined hrss range
        """
        log10h = segment(numpy.log10(self.h_values[0]), numpy.log10(self.h_values[1]))
        return 10**uniform_interval(log10h, 1)

    def log_egw(self):
        """
        Draw uniformly in the log of a predefined E_gw range
        """
        log10e = segment(numpy.log10(self.egw_range[0]), numpy.log10(self.egw_range[1]))
        return 10**uniform_interval(log10e, self.expnum)

    def uniform_sky(self):
        """
        Get a set of (RA, declination, polarization) randomized appopriately to astrophysical sources isotropically distributed in the sky.
        """
        expnum = self.expnum
        ra = uniform_phi(expnum)
        dec = uniform_dec(expnum)
        pol = uniform_phi(expnum)
        return ra, dec, pol

    def favorable_sky(net, time):
        """
        Wander through the skies, searching for a most favorable location --- draw extrinsic parameters as if the network antenna pattern magnitude were the PDF.
        """
        ndraws = len(time)
        ra_out, dec_out, psi_out = numpy.empty((3, len(time)))
        while ndraws > 0:
            ra = numpy.random.uniform(0, 2 * numpy.pi, 1000)
            dec = numpy.random.uniform(-numpy.pi / 2, numpy.pi / 2, 1000)
            psi = numpy.random.uniform(0, 2 * numpy.pi, 1000)
            rnd = numpy.random.uniform(0, 1, 1000)
            for r, d, p, n in zip(ra, dec, psi, rnd):
                if n < sky_params(net, time, r, d, p):
                    ndraws -= 1
                    ra_out[ndraws] = r
                    dec_out[ndraws] = d
                    psi_out[ndraws] = p
                    break
        return ra_out, dec_out, psi_out

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

    def uniform_time(self):
        """
        Get a set of randomized (integer) event times.
        """
        return random.randint(self.tstart, self.tstop, self.expnum) + random.rand(self.expnum)

    def _row(self, sim):
        """
        Produce a simburst table row for this waveform.

        Todo
        ----
        This can currently only make injections on a uniform sky. This should be fixed to take a generic distribution function.
        """
        self.row = sim.RowType()
        # Required columns not defined makes ligolw unhappy
        for a in lsctables.SimBurstTable.validcolumns.keys():
            setattr(self.row, a, None)
        self.row.waveform = self.waveform
        self.row.set_time_geocent(self.time)
        # Right now this only does uniform sky distributions, but we should provide a way to do /any/ distribution.
        self.row.ra, self.row.dec, self.row.psi = self.uniform_sky()
        self.row.simulation_id = sim.get_next_id()
        self.row.waveform_number = random.randint(0,int(2**32)-1)
        self.row.process_id = procrow.process_id
        self.row.time_slide_id = ilwd.ilwdchar("time_slide:time_slide_id:%d" % options.time_slide_id)




class SineGaussian(Waveform):
    """
    A class to represent a SineGaussian injection.
    """
    waveform = "SineGaussian"
    
    def __init__(self, q, frequency, hrss, polarisation, time, seed=0):
        """A class to represent a SineGaussian ad-hoc waveform.

        Parameters
        ----------
        q : float or list
           The quality factor.
           If a float is provded then the q-factor is fixed. If a list is provided then they are the
           minimum and maximum quality factor of the injections.

        frequency : float or list
           The frequency of the injection.
           If a float is provided the frequency will be fixed, if a list is provided then this will be the
           minimum and maximum frequencies.

        hrss : float or list
           The strain magnitude of the injection.
           If a float is provided then the hrss will be fixed, if a list is provided then this will be the 
           minimum and maximum hrss.

        polarisation : str {'linear', 'elliptical', 'circular'}
           The type of polarisation of the waveform.

        time : float or list 
           The time period over which the injection should be made. If
           a list is given they should be the start and end times, and
           the waveform will be produced at some random point in that
           time range. If a float is given then the injection will be
           made at that specific time.

        seed : float 
           The random seed used to make the injection time of the waveform.
           The default seed is 0.

        """

        if isinstance(q, list):
            self.q_values = segment(q[0], q[1])
        else:
            self.q_values = q
            
        if isinstance(frequency, list):
            self.f_values = segment(frequency[0], frequency[1])
        else:
            self.f_values = frequency

        if isinstance(hrss, list):
            self.h_values = segment(hrss[0], hrss[1])
        else:
            self.h_values = hrss

        if isinstance(time, list):
            self.tstart, self.tstop = time[0], time[1]
            self.seed = seed
            random.seed(seed)
            self.time = random.randint(self.tstart, self.tstop, 1) + random.rand(1)
        else:
            self.time = time

    def draw_params(self):
        """
        Draw a set of parameters (hrss, q, frequency) appropriate for a set of sinegaussian type injections.

        Returns
        -------
        q : float
           Quality factor
        f : float
           Frequency
        h : float
           HRSS
        """
        h, q, f = [], [], []
        if type(self.q_values) == segment:
            q = uniform_interval(self.q_valies, 1)
        else:
            q  = self.q_values
            #q = map(lambda i: self.q_range[i], random.randint(0, len(self.q_range), self.expnum) )

        if type(self.h_values) == segment:
            h = self.log_hrss()
        else:
            h = self.h_values
            #h = map(lambda i: self.h_range[i], random.randint(0, len(self.h_range), self.expnum) )

        if type(self.f_values) == segment:
            f = uniform_interval(self.f_values, self.expnum)
        else:
            f = self.f_values
            #f = map(lambda i: self.f_range[i], random.randint(0, len(self.f_range), self.expnum))

        return q, f, h

    def _row(self, sim):
        """
        Produce a simburst table row for this waveform.

        Todo
        ----
        This can currently only make injections on a uniform sky. This should be fixed to take a generic distribution function.
        """
        self.row = sim.RowType()
        # Required columns not defined makes ligolw unhappy
        for a in lsctables.SimBurstTable.validcolumns.keys():
            setattr(self.row, a, None)
        self.row.waveform = self.waveform
        self.row.set_time_geocent(self.time)
        # Right now this only does uniform sky distributions, but we should provide a way to do /any/ distribution.
        self.row.ra, self.row.dec, self.row.psi = self.uniform_sky()
        self.row.simulation_id = sim.get_next_id()
        self.row.waveform_number = random.randint(0,int(2**32)-1)
        self.row.process_id = procrow.process_id
        self.row.time_slide_id = ilwd.ilwdchar("time_slide:time_slide_id:%d" % options.time_slide_id)

        self.row.q, self.row.frequency, self.row.hrss = self.draw_params()
        self.row.pol_ellipse_e, self.row.pol_ellipse_angle = self.parse_polarisation(self.polarisation)

        return self.row

    #def _repr_html_(self):

    
