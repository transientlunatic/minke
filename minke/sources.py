import sys
import math
from optparse import OptionParser, Option, OptionGroup

from scipy import random
import numpy
np = numpy
from glue.ligolw import lsctables
from glue.ligolw import utils
from glue.ligolw import ligolw
from glue.ligolw import ilwd
from glue.segments import segment
from glue.lal import LIGOTimeGPS as GPS
from glue.ligolw.utils import process
from pylal.antenna import response

import scipy.signal as signal
import scipy.interpolate as interp

import os.path


import lal
import lalburst
import lalsimulation
lalsim = lalsimulation

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
        f, ax = plt.subplots(1,2)
        times = np.arange(0, self.hp.deltaT*len(self.hp.data.data), self.hp.deltaT)
        ax[0].plot(times, self.hp.data.data, label="+ polarisation")
        ax[0].plot(times, self.hx.data.data, label="x polarisation")
        ax[1].plot(self.hp.data.data, self.hx.data.data)

    def _generate(self, rate=16384.0, half=False):
        """
        Generate the burst described in a given row, so that it can be 
        measured.
        
        Parameters
        ----------
        rate : float
            The sampling rate of the signal, in Hz. Defaults to 16384.0Hz
            
        half : bool
           Only compute the hp and hx once if this is true; 
           these are only required if you need to compute the 
           cross products. Defaults to False.
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
        if not half :
            hp0, hx0 = lalburst.GenerateSimBurst(self.swig_row, 1.0/rate)
        else:
            hp0, hx0 = hp, hx
        self.hp, self.hx, self.hp0, self.hx0 = hp, hx, hp0, hx0
        return hp, hx, hp0, hx0

    def _row(self, sim=None, slide_id=1):
        """
        Produce a simburst table row for this waveform.

        Parameters
        ----------
        sim : table
           The table which the row should be made for.
           If this is left empty the table is assumed to be a 
           sim_burst_table.

        slide_id : int
           The timeslide id. Defaults to 1.
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
        self.row.time_slide_id = ilwd.ilwdchar("time_slide:time_slide_id:%d" % slide_id)

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
        self.params['seed'] = seed
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
           The duration, in seconds, of the Gaussian waveform.

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


class WhiteNoiseBurst(Waveform):
    """
    A class to represent a WNB injection.
    """

    waveform = "BTLWNB"

    def __init__(self, duration, bandwidth, frequency, time, hrss=None, egw=None, sky_dist=uniform_sky, seed=0):
        """A class to represent a white-noise burst ad-hoc waveform.

        Parameters
        ----------
        duration : float or list
           The duration, in seconds, of the WNB.

        bandwidth : float or list
           The bandwidth, in hertz, of the WNB.

        frequency : float or list
           The frequency, in hertz, of the WNB.

        hrss : float or list 
           The strain magnitude of the injection.
           If a float is provided then the hrss will be fixed, if a
           list is provided then this will be the minimum and maximum
           hrss. If the hrss is not provided then you should provide
           an EGW value instead.

        egw : float
           The gravitational wave energy. 
           You should provide this if you do not provide the Hrss.

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


        To Do
        -----
        Add ability to create a WNB by giving the EGW rather than the strain.

        Notes
        -----
        See 
        http://software.ligo.org/docs/lalsuite/lalsimulation/group___l_a_l_sim_burst__h.html#ga0419dc37e5b83f18cd3bb34722ddac54
        for what this calls "under the hood" in LALSuite. There are some important considerations here
        with respect to the differing sample rates used at LIGO and VIRGO, and so when creating the WNB it's important that the 
        burst is created at a single sampel rate, and then resampled appropriately, so that the same waveform is used.

        """
        self._clear_params()
        self.sky_dist = sky_dist
        if hrss:
            self.params['hrss'] = hrss
        elif egw:
            self.params['egw'] = egw
        else:
            raise valueError('You need to provide either an hrss or an egw to produce a WNB waveform')
        # The burst group describes WNBs by their lowest frequency, but LALInference wants them at the central frequency,
        # so add half the bandwidth to get the central freq
        self.params['frequency'] = frequency + bandwidth / 2.0
        # We need a minimum window size so that the whole burst can be contained within it,
        # so expand the duration if it's too small.
        min_len =  np.sqrt( 4 * (np.pi**(-2) / bandwidth**2) )
        if duration < min_len: 
            self.params['duration'] = min_len + 1e-5
        else:
            self.params['duration'] = duration
        self.params['bandwidth'] = bandwidth
        self.time = time
        self.params['pol_ellipse_e'], self.params['pol_ellipse_angle'] = 0.0, 0.0
        self.params['egw_over_rsquared'] = hrss**2 * np.pi**2 * frequency**2 * lal.C_SI / lal.G_SI / lal.MSUN_SI * lal.PC_SI**2
        # I'm really not sure if we need to do this, but apparently the 
        # hrss of the actual waveform is not exactly what we ask for
        # the old pyBurst code measured this by generating the waveform
        # which seems wasteful, but I'll replicate it here anyway, for
        # consistency with the method used for O1.
        #hp, hx, _, _ = self._generate(half=True)
        #self.params['hrss'] =  lalsimulation.MeasureHrss(hp, hx)


class Supernova(Waveform):
    """

    A superclass to handle the spherial harmonic decompositions which
    all supernova waveforms require.

    """

    waveform = "Supernova" # We shouldn't ever use this anyway

    def construct_Hlm(self, Ixx, Ixy, Ixz, Iyy, Iyz, Izz, l=2, m=2):
        """
        Construct the expansion parameters Hlm from T1000553.  Returns the expansion
        parameters for l=2, m=m 
        """

        if l!=2:
            print "l!=2 not supported"
            sys.exit()
            if abs(m)>2:
                print "Only l=2 supported, |m| must be <=2"
                sys.exit()

        if m==-2:
            Hlm = np.sqrt(4*lal.PI/5) * (Ixx - Iyy + 2*1j*Ixy)
        elif m==-1:
            Hlm = np.sqrt(16*lal.PI/5) * (Ixx + 1j*Iyz)
        elif m==0:
            Hlm = np.sqrt(32*lal.PI/15) * (Izz - 0.5*(Ixx + Iyy))
        elif m==1:
            Hlm = np.sqrt(16*lal.PI/5) * (-1*Ixx + 1j*Iyz)
        elif m==2:
            Hlm = np.sqrt(4*lal.PI/5) * (Ixx - Iyy - 2*1j*Ixy)

        return Hlm

    def interpolate(self, x_old, y_old, x_new):
        """
        Convenience funtion to avoid repeated code
        """
        interpolator = interp.interp1d(x_old, y_old)
        return interpolator(x_new)

    def decompose(self, numrel_file, sample_rate = 16384.0, step_back = 0.01, distance = 10e-3):
        """
        Produce the spherial harmonic decompositions of a numerical
        waveform.
        
        Parameters
        ----------
        numrel_file : str
           The location of the numerical relativity waveform file.
        
        sample_rate : float
           The sample rate of the NR file. Defaults to 16384.0 Hz.
        
        step_back : float
           The amount of time, in seconds, of the data which should be included
           before the peak amplitude. Defaults to 0.01 sec.

        distance : float
           The distance, in megaparsecs, from the observer at which the NR waveforms were
           simulated. Defaults to 10 kpc (i.e. 10e-3 Mpc).

        Returns
        -------
        decomposition : ndarray
           The l=2 mode spherical decompositions of the waveform. 
        """

        # Load the times from the file
        data = np.loadtxt(numrel_file)
        data = data.T
        times = data[0]
        times -= times[0]

        # Load the I components from the file        
        Ixx, Ixy, Ixz, Iyy, Iyz, Izz = data[5:]

        # Make the new time vector for the requried sample rate
        target_times = np.arange(0, times[-1], 1.0/sample_rate)

        # Prepare the output matrix
        output = np.zeros((len(target_times), 11))

        # Add the times in to the first column of said matrix
        output[:, 0] = target_times

        
        for i, m in enumerate([-2,-1,0,1,2]):
            Hlm = self.construct_Hlm(Ixx, Ixy, Ixz, Iyy, Iyz, Izz, l=2, m=m)
            #
            # Resample to uniform spacing at 16384 kHz
            #
            Hlm_real = self.interpolate(times, Hlm.real, target_times)
            Hlm_imag = self.interpolate(times, Hlm.imag, target_times)
            #
            # Make the output, and rescale it into dimensionless strain values
            #
            output[:,2*(i+1)-1] = Hlm_real * np.sqrt(lal.G_SI / lal.C_SI**4) #/lal.MRSUN_SI / ( distance * lal.PC_SI * 1e6)
            output[:,2*(i+1)] =   -Hlm_imag * np.sqrt(lal.G_SI / lal.C_SI**4)#/lal.MRSUN_SI / ( distance * lal.PC_SI * 1e6)

        return output


class Scheidegger2010(Supernova):
    """
    The Scheidegger2010 waveform.
    """

    waveform = "Scheidegger+10"

    def __init__(self, phi, incl, time, sky_dist=uniform_sky, filepath="R0E1CA.txt", decomposed_path=None):
        """

        Parameters
        ----------
        phi : float
           The internal phi parameter of the supernova injection.
        
        incl : float
           The internal inclination parameter of the supernova injection.

        time : float or list 
           The time period over which the injection should be made. If
           a list is given they should be the start and end times, and
           the waveform will be produced at some random point in that
           time range. If a float is given then the injection will be
           made at that specific time.

        sky_dist : func
           The function describing the sky distribution which the injections
           should be made over. Defaults to a uniform sky.

        filepath : str
           The filepath to the numerical relativity waveform.

        decomposed_path : str
           The location where the decomposed waveform file should be stored. Optional.
        """
        
        self._clear_params()
        self.time = time
        self.params['phi'] = phi
        self.params['incl'] = incl
        self.sky_dist = sky_dist
        if not decomposed_path : decomposed_path = filepath+".dec"
        if not os.path.isfile(decomposed_path) :
            decomposed = self.decompose(filepath, sample_rate = 16384.0, step_back = 0.01, distance = 10e-3)
            
            np.savetxt(decomposed_path, decomposed, header="# time (2,-2) (2,-1) (2,0) (2,1) (2,2)", fmt='%.8e')
        
        self.params['numrel_data'] = decomposed_path
        
        
class Dimmelmeier08(Supernova):
    """
    The Dimmelmeier08 waveform.
    """

    waveform = "Dimmelmeier+08"

    def __init__(self, time, sky_dist=uniform_sky, filepath="signal_s15a2o05_ls.dat", decomposed_path=None, ):
        """

        Parameters
        ----------
        time : float or list 
           The time period over which the injection should be made. If
           a list is given they should be the start and end times, and
           the waveform will be produced at some random point in that
           time range. If a float is given then the injection will be
           made at that specific time.

        sky_dist : func
           The function describing the sky distribution which the injections
           should be made over. Defaults to a uniform sky.

        filepath : str
           The filepath to the numerical relativity waveform.

        decomposed_path : str
           The location where the decomposed waveform file should be stored. Optional.
        """
        
        self._clear_params()
        self.time = time
        self.sky_dist = sky_dist
        if not decomposed_path : decomposed_path = filepath+".dec"
        if not os.path.isfile(decomposed_path) :
            decomposed = self.decompose(filepath, sample_rate = 16384.0, step_back = 0.01, distance = 10e-3)
            
            np.savetxt(decomposed_path, decomposed, header="# time (2,-2) (2,-1) (2,0) (2,1) (2,2)", fmt='%.8e')
        
        self.params['numrel_data'] = decomposed_path
        
    def decompose(self, numrel_file, sample_rate = 16384.0, step_back = 0.01, distance = 10e-3):
        """
        Produce the spherial harmonic decompositions of the Dimmelmeier numerical
        waveform. This is a special case since it is axisymmetric.
        
        Parameters
        ----------
        numrel_file : str
           The location of the numerical relativity waveform file.
        
        sample_rate : float
           The sample rate of the NR file. Defaults to 16384.0 Hz.
        
        step_back : float
           The amount of time, in seconds, of the data which should be included
           before the peak amplitude. Defaults to 0.01 sec.

        distance : float
           The distance, in megaparsecs, from the observer at which the NR waveforms were
           simulated. Defaults to 10 kpc (i.e. 10e-3 Mpc).

        Returns
        -------
        decomposition : ndarray
           The l=2 mode spherical decompositions of the waveform. 
        """

        # Load the times from the file
        data = np.loadtxt(numrel_file)
        data = data.T
        times = data[0]
        times -= times[0]
        
        # Load the hp components   
        strain = data[1]

        # Make the new time vector for the requried sample rate
        target_times = np.arange(0, times[-1], 1.0/sample_rate)

        # Prepare the output matrix
        output = np.zeros((len(target_times), 11))

        # Add the times in to the first column of said matrix
        output[:, 0] = target_times
        #
        # Resample to uniform spacing at 16384 kHz
        #
        strain_new = self.interpolate(times, strain, target_times)
        #
        # Make the output, and rescale it into dimensionless strain values
        #
        output[:,5] = strain_new # * np.sqrt(lal.G_SI / lal.C_SI**4) #/lal.MRSUN_SI / ( distance * lal.PC_SI * 1e6)

        return output
