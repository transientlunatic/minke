import sys
import scipy
from scipy import random
import numpy

import matplotlib.pyplot as plt

from igwn_ligolw import lsctables
from igwn_ligolw.utils import ilwd
#from glue.segments import segment
from ligotimegps import LIGOTimeGPS as GPS

import scipy.interpolate as interp

import os.path

import lal
import lalburst
import lalsimulation

from minke.distribution import *

lalsim = lalsimulation
np = numpy

try:
    import tkinter as tk
except ImportError:
    import matplotlib
    matplotlib.use("agg")



if "numrel_data" in lsctables.SimBurstTable.validcolumns.keys():
    NROK = True
else:
    NROK = False


class Waveform(object):
    """
    Generic container for different source types. 
    Currently, it checks for the waveform type and initializes itself appropriately. 
    In the future, different sources should subclass this and override the generation routines.
    """
    
    table_type = lsctables.SimBurstTable
    sim = lsctables.New(table_type)
    
    numrel_data = []
    waveform = "Generic"
    expnum = 1
    params = {}

    def _clear_params(self):
        self.params = {}
        for a in self.table_type.validcolumns.keys():
            self.params[a] = None
        
    def __getattr__(self, name):
        if name in self.params:
            return self.params[name]
        else:
            raise ValueError(f"The parameter {name} isn't located in this object.")

    def generate_tail(self, sampling=16384.0, length = 1, h_max = 1e-23, h_min = 0):
        """Generate a "low frequency tail" to append to the end of the
        waveform to overcome problems related to memory in the
        waveform.
        
        This code was adapted from an iPython notebook provided by
        Marek Szczepanczyk.

        The tail needs to be added to the waveform after all of the
        other corrections have been applied (DW: I think)

        Parameters
        ----------
        sampling : float
           The sample rate of the injection data. By default this is 16384 Hz, which is the standard advanced LIGO sampling rate.

        length : float
           The length of the tail to be added, in seconds; defaults to 1.

        h_max : float
           The strain at the beginning of the tail -- the strain at the end of the NR data.

        Notes
        -----

        * TODO Confirm that the tail is added-on after the waveform is
        convolved with the antenna pattern.

        """

        times = np.linspace(0, length, length * sampling)
        tail_f = 1.0 / length / 2.0 # Calculate the frequency for a half cosine function over the length of the tail

        tail = 0.5 * (h_max + (h_max-h_min) * np.cos( 2 * np.pi * tail_f * times) + h_min)

        tailout = lal.CreateREAL8Vector(len(tail))
        tailout.data = tail
        
        return tailout
    
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

    def plot(self, figsize=(10,5),):
        """
        Produce a plot of the injection.
        """
        hp, hx, _, _ = self._generate(half=True)
        f, ax = plt.subplots(1,2, figsize=figsize)
        times = np.arange(0, hp.deltaT*len(hp.data.data), hp.deltaT)
        ax[0].plot(times, hp.data.data, label="+ polarisation")
        ax[0].plot(times, hx.data.data, label="x polarisation")
        ax[1].plot(hp.data.data, hx.data.data)
        return f

    def _generate(self, rate=16384.0, half=False, distance=None): 
        """
        Generate the burst described in a given row, so that it can be
        measured.
        
        Parameters 
        ---------- 
        rate : float 
           The sampling rate of the signal, in Hz. 
           Defaults to 16384.0Hz
            
        half : bool 
           Only compute the hp and hx once if this is true;
           these are only required if you need to compute the cross
           products. Defaults to False.

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
        burstobj = self._burstobj()
                
        hp, hx = lalburst.GenerateSimBurst(burstobj, 1.0/rate) 
        if not half :
            hp0, hx0 = lalburst.GenerateSimBurst(burstobj, 1.0/rate) 
        else: 
            hp0, hx0 = hp, hx
        
        return hp, hx, hp0, hx0 

    def _burstobj(self):
        """
        Generate a SimBurst object for this waveform.
        """
        swig_row = self._row()
        burstobj = lalburst.CreateSimBurst()
        
        for a in self.table_type.validcolumns.keys():
            try:
                setattr(burstobj, a, getattr(swig_row,a))
            except AttributeError:
                continue
            except TypeError: 
                continue

        burstobj.waveform = str(self.waveform)

        if NROK:
            if swig_row.numrel_data:
                burstobj.numrel_data = str(swig_row.numrel_data)
            else:
                burstobj.numrel_data = str("")

        return burstobj
    
    def _generate_for_detector(self, ifos, sample_rate = 16384.0, nsamp = 2000):
        data = []
        # Loop through each interferometer
        for ifo in ifos:
            # Make the timeseries
            row = self._row()
            h_resp = lal.CreateREAL8TimeSeries("inj time series", lal.LIGOTimeGPS(0,0), 0, 1.0/sample_rate, lal.StrainUnit, nsamp)
            hp, hx = self._generate(half=True)[:2]
            # Get and apply detector response
            det = lalsimulation.DetectorPrefixToLALDetector(ifo)
            h_tot = lalsimulation.SimDetectorStrainREAL8TimeSeries(hp, hx, row.ra, row.dec, row.psi, det)
            # Inject the waveform into the overall timeseries
            lalsimulation.SimAddInjectionREAL8TimeSeries(h_resp, h_tot, None)
            return h_tot


    def _row(self, sim=None, slide_id=0):
        """
        Produce a simburst table row for this waveform.

        Parameters
        ----------
        sim : table
           The table which the row should be made for.
           If this is left empty the table is assumed to be a 
           sim_burst_table.

        slide_id : int
           The timeslide id. Defaults to 0.
        """
        if not sim: sim = self.sim
        row = sim.RowType()

        for a in self.table_type.validcolumns.keys():
            setattr(row, a, self.params[a])

        if NROK:
            if self.numrel_data:
                row.numrel_data = str(self.numrel_data)
            else:
                row.numrel_data = self.params['numrel_data']
            
        row.waveform = self.waveform
        # Fill in the time
        row.set_time_geocent(GPS(float(self.time)))
        # Get the sky locations
        if not row.ra:
            row.ra, row.dec, row.psi = self.sky_dist()
            row.ra = row.ra[0]
            row.dec = row.dec[0]
            row.psi = row.psi[0]
        row.simulation_id = sim.get_next_id()
        row.waveform_number = random.randint(0,int(2**32)-1)
        ### !! This needs to be updated.
        row.process_id = "process:process_id:0" #procrow.process_id
        row.time_slide_id = ilwd.ilwdchar("time_slide:time_slide_id:%d" % slide_id)

        return row
    
    def interpolate(self, x_old, y_old, x_new, method="linear"):
        """
        Convenience funtion to avoid repeated code
        """
        interpolator = interp.interp1d(x_old, y_old, method)
        return interpolator(x_new)

class StringCusp(Waveform):
    """
    A class to represent a StringCusp injection.
    """
    waveform = "StringCusp"
    
    def __init__(self, amplitude, f_max, time, sky_dist=uniform_sky,):
        """A class to represent a SineGaussian ad-hoc waveform.

        Parameters
        ----------
        amplitude : float
           The amplitude of the injection.

        f_max : float
           The maximum frequency of the injection.

        time : float
           The central time of the injection.

        sky_dist : func
           The function describing the sky distribution which the injections
           should be made over. Defaults to a uniform sky.
        """
        self._clear_params()
        self.sky_dist = sky_dist
        self.params['amplitude'] = amplitude
        self.params['frequency'] = f_max
        self.time = time


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
        self.hrss = self.params['hrss'] = hrss
        self.seed = self.params['seed'] = seed
        self.frequency = self.params['frequency'] = frequency
        self.q = self.params['q'] = q
        self.time = time
        self.polarisation = polarisation
        self.pol_ellipse_e, self.ellipse_angle = self.params['pol_ellipse_e'], self.params['pol_ellipse_angle'] = self.parse_polarisation(self.polarisation)    



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
        hp, hx, _, _ = self._generate(half=True)
        self.params['hrss'] =  lalsimulation.MeasureHrss(hp, hx)


class Numerical2Column(Waveform):
    """
    A superclass to handle ninja-based numerical relativity waveforms.
    """
    waveform = "Numerical" # We shouldn't ever use this anyway
    supernova = False

    extraction = None
    
    def _make_strain(self, sample_rate=16384):
        """
        Calculate the physical strain and time values which correspond to the natural unit
        values in the data file.

        Parameters
        ----------
        distance : float
           The distance, in megaparsec at which the waveform should be produced.
        mass : float
           The total mass, in solar masses of the system to be generated.
        sample_rate : float
           The desired sample rate for the waveform.

        Notes
        -----
        At the moment this only works for files with h+ and hx in the columns.
        """

        data = np.copy(self.data)

        time_scale = (self.total_mass * lal.MSUN_SI * (lal.G_SI / lal.C_SI**3 ))
        mass_geo = (self.total_mass * lal.MSUN_SI * (lal.G_SI / lal.C_SI**2 ))
        distance_geo = (self.distance * 1e6 * lal.PC_SI )#* (lal.C_SI**2/lal.G_SI))
        strain_scale = (distance_geo) / (mass_geo) #(self.total_mass* lal.MSUN_SI)

        if self.extraction:
            strain_scale /= (self.extraction)
        
        data[:,0] *= time_scale
        data[:, 1:] /= strain_scale

        times = data[:,0]

        target_times = np.arange(times[0], times[-1], 1./sample_rate)

        output = np.zeros((len(target_times), 3))
        output[:,0] = target_times
        output[:,1] = self.interpolate(times, data[:,1], target_times)
        output[:,2] = self.interpolate(times, data[:,2], target_times)

        return output
        
    
    def _generate(self, epoch="0.0", rate=16384.0, half=False, tail = True): 
        """
        Generate the burst described in a given row, so that it can be
        measured.
        
        Parameters 
        ---------- 
        rate : float 
           The sampling rate of the signal, in Hz. 
           Defaults to 16384.0Hz
            
        half : bool 
           Only compute the hp and hx once if this is true;
           these are only required if you need to compute the cross
           products. Defaults to False.

        epoch : str 
           The signal injection epoch. 
           This should be given as a string, which will then be 
           split at the decimal to preserve precision.

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

        epoch_sec, epoch_ms = list(map(int, epoch.split(".")))


        data = self._make_strain(rate)
        nsamp = len(data)
        hp = lal.CreateREAL8TimeSeries("inj time series", lal.LIGOTimeGPS(epoch_sec,epoch_ms), 0, 1.0/rate, lal.StrainUnit, nsamp)
        hx = lal.CreateREAL8TimeSeries("inj time series", lal.LIGOTimeGPS(epoch_sec,epoch_ms), 0, 1.0/rate, lal.StrainUnit, nsamp)
        hp.data.data = data[:,1]
        hx.data.data = data[:,2]

        return hp, hx, np.copy(hp), np.copy(hx)
        
    def interpolate(self, x_old, y_old, x_new):
        """
        Convenience funtion to avoid repeated code
        """
        interpolator = interp.interp1d(x_old, y_old)
        return interpolator(x_new)

class Hyperbolic(Numerical2Column):
    def __init__(self, datafile, total_mass, distance, extraction, sky_dist=uniform_sky, **kwargs):
        """
        A class to represent a hyperbolic or parabolic encounter waveform.

        total_mass : float
           The total mass of the system in solar masses.

        distance : float
           The distance, in megaparsecs, at which the waveform should be produced.

        extraction : float
           The extraction radius of the waveform.
        """
        self._clear_params()
        self.data = np.genfromtxt(datafile)
        self.total_mass = total_mass 
        self.distance = distance
        self.sky_dist = sky_dist
        self.extraction = extraction
        self.params.update(kwargs)
    
class Supernova(Waveform):
    """

    A superclass to handle the spherial harmonic decompositions which
    all supernova waveforms require.

    """

    waveform = "Supernova" # We shouldn't ever use this anyway
    supernova = True
    file_distance = 10e-3
    has_memory = False

    def construct_Hlm(self, Ixx, Ixy, Ixz, Iyy, Iyz, Izz, l=2, m=2):
        """
        Construct the expansion parameters Hlm from T1000553.  Returns the expansion
        parameters for l=2, m=m 
        """

        if l!=2:
            print("l!=2 not supported")
            sys.exit()
            if abs(m)>2:
                print("Only l=2 supported, |m| must be <=2")
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

    def _generate(self, rate=16384.0, half=False, distance=None, tail = True): 
        """
        Generate the burst described in a given row, so that it can be
        measured.
        
        Parameters 
        ---------- 
        rate : float 
           The sampling rate of the signal, in Hz. 
           Defaults to 16384.0Hz
            
        half : bool 
           Only compute the hp and hx once if this is true;
           these are only required if you need to compute the cross
           products. Defaults to False.

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
        burstobj = self._burstobj()
                
        hp, hx = lalburst.GenerateSimBurst(burstobj, 1.0/rate) 
        if not half :
            hp0, hx0 = lalburst.GenerateSimBurst(burstobj, 1.0/rate) 
        else: 
            hp0, hx0 = hp, hx


        random.seed(0)
            
        # detrend supernova waveforms
        if hasattr(self, "supernova"):
            hp.data.data, hx.data.data, hp0.data.data, hx0.data.data = scipy.signal.detrend(hp.data.data), scipy.signal.detrend(hx.data.data), scipy.signal.detrend(hp0.data.data), scipy.signal.detrend(hx0.data.data)
            # Rescale for a given distance 
        if burstobj.amplitude: 
            rescale = 1.0 / (self.file_distance / burstobj.amplitude)
            hp.data.data, hx.data.data, hp0.data.data, hx0.data.data = hp.data.data * rescale, hx.data.data * rescale, hp0.data.data * rescale, hx0.data.data * rescale

            if self.has_memory and tail:
                # Apply the tail correction for memory
                tail_hp = self.generate_tail(length = 1, h_max = hp.data.data[-1], h_min = hp.data.data[0])
                tail_hx = self.generate_tail(length = 1, h_max = hx.data.data[-1], h_min = hx.data.data[0])

                hp_data = np.append(hp.data.data,tail_hp.data)
                hx_data = np.append(hx.data.data,tail_hx.data)

                del tail_hp, tail_hx
                
                tail_hp = lal.CreateREAL8Vector(len(hp_data))
                tail_hp.data = hp_data
                tail_hx = lal.CreateREAL8Vector(len(hx_data))
                tail_hx.data = hx_data

                hp.data = tail_hp
                hx.data = tail_hx

                del tail_hp, tail_hx
        
        return hp, hx, hp0, hx0 
    
    def generate_tail(self, sampling=16384.0, length = 1, h_max = 1e-23, h_min = 0):
        """Generate a "low frequency tail" to append to the end of the
        waveform to overcome problems related to memory in the
        waveform.
        
        This code was adapted from an iPython notebook provided by
        Marek Szczepanczyk.

        The tail needs to be added to the waveform after all of the
        other corrections have been applied (DW: I think)

        Parameters
        ----------
        sampling : float
           The sample rate of the injection data. By default this is 16384 Hz, which is the standard advanced LIGO sampling rate.

        length : float
           The length of the tail to be added, in seconds; defaults to 1.

        h_max : float
           The strain at the beginning of the tail -- the strain at the end of the NR data.

        Notes
        -----

        * TODO Confirm that the tail is added-on after the waveform is
        convolved with the antenna pattern.

        """

        times = np.linspace(0, length, length * sampling)
        tail_f = 1.0 / length / 2.0 # Calculate the frequency for a half cosine function over the length of the tail

        tail = 0.5 * (h_max + (h_max-h_min) * np.cos( 2 * np.pi * tail_f * times) + h_min)

        tailout = lal.CreateREAL8Vector(len(tail))
        tailout.data = tail
        
        return tailout
        
        
    def interpolate(self, x_old, y_old, x_new):
        """
        Convenience funtion to avoid repeated code
        """
        interpolator = interp.interp1d(x_old, y_old)
        return interpolator(x_new)

    def decompose(self, numrel_file, sample_rate = 16384.0, step_back = 0.01):
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
        target_times = np.arange(times[0], times[-1], 1.0/sample_rate)

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

class Ott2013(Supernova):
    """
    The Ott+2013 supernova waveform
    """
    has_memory = True
    waveform = "Ott+13"
    def __init__(self, theta, phi, time, sky_dist=uniform_sky, distance = 10e-3,  filepath=None, family="s27fheat1p05", decomposed_path=None):
        """

        Parameters
        ----------
        phi : float
           The internal phi parameter of the supernova injection.
        
        theta : float
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

        distance : float 
           The distance, in megaparsecs, at which the injection should be made.  

        filepath : str
           The filepath to the folder containing the pre-rotated numerical relativity waveforms.

        family : str
           The family of waveforms which are to be used for the injection set.

        decomposed_path : str
           The location where the decomposed waveform file should be stored. Optional.
        """
        
        self._clear_params()
        self.time = time
        self.params['phi'] = phi
        self.params['incl'] = theta
        self.sky_dist = sky_dist
        #self.params['numrel_data'] = filepath #decomposed_path #self.numrel_data
        if not decomposed_path : decomposed_path = filepath+".dec"
        if not os.path.isfile(decomposed_path) :
            decomposed = self.decompose(filepath, sample_rate = 16384.0, step_back = 0.01)
            np.savetxt(decomposed_path, decomposed, header="time (2,-2) (2,-1) (2,0) (2,1) (2,2)", fmt='%.8e')
        self.numrel_data = self.params['numrel_data'] = decomposed_path
        self.params['amplitude'] = distance # We store the distance in the amplitude column because there isn't a distance column
        self.params['hrss'] = self.file_distance # Again the hrss value is the distance at which the files are scaled


class Mueller2012(Supernova):
    """
    The Mueller2012 waveform.
    """

    waveform = "Mueller+12"
    has_memory = True
    
    def __init__(self, theta, phi, time, distance = 10e-3, sky_dist=uniform_sky, filepath=None, family="L15-3", decomposed_path=None):
        """

        Parameters
        ----------
        phi : float
           The internal phi parameter of the supernova injection.
        
        theta : float
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
           The filepath to the folder containing the pre-rotated numerical relativity waveforms.

        family : str
           The family of waveforms which are to be used for the injection set.

        decomposed_path : str
           The location where the decomposed waveform file should be stored. Optional.
        """
        
        self._clear_params()
        self.time = time
        self.params['phi'] = phi
        self.params['incl'] = theta
        self.sky_dist = sky_dist
        if not decomposed_path : decomposed_path = filepath+".dec"
        if not os.path.isfile(decomposed_path) :
            decomposed = self.decompose(filepath, sample_rate = 16384.0, step_back = 0.01)
            np.savetxt(decomposed_path, decomposed, header="time (2,-2) (2,-1) (2,0) (2,1) (2,2)", fmt='%.8e')
        
        #self.numrel_data = filepath + "/" + family
        self.params['numrel_data'] = decomposed_path #self.numrel_data
        self.params['amplitude'] = distance # We store the distance in the amplitude column because there isn't a distance column
        self.params['hrss'] = self.file_distance # Again the hrss value is the distance at which the files are scaled

    def decompose(self, numrel_file, sample_rate = 16384.0, step_back = 0.01):
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
        times = data[1]
        times -= times[0]

        # Load the I components from the file        
        Ixx, Iyy, Izz, Ixy, Ixz, Iyz = data[6:]

        # Make the new time vector for the requried sample rate
        target_times = np.arange(times[0], times[-1], 1.0/sample_rate)

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



    # def _generate(self):
    #     """

    #     Generate the Mueller waveforms. This must be performed
    #     differently to other waveform morphologies, since we require
    #     the use of pre-generated text files.

    #     The filepath and the start of the filenames should be provided in
    #     the numrel_data column of the SimBurstTable, so we need to contruct
    #     the rest of the filename from the theta and phi angles, and then load 
    #     that file.

    #     """
    #     theta, phi = self.params['incl'], self.params['phi']
    #     numrel_file_hp = self.numrel_data + "_costheta{:.3f}_phi{:.3f}-plus.txt".format(theta, phi)
    #     numrel_file_hx = self.numrel_data + "_costheta{:.3f}_phi{:.3f}-cross.txt".format(theta, phi)

    #     data_hp = np.loadtxt(numrel_file_hp)
    #     data_hx = np.loadtxt(numrel_file_hx)

    #     return data_hp, data_hx, data_hp, data_hx

class Scheidegger2010(Supernova):
    """
    The Scheidegger2010 waveform.
    """

    waveform = "Scheidegger+10"

    def __init__(self, theta, phi, time, distance = 10e-3, sky_dist=uniform_sky, filepath=None, family="R1E1CA_L", decomposed_path=None):
        """

        Parameters
        ----------
        phi : float
           The internal phi parameter of the supernova injection.
        
        theta : float
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
           The filepath to the folder containing the pre-rotated numerical relativity waveforms.

        family : str
           The family of waveforms which are to be used for the injection set.

        decomposed_path : str
           The location where the decomposed waveform file should be stored. Optional.
        """
        
        self._clear_params()
        self.time = time
        self.params['phi'] = phi
        self.params['incl'] = theta
        self.sky_dist = sky_dist


        if not decomposed_path : decomposed_path = filepath+".dec"
        if not os.path.isfile(decomposed_path) :
            decomposed = self.decompose(filepath, sample_rate = 16384.0, step_back = 0.01)
            np.savetxt(decomposed_path, decomposed, header="time (2,-2) (2,-1) (2,0) (2,1) (2,2)", fmt='%.8e')
        
        #self.numrel_data = filepath + "/" + family
        self.params['numrel_data'] = decomposed_path #self.numrel_data
        self.params['amplitude'] = distance # We store the distance in the amplitude column because there isn't a distance column
        self.params['hrss'] = self.file_distance # Again the hrss value is the distance at which the files are scaled
        
class Dimmelmeier08(Supernova):
    """
    The Dimmelmeier08 waveform.
    """

    waveform = "Dimmelmeier+08"

    def __init__(self, time, distance = 10e-3, sky_dist=uniform_sky, filepath="signal_s15a2o05_ls.dat", decomposed_path=None, ):
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
            decomposed = self.decompose(filepath, sample_rate = 16384.0, step_back = 0.01)
            np.savetxt(decomposed_path, decomposed, header="time (2,-2) (2,-1) (2,0) (2,1) (2,2)", fmt='%.8e')
        self.params['phi']=0
        self.params['incl']=90
        self.params['numrel_data'] = decomposed_path#
        self.params['amplitude'] = distance # We store the distance in the amplitude column because there isn't a distance column
        self.params['hrss'] = self.file_distance # Again the hrss value is the distance at which the files are scaled
        
    def decompose(self, numrel_file, sample_rate = 16384.0, step_back = 0.01):
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

        Returns
        -------
        decomposition : ndarray
           The l=2 mode spherical decompositions of the waveform. 
        """
        extract_dist = 10e-3
        # Load the times from the file
        data = np.loadtxt(numrel_file)
        data = data.T
        times = data[0]*1e-3
        times -= times[0]
        
        # Load the hp components   
        strain = data[1]
        # Make the new time vector for the requried sample rate
        target_times = np.arange(times[0], times[-1], 1.0/sample_rate)

        # Prepare the output matrix
        output = np.zeros((len(target_times), 11))

        # Add the times in to the first column of said matrix
        output[:, 0] = target_times #/ lal.MTSUN_SI
        #
        # Resample to uniform spacing at 16384 kHz
        #
        strain_new = self.interpolate(times, strain, target_times)
        #
        # Make the output, and rescale it into dimensionless strain values
        #
        output[:,5] = strain_new #/*  ( extract_dist * lal.PC_SI * 1.0e6) 

        return output


class Ringdown(Waveform):
    """
    A class to handle Rindown waveforms.
    """
    table_type = lsctables.SimRingdownTable
    waveform = "GenericRingdown"

        
class Yakunin10(Supernova):
    """
    The Yakunin10 waveform.
    """

    waveform = "Yakunin+10"
    
    def __init__(self, time, distance = 10e-3, sky_dist=uniform_sky, filepath="Yakunin2010/hplus-B12-WH07_tail.txt", decomposed_path=None, ):
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
        self.params['amplitude'] = distance # We store the distance in the amplitude column because there isn't a distance column
        self.params['hrss'] = self.file_distance # Again the hrss value is the distance at which the files are scaled
        
        self.time = time
        self.sky_dist = sky_dist
        if not decomposed_path : decomposed_path = filepath+".dec"
        if not os.path.isfile(decomposed_path) :
            decomposed = self.decompose(filepath, sample_rate = 16384.0, step_back = 0.01)
            np.savetxt(decomposed_path, decomposed, header="time (2,-2) (2,-1) (2,0) (2,1) (2,2)", fmt='%.8e')
        self.params['phi']=0
        self.params['incl']=90
        self.params['numrel_data'] = decomposed_path
        
    def decompose(self, numrel_file, sample_rate = 16384.0, step_back = 0.01):
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

        Returns
        -------
        decomposition : ndarray
           The l=2 mode spherical decompositions of the waveform. 
        """
        extract_dist = 10e-3
        # Load the times from the file
        data = np.loadtxt(numrel_file)
        data = data.T
        times = data[0]
        times -= times[0]
        
        # Load the hp components   
        strain = data[1]
        # Make the new time vector for the requried sample rate
        target_times = np.arange(times[0], times[-1], 1.0/sample_rate)

        # Prepare the output matrix
        output = np.zeros((len(target_times), 11))

        # Add the times in to the first column of said matrix
        output[:, 0] = target_times #/ lal.MTSUN_SI
        #
        # Resample to uniform spacing at 16384 kHz
        #
        strain_new = self.interpolate(times, strain, target_times)
        #
        # Make the output, and rescale it into dimensionless strain values
        #
        output[:,5] = strain_new #/*  ( extract_dist * lal.PC_SI * 1.0e6) 

        return output

    
class LongDuration(Supernova):
    """

    A superclass to handle the spherial harmonic decompositions which
    long duration numerical relativity bursts may require.

    """

    waveform = "LongDuration" # We shouldn't ever use this anyway
    supernova = True

class ADI(LongDuration):
    """
    Accretion disk instability waveforms which are generated using the method described in 
    LIGO-T1100093, at https://dcc.ligo.org/LIGO-T1100093. The waveforms are based off a model
    by MH van Putten,
       M. H. van Putten, A. Levinson, H. K. Lee, T. Regimbau, M. Punturo, and G. M. Harry. Phys. Rev. D., 69(4), 044007, 2004.
       M. H. van Putten. Phys. Rev. Lett., 87(9), 091101, 2001.
    The waveforms are stored in .mat binary files which can be read-in by SciPy.


    """

    waveform = "ADI"

    def __init__(self, time, sky_dist=uniform_sky, filepath="stamp_adi_a_tapered.mat", decomposed_path=None, ):
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
            decomposed = self.decompose(filepath)
            np.savetxt(decomposed_path, decomposed, header="time\thplus\thcross", fmt='%.8e')
        #decomposed_path = filepath
        self.params['phi']=0
        self.params['incl']=90
        self.params['numrel_data'] = decomposed_path

    def _generate(self, rate = 16384.0, half=False):
        data = np.genfromtxt(self.params['numrel_data'])
        nsamp = len(data)
        hp = lal.CreateREAL8TimeSeries("inj time series", lal.LIGOTimeGPS(0,0), 0, 1.0/rate, lal.StrainUnit, nsamp)
        hx = lal.CreateREAL8TimeSeries("inj time series", lal.LIGOTimeGPS(0,0), 0, 1.0/rate, lal.StrainUnit, nsamp)

        hp.data.data = data[:,1]
        hx.data.data = data[:,2]

        return hp, hx, np.copy(hp), np.copy(hx)
        
    def decompose(self, numrel_file, sample_rate = 16384.0, step_back = 0.01):
        """
        Produce the spherial harmonic decompositions of the ADI
        waveform. This is a special case since it is axisymmetric.
        
        Parameters
        ----------
        numrel_file : str
           The location of the numerical relativity waveform file.
        
        sample_rate : float
           The sample rate of the output NR file. Defaults to 16384.0 Hz, and should
           be the same as the data rate of the detector. 
        
        step_back : float
           The amount of time, in seconds, of the data which should be included
           before the peak amplitude. Defaults to 0.01 sec.

        Returns
        -------
        decomposition : ndarray
           The re-interpolated file at the desired sample rate which is in the 
           <time hp hc> format which can be accepted by LALSimulation. 
        """
        from scipy import io
        # Load the matlab file
        data = io.matlab.loadmat(numrel_file)

        comment = data['comment'][0].split(';')
        comment_dict = {}
        for line in comment:
            sp = line.split("=")
            comment_dict[sp[0].strip()] = sp[1].strip()

        extract_dist = comment_dict['dist']
        # We actually want the extract distance as a float of megaparsecs
        if extract_dist == "1 Mpc": extract_dist = 1.0

        # Load the sample rate of the file from the file
        fs = int(data['fs'])
        # Determine the end time
        start = 0
        end = len(data['hp']) * 1.0 / fs
        # Make the time array
        times = np.arange(start, end, 1.0/fs)
        # Make the new time vector for the requried sample rate
        target_times = np.arange(times[0], times[-1], 1.0/sample_rate)

        #print len(target_times)
        # Load the hp components   
        strainp = data['hp'].T[0].astype(np.float32)
        strainc = data['hc'].T[0].astype(np.float32)
        #del data

        # Prepare the output matrix
        output = np.zeros((len(target_times), 3))

        # Add the times in to the first column of said matrix
        output[:, 0] = target_times
        #
        # Resample to uniform spacing at 16384 kHz
        #
        strainp_new = self.interpolate(times, strainp, target_times)
        strainc_new = self.interpolate(times, strainc, target_times)
        #
        # Make the output.
        #
        output[:,1] = strainp_new
        output[:,2] = strainc_new

        return output

class BBHRingdown(Ringdown):
    """
    A class to represent BBH ringdowns.
    """
    #lalsimfunction = SimBlackHoleRingdown
    waveform = "BBHRingdown"
    def __init__(self, time, phi0, mass, spin, massloss, distance, inclination, sky_dist=uniform_sky):
        """
        Binary Black Hole (BBH) Ringdown waveform

        Parameters
        ----------
        time : float
           The time that the waveform should be generated at, in gps seconds.
        phi0 : float
           The starting phase.
        mass : float
           The mass of the final black hole in solar masses.
        spin : float
           The dimensionless spin parameter for the final black hole.
        massloss : float
           The total mass loss of the system. (Also denoted epsilon).
        distance : float
           The effective luminosity distance at which the signal should be generated.
        inlination : float
            The inclination of the system in degrees.
        """
        self._clear_params()
        self.time = self.geocent_start_time = self.params['geocent_start_time'] =  time
        self.sky_dist = sky_dist
        self.params['simulation_id'] = self.simulation_id =  self.sim.get_next_id()
        self.params['phase'] = phi0
        self.params['mass'] = mass # in solar masses
        self.params['spin'] = spin
        self.params['epsilon'] = massloss
        self.params['eff_dist_l'] = self.eff_dist_l = distance # megaparsec
        self.params['inclination'] = self.inclination = float(inclination)

    def _generate(self, rate=16384.0, half=False, l = 2, m = 2):
        """
        Generate this BBH Ringdown waveform.

        Parameters
        ----------
        rate : float
           The signal sampling rate. Defaults to 16384.0 Hz.
        l : int
           The azimuthal number of the mode to be generated.
        m : int 
           The polar number of the mode to be generated.    
        half : bool 
           Only compute the hp and hx once if this is true;
           these are only required if you need to compute the cross
           products. Defaults to False.

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
        dt = 1.0 / rate
        hp, hx = lalsimulation.SimBlackHoleRingdown(self.params['geocent_start_time'],
                                                    self.params["phase"],
                                                    dt,
                                                    self.params['mass']*lal.MSUN_SI,
                                                    self.params['spin'],
                                                    self.params['epsilon'],
                                                    self.params['eff_dist_l'] *  1e6 * lal.PC_SI,
                                                    np.deg2rad(self.params['inclination']),
                                                    l, m)
        if not half:
            hp0, hx0 = lalsimulation.SimBlackHoleRingdown(self.params['geocent_start_time'],
                                                          self.params["phase"],
                                                          dt,
                                                          self.params['mass']*lal.MSUN_SI,
                                                          self.params['spin'],
                                                          self.params['epsilon'],
                                                          self.params['eff_dist_l'] *  1e6 * lal.PC_SI,
                                                          np.deg2rad(self.params['inclination']),
                                                          l, m)
        else:
            hp0, hx0 = hp, hx
        
        return hp, hx, hp0, hx0
