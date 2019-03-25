"""
88b           d88  88               88                    
888b         d888  ""               88                    
88`8b       d8'88                   88                    
88 `8b     d8' 88  88  8b,dPPYba,   88   ,d8   ,adPPYba,  
88  `8b   d8'  88  88  88P'   `"8a  88 ,a8"   a8P_____88  
88   `8b d8'   88  88  88       88  8888[     8PP"""""""  
88    `888'    88  88  88       88  88`"Yba,  "8b,   ,aa  
88     `8'     88  88  88       88  88   `Y8a  `"Ybbd8"'  
                                                          
--------------------------------------------------------

This file is a part of Minke, a tool for generating simulated
gravitational wave signals, used for characterising and training
search algorithms.

Minke was created by Daniel Williams, based on work started by Chris
Pankow and others, and is built around the LALSimulation library.



"""
from glue.ligolw import ligolw, utils, lsctables
lsctables.use_in(ligolw.LIGOLWContentHandler);
import numpy
import lalburst, lalsimulation, lalmetaio
from minke.antenna import response

from lal import TimeDelayFromEarthCenter as XLALTimeDelayFromEarthCenter
#from pylal.xlal.datatypes.ligotimegps import LIGOTimeGPS
from lal import LIGOTimeGPS
from glue.ligolw.utils import process
import glue

import glue.ligolw
import gzip 

import lal, lalframe

import numpy as np
import pandas as pd
import os
import os.path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import re
import random
import minke
from minke import sources
sourcemap = {}
for classin in dir(sources):
    classin = sources.__dict__[classin]
    if hasattr(classin, "waveform"):
        sourcemap[classin.waveform] = classin
        
def source_from_row(row):
    waveform = row.waveform
    sourceobj = sourcemap[row.waveform].__new__(sourcemap[row.waveform])
    sourceobj.numrel_data = str("")
    params = {}
    for attr in dir(row):
        if not attr[0] == "_" and not attr[:3] =="get":
            #print attr
            try:
                params[attr] = getattr(row, attr)
                setattr(sourceobj, attr, getattr(row, attr))
            except AttributeError:
                print("Error processing the {} column".format(attr))
    sourceobj.params = params
    try:
        sourceobj.time = row.time_geocent_gps
    except:
        sourceobj.time = row.geocent_start_time
        pass
    return sourceobj

def source_from_dict(params):
    sourceobj = sourcemap[params['morphology']].__new__(sourcemap[params['morphology']])
    sourceobj.numrel_data = str("")
    params = {}
    for attr in dir(row):
        if not attr[0] == "_" and not attr[:3] =="get":
            #print attr
            params[attr] = getattr(row, attr)
            setattr(sourceobj, attr, getattr(row, attr))
    sourceobj.params = params
    try:
        sourceobj.time = row.time_geocent_gps
    except:
        sourceobj.time = row.geocent_start_time
        pass
    return sourceobj


table_types = {
    # Ad-Hoc
    "ga" : lsctables.SimBurstTable,
    "sg" : lsctables.SimBurstTable,
    "wnb" : lsctables.SimBurstTable,
    "sc" : lsctables.SimBurstTable,
    # Supernova Families
    "d08" : lsctables.SimBurstTable,
    "s10" : lsctables.SimBurstTable,
    "m12" : lsctables.SimBurstTable,
    "o13" : lsctables.SimBurstTable,
    "y10" : lsctables.SimBurstTable,
    # Long Duration
    "adi" : lsctables.SimBurstTable,
    # Ringdown
    "rng" : lsctables.SimRingdownTable,
    "gng" : lsctables.SimRingdownTable,
    }

tables = {
    "burst" : lsctables.SimBurstTable,
    "ringdown" : lsctables.SimRingdownTable
    }

def mkdir(path):
    """
    Make all of the tree of directories in a given path if they don't
    already exist.

    Parameters
    ----------
    path : str
       The path to the desired directory.

    """
    sub_path = os.path.dirname(path)
    if not os.path.exists(sub_path):
        mkdir(sub_path)
    if not os.path.exists(path):
        os.mkdir(path)


class TableTypeError(Exception):
    pass
        
class MDCSet():

    inj_families_names = {'ga' : 'Gaussian',
                          'sg' : 'SineGaussian',
                          'wnb': 'BTLWNB',
                          "sc" : "StringCusp",
                          # Supernova families
                          'd08' : 'Dimmelmeier+08',
                          's10' : 'Scheidegger+10',
                          'm12' : 'Mueller+12',
                          'o13' : 'Ott+13',
                          'y10' : "Yakunin+10",
                          # Long-duration
                          'adi' : 'ADI',
                          # Ringdown
                          'rng' : "BBHRingdown",
                          'gng' : "GenericRingdown",
                          }

    inj_families_abb = dict((v,k) for k,v in inj_families_names.iteritems())

    hist_parameters = {
        "StringCusp": ["amplitude", "ra", "dec"],
        "SineGaussian": ["hrss", "psi", "ra", "dec"],
        "Gaussian": ["hrss", "psi", "ra", "dec"],
        "BTLWNB": ["hrss", "ra", "dec"],
        "Dimmelmeier+08": ['hrss', 'ra', 'dec']
    }

    waveforms = []

    def __init__(self, detectors, name='MDC Set', table_type = "burst"):
        """
        Represents an MDC set, stored in an XML SimBurstTable file.
        
        Parameters
        ----------
        detectors : list 
            A list of detector names where the injections should be made.

        name : str
            A name for the MDC Set. Defaults to 'MDC Set'.

        table_type : str
            The type of table which should be generated. Default is `burst`, 
            which generates a SimBurstTable.
        """
        self.detectors = detectors
        
        self.waveforms = []
        self.strains = []
        self.egw = []
        self.times = []
        self.name = name
            
        self.times = np.array(self.times)

        self.table_type = tables[table_type]

    def __add__(self, waveform):
        """
        Handle a waveform being added to the MDC set.

        Parameters
        ----------
        waveform : Waveform object
           The waveform which should be added to the MDC set.

        """

        # Check that this type of waveform can go into this type of
        # XML file.
        if not table_types[self.inj_families_abb[waveform.waveform]] == self.table_type:
            raise TableTypeError()
        
        self.waveforms.append(waveform)
        self.times = np.append(self.times, waveform.time)

    def save_xml(self, filename):
        """
        Save the MDC set as an XML SimBurstTable.

        Parameters
        ----------
        filename : str
           The location to save the xml file. The output is gzipped, so ending it with 
           a ".gz" would stick with convention.
        """
        xmldoc = ligolw.Document()
        lw = xmldoc.appendChild(ligolw.LIGO_LW())
        sim = lsctables.New(self.table_type)
        lw.appendChild(sim)
        # This needs to be given the proper metadata once the package has the maturity to
        # write something sensible.
        for waveform in self.waveforms:
            procrow = process.register_to_xmldoc(xmldoc, "minke_burst_mdc+{}".format(minke.__version__), {}) # waveform.params)
            try:
                waveform_row = waveform._row(sim)
                waveform_row.process_id = procrow.process_id
            except:
                row = sim.RowType()
                for a in self.table_type.validcolumns.keys():
                    if a in waveform.params.keys():
                        setattr(row, a, waveform.params[a])
                    else:
                        if not hasattr(waveform, a):
                            setattr(row, a, 0)
                        else:
                            setattr(row, a, getattr(waveform, a))

                row.waveform = waveform.waveform
                if self.table_type == lsctables.SimBurstTable:
                    # Fill in the time
                    row.set_time_geocent(GPS(float(waveform.time)))
                    # Get the sky locations
                    row.ra, row.dec, row.psi = waveform.ra, waveform.dec, waveform.psi
                row.simulation_id = waveform.simulation_id
                row.waveform_number = random.randint(0,int(2**32)-1)
                ### !! This needs to be updated.
                row.process_id = "process:process_id:0" #procrow.process_id

                waveform_row = row
            
            sim.append(waveform_row)
            #del waveform_row
        # Write out the xml and gzip it.
        utils.write_filename(xmldoc, filename, gz=True)

    def load_xml(self, filename, full=True, start=None, stop=None):
        """Load the MDC Set from an XML file containing the SimBurstTable.

        Parameters
        ----------
        filename : str
           The filename of the XML file.

        full : bool 
           If this is true (which is the default) then all of
           the calculated parameters are computed from the waveform
           definintion.

        start : float 
           The time at which the xml read-in should
           start. The default is "None", in which case the xml file
           will be read-in from the start.

        end : float
           The last time to be read from the xml file. The default is None, 
           which causes the xml to be read right-up to the last time in the 
           file.

        To Do
        -----
        A the moment this loads the information in to the object, but it 
        doesn't produce waveform objects for each of the injections in the
        file. This should be fixed so that the object works symmetrically.
        """
        i = 0
        #sim_burst_table = lalburst.SimBurstTableFromLIGOLw(filename, start, stop)

        xml = glue.ligolw.utils.load_filename(filename, 
                                              contenthandler = glue.ligolw.ligolw.LIGOLWContentHandler,
                                              verbose = True)
        sim_burst_table = glue.ligolw.table.get_table(xml, self.table_type.tableName)
        
        for i,simrow in enumerate(sim_burst_table):
            # This is an ugly kludge to get around the poor choice of wavform name in the xmls, and
            if simrow.waveform[:3]=="s15": 
                self.numrel_file = str(sim_burst_table.waveform)
                sim_burst_table.waveform = "Dimmelmeier+08"

            self.waveforms.append(source_from_row(simrow))
            
            if full:
                self._measure_hrss(i)
                self._measure_egw_rsq(i)

            if self.table_type == tables["burst"]:
                self.times = np.append(self.times, float(simrow.time_geocent_gps))
            
    def _generate_burst(self,row,rate=16384.0):
        """
        Generate the burst described in a given row, so that it can be 
        measured.
        
        Parameters
        ----------
        row : SimBurst Row
            The row of the waveform to be measured
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
        row = self.waveforms[row]
        hp, hx, hp0, hx0 = row._generate()
        return hp, hx, hp0, hx0
    
    def _getDetector(self, det):
        """
        A method to return a LALDetector object corresponding to a detector's
        X#-style name, e.g. 'H1' as the Hanford 4km detector.
        
        Parameters
        ----------
        det : str
            A string describing the detector in the format letter-number, e.g
            "H1" would be the Hanford 4km detector, "L1" would be the 
            Livingston 4km, and so-forth.
            
        Returns
        -------
        detector : LALDetector
            The LAL object describing the detector
        """
        # get detector
        return lalsimulation.DetectorPrefixToLALDetector(det)
        #if det not in lal.cached_detector_by_prefix.keys(): 
        #      raise ValueError, "%s is not a cached detector.  "\
        #            "Cached detectors are: %s" % (det, inject.cached_detector.keys()) 
        #return lal.cached_detector_by_prefix[det]

    def _timeDelayFromGeocenter(self, detector, ra, dec, gpstime):
        """
        Calculate the time delay between the geocentre and a given detector
        for a signal from some sky location.
        
        Parameters
        ----------
        detector : str
            A string describing the detector, e.g. H1 is the Hanford 4km 
            detector.
        ra : float
            The right-ascension of the observation in radians
        dec : float
            The declination of the obser
        """
        if isinstance(detector, str): detector = self._getDetector(detector)
        gpstime = LIGOTimeGPS(float(gpstime))
        return XLALTimeDelayFromEarthCenter(detector.location, ra, dec, gpstime)
    
    def directory_path(self):
        """
        Generate the directory where the frames from this MDC should be stored, 
        so, e.g. Gaussians 0d100 would go in "ga/ga0d100/"

        Returns
        -------
        str 
           the folder structure
        """
        name = self._simID(0)
        abb = self.inj_families_abb[self.waveforms[0].waveform].lower()
        return "{}/{}".format(abb, name)
        
        
    def _simID(self, row):
        """
        Generate a name for an injection set in the format expected by cWB
        
        Parameters
        ----------
        row : SimBurst
            The simburst table row describing the injection

        Returns
        -------
        str
           The name of the injection in the cWB format
        """
        row = self.waveforms[row]
        
        name = ''
        numberspart = ''
        if row.waveform in ("Dimmelmeier+08", "Scheidegger+10", "Mueller+12", "Ott+13", "Yakunin+10"):
            #print row
            numberspart = os.path.basename(row.params['numrel_data']).split('.')[0]

        if row.waveform == "Gaussian":
            numberspart = "{:.3f}".format(row.duration * 1e3)
        elif row.waveform == "SineGaussian":
            if row.pol_ellipse_e==1.0: 
                pol="linear"
            elif row.pol_ellipse_e==0.0:
                pol="circular"
            elif 0.0<row.pol_ellipse_e<1.0: 
                pol = "elliptical"
            else:
                pol = "inclined"
            numberspart = "f{:.0f}_q{:.0f}_{}".format(row.frequency, row.q, pol)
        elif row.waveform == "BTLWNB":
            numberspart = "{}b{}tau{}".format(row.frequency, row.bandwidth, row.duration)
        name += '{}_{}'.format(self.inj_families_abb[row.waveform].lower(), numberspart).replace('.','d')

        return name
    
    def _measure_hrss(self, row, rate=16384.0):
        """
        Measure the various components of hrss (h+^2, hx^2, hphx) for a given 
        input row. This is accomplished by generating the burst and calling 
        the SWIG wrapped  XLALMeasureHrss in lalsimulation. 
        
        Parameters
        ----------
        row : int
            The row number of the waveforms to be measured
        rate : float
            The sampling rate of the signal, in Hz. Defaults to 16384.0Hz
            
        Returns
        -------
        hrss : float
            The measured hrss of the waveform amplitude: sqrt(|Hp|^2 + |Hx|^2)
        hphp : float
            The hrss of the + polarisation only.
        hxhx : float
            The hrss of the x polarisation only.
        hphx : float
            The hrss of |HpHx| 
        """
        row = self.waveforms[row]
        hp, hx, hp0, hx0 = row._generate() #self._generate_burst(row)# self.hp, self.hx, self.hp0, self.hx0

        hp0.data.data *= 0
        hx0.data.data *= 0
        
        # H+ hrss only
        hphp = lalsimulation.MeasureHrss(hp, hx0)**2
        # Hx hrss only
        hxhx = lalsimulation.MeasureHrss(hp0, hx)**2
        # sqrt(|Hp|^2 + |Hx|^2)
        hrss = lalsimulation.MeasureHrss(hp, hx)

        hp.data.data = numpy.abs(hx.data.data) + numpy.abs(hp.data.data)
        # |H+Hx|
        hphx = (lalsimulation.MeasureHrss(hp, hx0)**2 - hrss**2)/2
        #print hrss
        self.strains.append([hrss, hphp, hxhx, hphx])
    
    def _measure_egw_rsq(self, row,  rate=16384.0):
        """
        Measure the energy emitted in gravitational waves divided 
        by the distance squared in M_solar / pc^2. This is accomplished 
        by generating the burst and calling the SWIG wrapped  
        XLALMeasureHrss in lalsimulation. 
        
        Parameters
        ----------
        row : int
            The row number of the waveforms to be measured
        rate : float
            The sampling rate of the signal, in Hz. Defaults to 16384.0Hz
            
        Returns 
        -------
        egw : float
            The energy emitted in gravitational waves divided 
            by the distance squared in M_solar / pc^2.
        """
        hp, hx, _, _ =  self._generate_burst(row)
        self.egw.append(lalsimulation.MeasureEoverRsquared(hp, hx))
    
    def _responses(self, row):
        """
        Calculate the antenna repsonses for each detector to the waveform.
        
        Parameters
        ----------
        row : int
            The row number of the waveforms to be measured
            
        Returns
        -------
        responses : list of lists
            A list containing the lists of antenna responses, with the first 
            element of each list containing the detector acronym.
        """
        output = []
        row = self.waveforms[row]
        for detector in self.detectors:
            time = row.time_geocent_gps + self._timeDelayFromGeocenter(detector, row.ra, row.dec, row.time_geocent_gps)
            time = np.float64(time)
            rs = response(time, row.ra, row.dec, 0, row.psi, 'radians', detector)
            output.append([detector, time, rs[0], rs[1]]   )
        return output
    
    def plot_skymap(self):
        """
        Plot a skymap of the injections distribution in RA and DEC on a Hammer projection.
        
        Returns
        -------
        matplotlib figure
        """
        fig = plt.figure()
        # Load the ra and dec numbers out of the waveforms 
        dec = [getattr(s, 'dec') for s in self.waveforms]
        ra = [getattr(s, 'ra') for s in self.waveforms]
        
        # Make the plot on a hammer projection
        plt.subplot(111, projection='hammer')
        H, x, y = np.histogram2d(ra, dec, [50, 25], range=[[0, 2*np.pi], [-np.pi/2, np.pi/2]])
        dist = plt.pcolormesh(x-np.pi,y, H.T, cmap="viridis")
        plt.title("Sky distribution")
        plt.colorbar(dist, orientation='horizontal')
        return fig
    
    def plot_hist(self, parameter):
        """
        Plot a histogram of a waveform parameter.

        Parameters
        ----------
        parameter : str
           The name of the simburst table parameter which is desired for the plot.

        Returns
        -------
        matplotlib figure
        """           
        fig = plt.figure()
        prms = [getattr(s, parameter) for s in self.waveforms]
        ax2 = plt.subplot(111)
        ax2.set_title("{} distribution".format(parameter))
        ax2.set_xlabel(parameter)
        ax2.hist(prms, bins=100, log=True, histtype="stepfilled", alpha=0.6);
        return fig

    def gravEn_row(self, row, frame):
        """
        Produces a gravEn-style log row for a row of the simBurstTable.
        
        Parameters
        ----------
        row : int
            The row number of the waveforms to be measured
            
        Returns
        -------
        str
            A string in the gravEn format which describes the injection.
        """
        strains = self.strains[row]
        rowname = self._simID(row)
        responses = self._responses(row)
        energy = self.egw[row]
        row = self.waveforms[row]
        output = []
        if not row.incl:
            cosincl = ""
        else:
            cosincl = np.cos(row.incl)

        output.append(self.name)                  # GravEn_SimID
        output.append(strains[0])                 # SimHrss
        output.append(energy)                     # SimEgwR2
        output.append(strains[0])                 # GravEn_Ampl
        output.append(cosincl)           # Internal_x the cosine of the angle the LOS makes with axis of angular momentum
        output.append(row.phi)                    # Intenal_phi angle between source x-axis and the LOS
        output.append(np.cos(np.pi/2.0 - row.dec)) # cos(External_x) # this needs to be the co-declination
        output.append(row.ra if row.ra < np.pi else row.ra - 2*np.pi)
        # ^ External_phi # This is the RA projected onto an Earth-based coordinate system
        output.append(row.psi)                    # External_psi # source's polarisation angle
        output.append(frame.start)                # FrameGPS
        output.append(row.time_geocent_gps)       # EarthCtrGPS
        output.append(rowname)                    # SimName
        output.append(strains[1])                 # SimHpHp
        output.append(strains[2])                 # SimHcHc
        output.append(strains[3])                 # SimHpHp
        output.append(" ".join(" ".join(map(str,l)) for l in responses))
        return ' '.join(str(e) for e in output)

class Frame():
    """
    Represents a frame, in order to prepare the injection frames
    """
    def __init__(self, start, duration, ifo, number = -1):
        """

        Parameters
        ----------
        number : int
           The frame's number within the project. Defaults to -1.
        """
        self.start = start
        self.duration = duration
        self.end = self.start + duration
        self.ifos = ifo
        self.number = -1
        
    def __repr__(self):
        out = ''
        out += "MDC Frame \n"
        for ifo in self.ifos:
            out += "{} {} {} \n".format(ifo, self.start, self.duration)
        return out
    
    def get_rowlist(self,mdcs):
        """
        Return the rows from an MDC set which correspond to this frame.
        
        Parameters
        ----------
        mdcs : MDCSet object
            The set of MDCs from which the rows are to be found.
        """
        return np.where((mdcs.times<self.end)&(mdcs.times>self.start))[0]
    
    def calculate_n_injections(self, mdcs):
        return len(mdcs.times[(mdcs.times<self.end)&(mdcs.times>self.start)])
    
    def generate_log(self,mdc):
        log = '#  GravEn_SimID  SimHrss  SimEgwR2  GravEn_Ampl  Internal_x  Internal_phi  External_x  External_phi External_psi  FrameGPS  EarthCtrGPS  SimName  SimHpHp  SimHcHc  SimHpHc  H1       H1ctrGPS        H1fPlus        H1fCross    L1       L1ctrGPS        L1fPlus        L1fCross\n'
        rowlist = self.get_rowlist(mdc)
        for row in rowlist:
            log += mdc.gravEn_row(row, self)
            log += "\n"
        return log

    def generate_gwf(self, mdc, directory, project = "Minke", channel="SCIENCE", force=False, rate=16384.0):
        """
        Produce the gwf file which corresponds to the MDC set over the period of this frame.

        Parameters
        ----------
        mdc : MDCSet object
           The MDC set which should be used to produce this frame.
        directory : str
           The root directory where all of the frames are to be stored, for example
           "/home/albert.einstein/data/mdc/frames/"
           would cause the SineGaussian injections to be made in the directories under
           "/home/albert.einstein/data/mdc/frames/sg"
        project : str
           The name of the project which this frame is a part of. Defaults to 'Minke'.
        channel : str
           The name of the channel which the injections should be made into. This is prepended by the initials
           for each interferometer, so there will be a channel for each interferometer in the gwf.
        force : bool
           If true this forces the recreation of a GWF file even if it already exists.

        Outputs
        -------
        gwf
           The GWF file for this frame.
        """
        ifosstr = "".join(set(ifo[0] for ifo in self.ifos))
        family = mdc.waveforms[0].waveform
        epoch = lal.LIGOTimeGPS(self.start)
        filename = "{}-{}-{}-{}.gwf".format(ifosstr, family, self.start, self.duration)

        self.frame = lalframe.FrameNew(epoch = epoch,
                                       duration = self.duration, project='', run=1, frnum=1,
                                       detectorFlags=lal.LALDETECTORTYPE_ABSENT)


        ifobits = np.array([getattr(lal,"{}_DETECTOR_BIT".format(lal.cached_detector_by_prefix[ifo].frDetector.name.upper()))
                   for ifo in self.ifos])
        ifoflag = numpy.bitwise_or.reduce(ifobits)
        RUN_NUM = -1 # Simulated data should have a negative run number

        head_date = str(self.start)[:5]
        frameloc = directory+"/"+mdc.directory_path()+"/"+head_date+"/"
        mkdir(frameloc)
        if not os.path.isfile(frameloc + filename) or force:
            epoch = lal.LIGOTimeGPS(self.start)
            frame = lalframe.FrameNew(epoch, self.duration, project, RUN_NUM, self.number, ifoflag)
            data = []
            # Loop through each interferometer
            
            for ifo in self.ifos:
                # Calculate the number of samples in the timeseries
                nsamp = int((self.end-self.start)*rate)
                # Make the timeseries
                h_resp = lal.CreateREAL8TimeSeries("{}:{}".format(ifo, channel), epoch, 0, 1.0/rate, lal.StrainUnit, nsamp)
                # Loop over all of the injections corresponding to this frame
                rowlist = self.get_rowlist(mdc)
                
                if len(rowlist)==0: return
                for row in rowlist:
                    sim_burst = mdc.waveforms[row]._row()

                    if sim_burst.hrss > 1:
                        distance = sim_burst.amplitude
                    else:
                        distance = None
                    
                    #hp, hx = lalburst.GenerateSimBurst(sim_burst, 1.0/rate);
                    hp, hx, _, _ = mdc.waveforms[row]._generate(rate=rate, half=True, distance=distance)
                    # Apply detector response
                    det = lalsimulation.DetectorPrefixToLALDetector(ifo)
                    # Produce the total strains
                    
                    h_tot = lalsimulation.SimDetectorStrainREAL8TimeSeries(hp, hx,
                                                                           sim_burst.ra, sim_burst.dec, sim_burst.psi, det)
                    # Inject the waveform into the overall timeseries
                    lalsimulation.SimAddInjectionREAL8TimeSeries(h_resp, h_tot, None)

                lalframe.FrameAddREAL8TimeSeriesSimData(frame, h_resp)

            # Make the directory in which to store the files
            # if it doesn't exist already

            # Write out the frame file
            lalframe.FrameWrite(frame, frameloc+filename)
        

class HWInj(Frame):
    """
    Represents a hardware injection frame.
    
    Injection frames must be an ASCII file of the hoft sampled at 
    the antenna sampling rate, appropriately convolved with an 
    antenna response function.

    As a result of the simplicity of this specific output format
    we do not need information such as start-time in the file itself,
    however we should have a sensible naming scheme for the ASCII files
    since they will need to be produced as sidecars for an xml file.

    
    """
    def __init__(self, ifos):
        """We'll need to know the start-time, the duration, and the ifo
        for each which is to be used for hardware injections in order
        to keep consistency with the data in the xml file, and so that the 
        appropriate waveform is injected into the appropriate detector.

        Parameters
        ----------
        ifos : list
           The name of the interferometers, e.g. "L1" for the Livingston, LA LIGO detector.

        """
        self.ifos = ifos

    def __repr__(self):
        """
        The printable representation of this object.
        """
        out = ""
        out += "Hardware MDC Frame \n"
        for ifo in self.ifos:
            out += "{} \n".format(ifo)
        return out

    def generate_pcal(self, mdc, directory, force = False, rate=16384):
        """
        Produce the PCAL-ready hardware injection files as an ASCII list
        sampled at the detector's sample rate.

        Parameters
        ----------
        mdc : MDCSet object
           The signal set which should be used to generate the frame.
        directory : str
           The root directory where all of the frames are to be stored, for example
           "/home/albert.einstein/data/mdc/frames/"
           would cause the SineGaussian injections to be made in the directories under
           "/home/albert.einstein/data/mdc/frames/sg"
        force : bool
           If true this forces the regeneration of the file, even if it
           already exists.
        
        Outputs
        -------
        ascii file
           The ASCII file containing the correctly sampled waveform convolved with
           the antenna pattern.
        """
        
        family = mdc.waveforms[0].waveform

        frameloc = os.path.join(directory, (mdc.directory_path()))
        #rowlist = self.get_rowlist(mdc)
        # Unlike with a conventional frame, we need to produce a separate file
        # for each IFO.
        for ifo in self.ifos:
            for sim_burst in mdc.waveforms:
                #sim_burst = mdc.waveforms[row]
                # Check if the file exists, or if we're forcing the creation
                filename = "{}_{}_{}.txt".format(family, 
                                                 sim_burst.time, 
                                                 ifo)
                if not os.path.isfile(frameloc + filename) or force:
                    data = []
                    epoch = lal.LIGOTimeGPS(sim_burst.time)
                    duration = 10
                    nsamp = duration*rate

                    h_tot = sim_burst._generate_for_detector([ifo], sample_rate=rate)
                    
                    data = np.array(h_tot.data.data)
                    np.savetxt(filename, data)

class HWFrameSet():
    def __init__(self, ifos=["H1", "L1"]):
        """
        A collection of hardware injection frames.

        Parameters
        ----------
        frame_list : str
            The filespath of a CSV file containing the list of frames, 
            and the parameters required to produce them: the start and 
            duration times, and the interferometers they describe.
        """
    
        self.frames = []
        self.frames = [HWInj(ifos)]
        #self.frames.append(frame)

    def full_frameset(self, mdc, directory, force=False):
        """
        Produce the gwf files which corresponds to the MDC set over the period of the frames in this collection.

        Parameters
        ----------
        mdc : MDCSet object
           The MDC set which should be used to produce this frame.
        directory : str
           The root directory where all of the frames are to be stored, for example
           "/home/albert.einstein/data/mdc/frames/"
           would cause the SineGaussian injections to be made in the directories under
           "/home/albert.einstein/data/mdc/frames/sg"
        force : bool
           If true this forces the recreation of a GWF file even if it already exists.

        Outputs
        -------
        ascii files
           The ASCII files for these hardware injections.
        """
        for frame in self.frames:
            frame.generate_pcal(mdc, directory, force)

class FrameSet():

    def __init__(self, frame_list):
        """
        A collection of frames.

        Parameters
        ----------
        frame_list : str
            The filespath of a CSV file containing the list of frames, 
            and the parameters required to produce them: the start and 
            duration times, and the interferometers they describe.
        """

        self.frames = []
        self.frame_list = frame_list = pd.read_csv(frame_list)
        for frame in frame_list.iterrows():
            frame = frame[1]
            ifos = frame['ifo'].replace("['",'').replace("']",'').replace("'",'').split(' ')
            frame = Frame(frame['start time'],frame['duration'],ifos)
            self.frames.append(frame)
        
    def full_frameset(self, mdc, directory, channel="SCIENCE", force=False):
        """
        Produce the gwf files which corresponds to the MDC set over the period of the frames in this collection.

        Parameters
        ----------
        mdc : MDCSet object
           The MDC set which should be used to produce this frame.
        directory : str
           The root directory where all of the frames are to be stored, for example
           "/home/albert.einstein/data/mdc/frames/"
           would cause the SineGaussian injections to be made in the directories under
           "/home/albert.einstein/data/mdc/frames/sg"
        channel : str
           The name of the channel which the injections should be made into. This is prepended by the initials
           for each interferometer, so there will be a channel for each interferometer in the gwf.
        force : bool
           If true this forces the recreation of a GWF file even if it already exists.

        Outputs
        -------
        gwf files
           The GWF files for these frames.
        """
        for frame in self.frames:
            frame.generate_gwf(mdc, directory, channel, force)


    def full_logfile(self, mdc, location):
        """
        Produce a log file for the entire frame set 
        """
        full_log = ''
        for frame in self.frames:
            full_log += frame.generate_log(mdc)
            
        with open(location, "w") as text_file:
            text_file.write(full_log)

