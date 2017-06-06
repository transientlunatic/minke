import unittest
from minke import sources
from minke import mdctools
from minke import distribution
import numpy as np


ott_data_url = "http://www.stellarcollapse.org/gwdata/ottetal2012b/s27WHW02_LS220_j0_rx3_full_cc_fheat1.05_HlmD.dat.gz"

def download_nr(url):
    """
    Download the datafile for an NR waveform, and extract it from a gzip if relevent
    """
    import urllib2, gzip, StringIO

    attempts = 0

    while attempts < 3:
        try:
            response = urllib2.urlopen(url, timeout = 5)
            fname = url.split('/')[-1]
            print "Getting {}".format(url)
            content = response.read()
            with open( fname, 'w' ) as f:
                f.write( content )
            with open (fname, 'rb') as f:
                decomp = gzip.GzipFile(fileobj = f, mode='rb')
                with open(fname.strip(".gz"), 'w') as outfile:
                    outfile.write(decomp.read())
            break
        except urllib2.URLError as e:
            attempts += 1
            print type(e)

    return fname.strip(".gz")

class TestMinkeSources(unittest.TestCase):
    def setUp(self):
        """
        Set things up for the tests by making the MDC set, and defining the various parameter distributions.
        """

        self.datafiles = {}
        self.datafiles['Ott13'] = download_nr(ott_data_url)

        self.mdcset = mdctools.MDCSet(['L1', 'H1'])
        self.times = distribution.even_time(start = 1126620016, stop = 1136995216, rate = 630720, jitter = 20)
        self.angles = distribution.supernova_angle(len(self.times))

    #### TEST OTT

    def test_OttWaveform(self):
        sn = sources.Ott2013(theta = 0 , phi = 0, time=1126620016,
                             filepath=self.datafiles['Ott13'],
                             family = "s27fheat1p05",
        )

    def test_OttWaveform_distance(self):
        sn_10 = sources.Ott2013(theta = 0 , phi = 0, time=1126620016,
                             filepath=self.datafiles['Ott13'],
                             family = "s27fheat1p05",
        )
        sn_100 = sources.Ott2013(theta = 0 , phi = 0, time=1126620016,
                             filepath=self.datafiles['Ott13'],
                             family = "s27fheat1p05",
                             distance = 100e-3,
        )
        self.assertAlmostEqual(sn_10.params['amplitude'],float(0.1 * sn_100.params['amplitude']))
        sn_10_hp, _, _, _ = sn_10._generate()
        sn_100_hp, _, _, _ = sn_100._generate()
        assert np.all(sn_100_hp.data.data == 10*sn_10_hp.data.data)

    def test_OttXML(self):
        mdcset = mdctools.MDCSet(['L1', 'H1'])
        sn = sources.Ott2013(theta = 0 , phi = 0, time=1126620016,
                             filepath=self.datafiles['Ott13'],
                             family = "s27fheat1p05",)
        mdcset + sn

        mdcset.save_xml('ott13.xml.gz')

    def test_OttFrame(self):
        mdcset = mdctools.MDCSet(['L1', 'H1'])
        mdcset.load_xml('ott13.xml.gz')
        o1 = mdctools.FrameSet('tests/data/frame_list.dat')
        mdc_folder = "./frames"
        for o1frame in o1.frames:
            o1frame.generate_gwf(mdcset, mdc_folder, 'SCIENCE')

    def test_OttGraven(self):
        mdcset = mdctools.MDCSet(['L1', 'H1'])
        mdcset.load_xml('ott13.xml.gz')
        o1 = mdctools.FrameSet('tests/data/frame_list.dat')
        o1.full_logfile(mdcset, './frames/logfile.txt')


        
    def TestSupernovaXML(self):
        
        pass

if __name__ == '__main__':
    import sys
    sys.exit(unittest.main())
