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
        self.datafiles['Ott13'] = 'tests/data/ott_test.dat' #download_nr(ott_data_url)

        mdctools.mkdir("./testout/frames")

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
        
        np.testing.assert_almost_equal(10*sn_100_hp.data.data, sn_10_hp.data.data)

    def test_OttXML(self):
        mdcset = mdctools.MDCSet(['L1', 'H1'])
        sn = sources.Ott2013(theta = 0 , phi = 0, time=1126620016,
                             filepath=self.datafiles['Ott13'],
                             family = "s27fheat1p05",)
        mdcset + sn

        mdcset.save_xml('testout/ott13.xml.gz')

    def test_OttFrame(self):
        mdcset = mdctools.MDCSet(['L1', 'H1'])
        mdcset.load_xml('tests/data/ott_test.xml.gz')
        o1 = mdctools.FrameSet('tests/data/frame_list.dat')
        mdc_folder = "testout/frames"
        for o1frame in o1.frames:
            o1frame.generate_gwf(mdcset, mdc_folder, 'SCIENCE')

    def test_OttGraven(self):
        mdcset = mdctools.MDCSet(['L1', 'H1'])
        mdcset.load_xml('tests/data/ott_test.xml.gz')
        o1 = mdctools.FrameSet('tests/data/frame_list.dat')
        o1.full_logfile(mdcset, 'testout/frames/logfile.txt')

    def TestRingdown(self):

        testdata = np.array(
            [ -4.14099902e-20,   4.05706892e-21,   9.44712664e-22,
              -2.70287547e-22,   1.98245659e-24,   9.41028314e-24,
              -1.31703741e-24,  -1.62366854e-25,   6.86331926e-26,
              -3.27763986e-27,  -2.02216144e-27,   3.85067715e-28,
              2.13773021e-29,  -1.66111808e-29,   1.43458903e-30,
              4.04499843e-31,  -1.04902792e-31,  -5.84838741e-34,
              3.83160548e-33,  -4.86443978e-34,  -7.27079363e-35,
              2.70362403e-35]
        )
        
        mdcset = mdctools.MDCSet(["L1"], table_type = "ringdown")
        waveform = sources.BBHRingdown(1000, 1e-22, 0 ,0.1, 10, 0.1, 0.01, 10, 0, 2, 2,)

        mdcset + waveform

        hp, hx = waveform._generate()

        np.testing.assert_array_equal(hp, testdata)
        
    def TestSupernovaXML(self):
        
        pass

if __name__ == '__main__':
    import sys
    sys.exit(unittest.main())
