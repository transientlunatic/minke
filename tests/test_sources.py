import unittest
from minke import sources
from minke import mdctools
from minke import distribution
import numpy as np

from igwn_ligolw import lsctables

ott_data_url = "http://www.stellarcollapse.org/gwdata/ottetal2012b/s27WHW02_LS220_j0_rx3_full_cc_fheat1.05_HlmD.dat.gz"

def check_no_nr_support():
    """
    Check to see if this installation does not come with the LALSuite support for NR waveforms.
    """
    if "numrel_data" in lsctables.SimBurstTable.validcolumns.keys():
        NROK = True
    else:
        NROK = False

    return ~NROK

NO_NR_REASON = "This test was skipped as the installation has no support for NR waveforms."


@unittest.skipIf(check_no_nr_support(), NO_NR_REASON)
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
            print("Getting {}".format(url))
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

    return fname.strip(".gz")

class TestMinkeAdHocSources(unittest.TestCase):
    """
    Tests for the adhoc analytical waveforms
    """

    def setUp(self):
        """
        Set everything up to make things work.
        """
        self.mdcset = mdctools.MDCSet(['L1', 'H1'])
        self.times = distribution.even_time(start = 1126620016, stop = 1136995216, rate = 630720, jitter = 20)
        self.angles = distribution.supernova_angle(len(self.times))

    def test_Gaussian_Waveform_generation(self):
        """Test that Gaussian waveforms are generated sensibly"""
        ga = sources.Gaussian(1,1,1)
        data = ga._generate()

        gadata = np.array([  3.22991378e-57,   1.68441642e-25,   1.43167033e-23,
         6.20128276e-22,   1.92272388e-20,   4.74643578e-19,
         9.78126889e-18,   1.72560021e-16,   2.64547116e-15,
         3.55831241e-14,   4.22637998e-13,   4.45285513e-12,
         4.17503728e-11,   3.49184250e-10,   2.60952919e-09,
         1.74465929e-08,   1.04437791e-07,   5.60038118e-07,
         2.70506254e-06,   1.18998283e-05,   4.76932871e-05,
         1.74151408e-04,   5.79361856e-04,   1.75600595e-03,
         4.84903412e-03,   1.21993775e-02,   2.79623250e-02,
         5.83931709e-02,   1.11097431e-01,   1.92574662e-01,
         3.04121728e-01,   4.37571401e-01,   5.73592657e-01,
         6.85032870e-01,   7.45370866e-01,   7.38901572e-01,
         6.67350426e-01,   5.49129112e-01,   4.11669004e-01,
         2.81173926e-01,   1.74966576e-01,   9.91946657e-02,
         5.12359410e-02,   2.41109482e-02,   1.03372980e-02,
         4.03787583e-03,   1.43698441e-03,   4.65912462e-04,
         1.37628942e-04,   3.70397795e-05,   9.08197260e-06,
         2.02882766e-06,   4.12393041e-07,   7.53383218e-08,
         1.23283052e-08,   1.80606751e-09,   2.36657746e-10,
         2.77012812e-11,   2.89124066e-12,   2.68404945e-13,
         2.20863500e-14,   1.60321302e-15,   1.01948770e-16,
         5.62074344e-18,   2.64306001e-19,   1.03069500e-20,
         3.15739881e-22,   6.68755767e-24,   6.21121403e-26])

        np.testing.assert_array_almost_equal(data[0].data.data[::5000], gadata)
        
    def test_SG_Waveform_generation(self):
        """
        Regression test for SineGaussian Waveforms.
        """
        
        sg = sources.SineGaussian(10,1,1,"linear",1, seed = 0)
        data = sg._generate()

        sgdata = np.array([ -8.73180626e-58,  -1.93967534e-26,   5.31448156e-25,
         2.85625032e-24,  -1.16897801e-22,   4.54221716e-22,
         7.30222266e-21,  -7.05284947e-20,  -7.18443386e-20,
         3.83540555e-18,  -1.35146997e-17,  -9.19982707e-17,
         8.18155341e-16,  -1.75844610e-16,  -2.20122350e-14,
         7.73708586e-14,   2.54778729e-13,  -2.36261229e-12,
         2.01404378e-12,   3.49131470e-11,  -1.23199563e-10,
        -1.89419681e-10,   2.06578175e-09,  -2.48373900e-09,
        -1.71767992e-08,   6.11421635e-08,   3.47780769e-08,
        -5.83248940e-07,   8.10839193e-07,   2.72782093e-06,
        -1.00820213e-05,  -1.30630395e-07,   5.80605496e-05,
        -8.87448223e-05,  -1.55818895e-04,   6.31437341e-04,
        -2.04866297e-04,  -2.23611311e-03,   3.62241631e-03,
         3.19299434e-03,  -1.53446556e-02,   8.16204089e-03,
         3.31744792e-02,  -5.62402740e-02,  -2.09500348e-02,
         1.44826566e-01,  -9.53929854e-02,  -1.87967432e-01,
         3.35588750e-01,   2.19275523e-02,  -5.30630637e-01,
         3.90608881e-01,   4.00367320e-01,  -7.74119199e-01,
         9.53218059e-02,   7.53180228e-01,  -5.91727767e-01,
        -3.10327911e-01,   6.92512672e-01,  -1.65601830e-01,
        -4.12458030e-01,   3.39597507e-01,   8.05475578e-02,
        -2.40582205e-01,   7.46124228e-02,   8.65034217e-02,
        -7.47429394e-02,  -4.86389437e-03,   3.24540997e-02,
        -1.14700762e-02,  -6.85576126e-03,   6.35159480e-03,
        -2.55523511e-04,  -1.69725819e-03,   6.45883694e-04,
         1.99966867e-04,  -2.09187908e-04,   2.19763186e-05,
         3.42911180e-05,  -1.37134844e-05,  -2.01577460e-06,
         2.67493920e-06,  -3.87356579e-07,  -2.65313667e-07,
         1.10113354e-07,   5.46563207e-09,  -1.26429809e-08,
         2.07544399e-09,   6.97707443e-10,  -2.97465468e-10,
         2.38449797e-12,   1.95277082e-11,  -3.35619371e-12,
        -5.64458960e-13,   2.55522574e-13,  -1.03406456e-14,
        -9.34634494e-15,   1.60335007e-15,   1.24686211e-16,
        -6.41126505e-17,   3.59252308e-18,   1.20462046e-18,
        -1.94100550e-19,  -5.20175883e-21,   3.52788114e-21,
        -1.92473065e-22,  -2.39883058e-23,   2.65759767e-24,
         2.11709555e-27,  -2.28657022e-27])

        np.testing.assert_array_almost_equal(data[0].data.data[::5000], sgdata)

    def test_SC_Waveform_generation(self):
        """
        Regression test for StringCusp waveforms
        """
        
        sc = sources.StringCusp(1e-23, 150, 100)
        data = sc._generate()

        scdata = np.array([  5.70001763e-63,   5.92586457e-32,   2.59625621e-31,
         1.29668495e-30,   4.86418685e-30,   1.25517562e-29,
         1.95684680e-29,  -2.79296626e-29,  -3.97865144e-28,
        -1.98228555e-27,  -5.83672705e-27,  -1.84365082e-27,
         7.98894273e-26,   2.10292180e-25,  -2.26513533e-24,
        -1.58525748e-24,  -2.79263187e-25,   1.76952837e-25,
         2.03163416e-26,  -6.90900890e-27,  -3.64716334e-27,
        -9.34585595e-28,  -1.33656482e-28,   1.26838192e-29,
         1.67225169e-29,   7.30560692e-30,   2.67365656e-30,
         8.41779831e-31,   1.43679337e-31,   6.07001624e-33])

        np.testing.assert_array_almost_equal(data[0].data.data[::5000], scdata)


class TestMinkeWNBSources(unittest.TestCase):
    """
    Run tests on the creation of white noise bursts to 
    ensure that these are produced as expected.
    """

    def setUp(self):
        """
        Set everything up to make things work.
        """
        self.mdcset = mdctools.MDCSet(['L1', 'H1'])
        self.times = distribution.even_time(start = 1126620016, stop = 1136995216, rate = 630720, jitter = 20)

    def test_random_seed(self):
        """Check that the WNB source produces the same waveform each time it is
        called with the same seed."""

        sc = sources.WhiteNoiseBurst(0.1, 10, 1024, 1126620016, hrss=1e-22)

        hp, hx, _, _ = sc._generate()
        hp2, hx2, _, _ = sc._generate()

        np.testing.assert_array_almost_equal(hp.data.data, hp2.data.data)

    
        
class TestMinkeADISources(unittest.TestCase):
    def setUp(self):
        """
        Set things up for the tests by making the MDC set, and defining the various parameter distributions.
        """

        self.datafiles = {}
        self.datafiles['ADI'] = 'tests/data/adi_test.mat'

        mdctools.mkdir("./testout/frames")

        self.mdcset = mdctools.MDCSet(['L1', 'H1'])
        self.times = distribution.even_time(start = 1126620016, stop = 1136995216, rate = 630720, jitter = 20)
        self.angles = distribution.supernova_angle(len(self.times))

    def test_ADIWaveform(self):
        """Check the generation of the ADI waveform."""

        adi  = sources.ADI(time = 100, filepath=self.datafiles['ADI'])
        adiwave = adi._generate()[0].data.data[::20000]
        reg = np.array([  0.00000000e+00,   6.62293373e-22,   1.52831641e-21,
         1.54044027e-21,   1.27707122e-21,  -2.19617692e-22,
        -1.52414321e-21,   4.69396022e-22,   1.21479442e-21,
        -1.26396697e-21,  -5.87290821e-22,   1.36010407e-21,
         1.18247173e-21,   3.76119821e-22,   3.08731444e-22,
         1.14768993e-21,   1.08488237e-21,  -1.50088643e-21,
         1.45463205e-21,   1.67209738e-23,  -8.44421117e-22,
         9.54960410e-23,   1.29069231e-21,  -9.77430178e-22,
        -1.34761100e-21,  -1.47183099e-21,   4.61589317e-22,
        -1.18249701e-21,  -3.27337305e-22,  -1.45357383e-21,
         1.45374218e-21,   4.60667061e-22])

        np.testing.assert_array_almost_equal(adiwave, reg)

class TestMinkeSupernovaSources(unittest.TestCase):
    def setUp(self):
        """
        Set things up for the tests by making the MDC set, and defining the various parameter distributions.
        """

        self.datafiles = {}
        self.datafiles['Ott13'] = 'tests/data/ott_test.dat'
        self.datafiles['Yakunin10'] = "tests/data/yakunin_test.dat"

        mdctools.mkdir("./testout/frames")

        self.mdcset = mdctools.MDCSet(['L1', 'H1'])
        self.times = distribution.even_time(start = 1126620016, stop = 1136995216, rate = 630720, jitter = 20)
        self.angles = distribution.supernova_angle(len(self.times))

    #### TEST YAKUNIN
    @unittest.skipIf(check_no_nr_support(), NO_NR_REASON)
    def test_YakuninWaveform(self):
        """Test generation of a Yakunin waveform."""

        reg = np.array([
            -8.15750499e-23,  -5.99101709e-23,  -8.85420897e-23,
            -7.56873580e-23,  -4.13010667e-23,  -2.76929632e-23,
            -6.17073634e-23,  -7.84874142e-23,   1.27751700e-23,
            2.29718421e-23,   4.42735947e-23,   5.19635266e-23,
            4.74454890e-23,   4.91091850e-23,   4.82164382e-23,
            4.64586078e-23,   4.59085183e-23,   4.48057968e-23,
            4.30575525e-23,   4.15572222e-23,   4.61748203e-23,
            4.09097191e-23,   4.09360107e-23,   3.93016465e-23,
            3.57392283e-23,   3.03996763e-23,   2.34990480e-23,
            1.53106213e-23,   6.15487802e-24,  -3.61224433e-24,
            -1.36123962e-23,  -2.34586878e-23,  -3.27698685e-23,
            -4.11843030e-23,  -4.83732254e-23,  -5.40527874e-23,
            -5.79944705e-23,  -6.00334561e-23
        ])
        
        sn = sources.Yakunin10(time=1000, 
                               filepath = self.datafiles['Yakunin10'],
        )
        np.testing.assert_array_almost_equal(sn._generate()[0].data.data[::1000], reg)

    @unittest.skipIf(check_no_nr_support(), NO_NR_REASON)
    def test_YakuninWaveform_distance(self):
        sn_10 = sources.Yakunin10(time=1126620016,
                             filepath=self.datafiles['Yakunin10'],
        )
        sn_100 = sources.Yakunin10(time=1126620016,
                             filepath=self.datafiles['Yakunin10'],
                             distance = 100e-3,
        )
        self.assertAlmostEqual(sn_10.params['amplitude'],float(0.1 * sn_100.params['amplitude']))
        sn_10_hp, _, _, _ = sn_10._generate()
        sn_100_hp, _, _, _ = sn_100._generate()
        
        np.testing.assert_almost_equal(10*sn_100_hp.data.data, sn_10_hp.data.data)                           
        
    #### TEST OTT
    @unittest.skipIf(check_no_nr_support(), NO_NR_REASON)
    def test_OttWaveform(self):
        sn = sources.Ott2013(theta = 0 , phi = 0, time=1126620016,
                             filepath=self.datafiles['Ott13'],
                             family = "s27fheat1p05",
        )

    
    @unittest.skipIf(check_no_nr_support(), NO_NR_REASON)
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
    
    @unittest.skipIf(check_no_nr_support(), NO_NR_REASON)    
    def test_OttXML(self):
        mdcset = mdctools.MDCSet(['L1', 'H1'])
        sn = sources.Ott2013(theta = 0 , phi = 0, time=1126620016,
                             filepath=self.datafiles['Ott13'],
                             family = "s27fheat1p05",)
        mdcset + sn

        mdcset.save_xml('testout/ott13.xml.gz')

    @unittest.skipIf(check_no_nr_support(), NO_NR_REASON)
    def test_OttFrame(self):
        mdcset = mdctools.MDCSet(['L1', 'H1'])
        mdcset.load_xml('tests/data/ott_test.xml.gz')
        o1 = mdctools.FrameSet('tests/data/frame_list.dat')
        print(o1.frames[0].generate_gwf(mdcset, "/home/daniel"))
        mdc_folder = "testout/frames"
        for o1frame in o1.frames:
            o1frame.generate_gwf(mdcset, mdc_folder, 'SCIENCE')

    @unittest.skipIf(check_no_nr_support(), NO_NR_REASON)
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
