#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
test_minke
----------------------------------

Tests for `minke` module.
"""

import unittest

import minke
from minke import mdctools, sources
import numpy as np

class TestMinke(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    # Check that failure occurs if the wrong waveforms are inserted into the wrong tables

    def test_insert_wrong_waveform_to_table(self):
        """Test whether adding a ringdown waveform to a SimBurstTable throws an error."""
        mdcset = mdctools.MDCSet(["L1"], table_type = "burst")
        ring = sources.Ringdown()

        with self.assertRaises(mdctools.TableTypeError):
            mdcset + ring

    def test_insert_burst_waveform_to_burst_table(self):
        """Test whether inserting the correct type of waveform into a table works."""
        mdcset = mdctools.MDCSet(["L1"], table_type = "burst")
        waveform = sources.Gaussian(0,0,0)

        mdcset + waveform

        self.assertEqual(len(mdcset.waveforms), 1)

    # Check that the correct XML files are produced for different table types

    def test_readwrite_simbursttable(self):
        """Test writing out a simburst table xml file."""
        mdcset = mdctools.MDCSet(["L1"], table_type = "burst")
        waveform = sources.Gaussian(0.1,1e-23,1000)

        mdcset + waveform

        mdcset.save_xml("testout/test_simbursttable.xml.gz")

        mdcset = mdctools.MDCSet(["L1"])
        mdcset.load_xml("testout/test_simbursttable.xml.gz", full = False)

        self.assertEqual(len(mdcset.waveforms), 1)

    def test_readwrite_simringdowntable(self):
        """Write out a simringdown table xml file
        """
        mdcset = mdctools.MDCSet(["L1"], table_type = "ringdown")
        waveform = sources.BBHRingdown(1000,  0 ,0.1, 10, 0.1, 0.01, 10, 0)

        mdcset + waveform

        mdcset.save_xml("testout/test_simringdowntable.xml.gz")

        mdcset = mdctools.MDCSet(["L1"], table_type = "ringdown",)
        mdcset.load_xml("testout/test_simringdowntable.xml.gz", full = False)

        self.assertEqual(len(mdcset.waveforms), 1)


class TestMDC(unittest.TestCase):
    def setUp(self):
        mdctools.mkdir("./testout")
        
        rng = minke.sources.BBHRingdown(100.23, np.rad2deg(.1), 10., 0.97, 0.01, 1.0,45.)
        
        mdcset = mdctools.MDCSet(["L1"], table_type="ringdown")
        mdcset + rng
        mdcset.save_xml("testout/ringdowngentest.xml.gz")
    
    def test_Gaussian_Waveform_generation_in_MDC(self):
        """Check that waveforms inside an MDC are generated. """
        ga = sources.Gaussian(1,1,1)

        mdcset = mdctools.MDCSet(["L1"])
        mdcset + ga
        
        data = mdcset._generate_burst(0)

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


    def test_Ringdown_Waveform_generation_in_MDC(self):
        """Check that waveforms inside an MDC are generated. """       
        mdcset = mdctools.MDCSet(["L1"], table_type="ringdown")
        mdcset.load_xml("testout/ringdowngentest.xml.gz")
        
        data = mdcset._generate_burst(0)

        testdata = np.array([
            -4.14099902e-20,   4.05706892e-21,   9.44712664e-22,
            -2.70287547e-22,   1.98245659e-24,   9.41028314e-24,
            -1.31703741e-24,  -1.62366854e-25,   6.86331926e-26,
            -3.27763986e-27,  -2.02216144e-27,   3.85067715e-28,
            2.13773021e-29,  -1.66111808e-29,   1.43458903e-30,
            4.04499843e-31,  -1.04902792e-31,  -5.84838741e-34,
            3.83160548e-33,  -4.86443978e-34,  -7.27079363e-35,
            2.70362403e-35])

        np.testing.assert_array_almost_equal(data[0].data.data[::30], testdata)

class TestMDCEdit(unittest.TestCase):
    def setUp(self):
        ga = sources.Gaussian(1,1,1)

        mdcset = mdctools.MDCSet(["L1"])
        mdcset + ga
        
        mdctools.mkdir("./testout")
        mdcset.save_xml("testout/gaussian_edit_test.xml")

    def testXMLLoad(self):
        """Check that XMLs defining an MDC can be loaded."""
        mdcset = mdctools.MDCSet(["L1"])
        mdcset.load_xml("testout/gaussian_edit_test.xml")

        self.assertEqual(len(mdcset.waveforms), 1)

        self.assertEqual(type(mdcset.waveforms[0]), minke.sources.Gaussian)
        
    def testXMLAddition(self):
        """Check that rows can be added to an XML."""
        mdcset = mdctools.MDCSet(["L1"])
        mdcset.load_xml("testout/gaussian_edit_test.xml")
        ga = sources.Gaussian(2,2,2)
        mdcset + ga
        self.assertEqual(len(mdcset.waveforms), 2)

    def testXMLRowEdit(self):
        """Check that you can edit the contents of an XML row."""
        mdcset = mdctools.MDCSet(["L1"])
        mdcset.load_xml("testout/gaussian_edit_test.xml")
        mdcset.waveforms[0].amplitude = 3

        self.assertEqual(mdcset.waveforms[0].amplitude, 3)
        
if __name__ == '__main__':
    import sys
    sys.exit(unittest.main())
