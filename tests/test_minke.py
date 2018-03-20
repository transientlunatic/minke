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

class TestMinke(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    # Check that failure occurs if the wrong waveforms are inserted into the wrong tables

    def test_insert_wrong_waveform_to_table(self):
        """Test whether adding a ringdown waveform to a SimBurstTable throws
        an error.
        """
        mdcset = mdctools.MDCSet(["L1"], table_type = "burst")
        ring = sources.Ringdown()

        with self.assertRaises(mdctools.TableTypeError):
            mdcset + ring

    def test_insert_burst_waveform_to_burst_table(self):
        """
        Test whether inserting the correct type of waveform into a table works
        """
        mdcset = mdctools.MDCSet(["L1"], table_type = "burst")
        waveform = sources.Gaussian(0,0,0)

        mdcset + waveform

        self.assertEqual(len(mdcset.waveforms), 1)

    # Check that the correct XML files are produced for different table types

    def test_write_simbursttable(self):
        """
        Write out a simburst table xml file
        """
        mdcset = mdctools.MDCSet(["L1"], table_type = "burst")
        waveform = sources.Gaussian(0.1,1e-23,1000)

        mdcset + waveform

        mdcset.save_xml("test_simbursttable.xml.gz")

    def test_verify_simbursttable(self):
        """
        Read-in the xml simburst table.
        """
        mdcset = mdctools.MDCSet(["L1"])
        mdcset.load_xml("test_simbursttable.xml.gz", full = False)

        self.assertEqual(len(mdcset.waveforms), 1)

    def test_write_simringdowntable(self):
        """
        Write out a simburst table xml file
        """
        mdcset = mdctools.MDCSet(["L1"], table_type = "ringdown")
        waveform = sources.BBHRingdown(1000, 1e-22, 0 ,0.1, 10, 0.1, 0.01, 10, 0, 2, 2,)

        mdcset + waveform

        mdcset.save_xml("test_simringdowntable.xml.gz")

    def test_verify_simringdowntable(self):
        """
        Read-in the xml simburst table.
        """
        mdcset = mdctools.MDCSet(["L1"], table_type = "ringdown",)
        mdcset.load_xml("test_simringdowntable.xml.gz", full = False)

        self.assertEqual(len(mdcset.waveforms), 1)

        
        
if __name__ == '__main__':
    import sys
    sys.exit(unittest.main())
