import unittest

import astropy.units as u
import torch
import numpy as np
import gpytorch

from heron.training.makedata import make_manifold, make_optimal_manifold
from heron.models.gpytorch import HeronNonSpinningApproximant
from heron.models.lalsimulation import SEOBNRv3, IMRPhenomPv2
from heron.models.lalnoise import AdvancedLIGO
from heron.filters import Overlap
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use("agg")


class _GenericWaveform(unittest.TestCase):
    """
    Test a generic waveform interface.

    This should be inherited by other specific classes/
    """

    @classmethod
    def setUpClass(cls):
        # initialize likelihood and model
        cls.model = IMRPhenomPv2()

    
    def test_time_series(self):
        """Check that the model produces a time series at all."""

        time_domain = self.model.time_domain(parameters={"m1": 10 * u.solMass, "m2": 10 * u.solMass})

        self.assertTrue("plus" in time_domain.waveforms)
        self.assertTrue("cross" in time_domain.waveforms)

    def test_time_axis_epoch_not_set(self):
        """Check that the correct time axis is produced"""
        time_domain = self.model.time_domain(parameters={"m1": 10 * u.solMass, "m2": 10 * u.solMass})
        epoch = 0
        self.assertTrue(np.min(np.abs(time_domain['plus'].times.value - epoch)) < 0.001)
        # Test that the epoch isn't at the start of the timeseries
        self.assertTrue(np.argmin(np.abs(time_domain['plus'].times.value - epoch)) > 10)
        f, ax = plt.subplots(1,1)
        ax.plot(np.sqrt(np.array(time_domain['plus'].data)**2 + np.array(time_domain['cross'].data)**2))
        f.savefig("test_time_axis_epoch_not_set.png")

        
    def test_time_axis_epoch_set(self):
        """Check that the correct time axis is produced"""
        epoch = 4000
        time_domain = self.model.time_domain(parameters={"m1": 10 * u.solMass, "m2": 10 * u.solMass, "gpstime": epoch})
        self.assertTrue(np.min(np.abs(time_domain['plus'].times.value - epoch)) < 0.001)
        # Test that the epoch isn't at the start of the timeseries
        self.assertTrue(np.argmin(np.abs(time_domain['plus'].times.value - epoch)) > 10)        

    def test_time_axis_axis_specified(self):
        """Check that the axis is correct if it has been provided to the call."""
        times = np.linspace(3999.95, 4000.01, 150)
        epoch = 4000
        time_domain = self.model.time_domain(parameters={"m1": 10 * u.solMass, "m2": 10 * u.solMass, "gpstime": epoch}, times=times)
        self.assertTrue(np.min(np.abs(time_domain['plus'].times.value - epoch)) < 0.001)
        # Test that the epoch isn't at the start of the timeseries
        self.assertTrue(np.argmin(np.abs(time_domain['plus'].times.value - epoch)) > 10)        
        # Check that the time axis is the same as the input one
        self.assertTrue(np.all(times == time_domain['plus'].times.value))

    def test_time_axis_axis_constructed(self):
        """Check that the axis is correct if it has been provided to the call in the parameters."""
        times = np.linspace(3999.95, 4000.01, 150)
        epoch = 4000
        time_domain = self.model.time_domain(parameters={"m1": 10 * u.solMass, "m2": 10 * u.solMass, "gpstime": epoch,
                                                         "time": {"upper": 4000.01, "lower": 3999.95, "number": 150}})
        self.assertTrue(np.min(np.abs(time_domain['plus'].times.value - epoch)) < 0.001)
        # Test that the epoch isn't at the start of the timeseries
        self.assertTrue(np.argmin(np.abs(time_domain['plus'].times.value - epoch)) > 10)        
        # Check that the time axis is the same as the input one
        self.assertTrue(np.all(times == time_domain['plus'].times.value))

    def test_waveform_epoch_not_set(self):
        """Test that the waveform looks sane if the epoch is not set."""
        
        time_domain = self.model.time_domain(parameters={"m1": 10 * u.solMass, "m2": 10 * u.solMass})
        epoch = 0
        t0_ix = np.argmin(np.abs(time_domain['plus'].times.value - epoch))
        hp_ix = np.argmax(np.sqrt(np.array(time_domain['plus'].data)**2 + np.array(time_domain['cross'].data)**2))

        f, ax = plt.subplots(1,1)
        ax.plot(np.array(time_domain['plus'].data))
        f.savefig("test_waveform_epoch_not_set.png")
        
        self.assertTrue(np.abs(t0_ix - hp_ix)<10)
        self.assertTrue(hp_ix + 30 < len(time_domain['plus'].times.value))

    def test_waveform_epoch_set(self):
        """Test that the waveform looks sane if the epoch is not set."""
        epoch = 5009.2        
        time_domain = self.model.time_domain(parameters={"m1": 10 * u.solMass, "m2": 10 * u.solMass, "gpstime": epoch})

        t0_ix = np.argmin(np.abs(time_domain['plus'].times.value - epoch))
        hp_ix = np.argmax(np.sqrt(np.array(time_domain['plus'].data)**2 + np.array(time_domain['cross'].data)**2))

        f, ax = plt.subplots(1,1)
        ax.plot(np.array(time_domain['plus'].data))
        f.savefig("test_waveform_epoch_not_set.png")
        
        self.assertTrue(np.abs(t0_ix - hp_ix)<10)
        self.assertTrue(hp_ix + 30 < len(time_domain['plus'].times.value))

    def test_waveform_axis_specified(self):
        """Check that the axis is correct if it has been provided to the call."""
        times = np.linspace(3999.95, 4000.01, 150)
        epoch = 4000
        time_domain = self.model.time_domain(parameters={"m1": 10 * u.solMass, "m2": 10 * u.solMass, "gpstime": epoch}, times=times)
        t0_ix = np.argmin(np.abs(time_domain['plus'].times.value - epoch))
        hp_ix = np.argmax(np.sqrt(np.array(time_domain['plus'].data)**2 + np.array(time_domain['cross'].data)**2))
        f, ax = plt.subplots(1,1)
        ax.plot(np.array(time_domain['plus'].data))
        f.savefig("test_waveform_axis_specified.png")

        
        self.assertTrue(np.abs(t0_ix - hp_ix)<10)
        self.assertTrue(hp_ix + 5 < len(time_domain['plus'].times.value))

    def test_waveform_axis_constructed(self):
        """Check that the axis is correct if it has been provided to the call."""
        times = np.linspace(3999.95, 4000.01, 150)
        epoch = 4000
        time_domain = self.model.time_domain(parameters={"m1": 10 * u.solMass, "m2": 10 * u.solMass, "gpstime": epoch,
                                                         "time": {"upper": 4000.01, "lower": 3999.95, "number": 150}})
        t0_ix = np.argmin(np.abs(time_domain['plus'].times.value - epoch))
        hp_ix = np.argmax(np.sqrt(np.array(time_domain['plus'].data)**2 + np.array(time_domain['cross'].data)**2))
        f, ax = plt.subplots(1,1)
        ax.plot(np.array(time_domain['plus'].data))
        f.savefig("test_waveform_axis_constructed.png")
        
        self.assertTrue(np.abs(t0_ix - hp_ix)<10)
        self.assertTrue(hp_ix + 5 < len(time_domain['plus'].times.value))

    def test_waveform_mass_ratio_total_mass(self):
        """Test that a waveform can be constructed if the mass ratio and total mass are provided."""
        epoch = 4000
        time_domain = self.model.time_domain(parameters={"total_mass": 20 * u.solMass, "mass_ratio": 0.5})
        # print(self.model.args)
        t0_ix = np.argmin(np.abs(time_domain['plus'].times.value - epoch))
        hp_ix = np.argmax(np.sqrt(np.array(time_domain['plus'].data)**2 + np.array(time_domain['cross'].data)**2))

        f, ax = plt.subplots(1,1)
        ax.plot(time_domain['plus'].times, np.array(time_domain['plus'].data))
        f.savefig("test_waveform_mass_ratio.png")
        self.assertFalse(self.model.args['m1'] == 10)
        self.assertFalse(self.model.args['m1'] == self.model.args['m2'])
        self.assertTrue(hp_ix + 100 < len(time_domain['plus'].times.value))

    def test_waveform_mass_ratio_chirp_mass(self):
        """Test that a waveform can be constructed if the mass ratio and total mass are provided."""
        epoch = 4000
        time_domain = self.model.time_domain(parameters={"chirp_mass": 20 * u.solMass, "mass_ratio": 0.5})
        t0_ix = np.argmin(np.abs(time_domain['plus'].times.value - epoch))
        hp_ix = np.argmax(np.sqrt(np.array(time_domain['plus'].data)**2 + np.array(time_domain['cross'].data)**2))

        f, ax = plt.subplots(1,1)
        ax.plot(time_domain['plus'].times, np.array(time_domain['plus'].data))
        f.savefig("test_waveform_mass_ratio.png")
        self.assertFalse(self.model.args['m1'] == 10)
        self.assertFalse(self.model.args['m1'] == self.model.args['m2'])
        self.assertTrue(hp_ix + 100 < len(time_domain['plus'].times.value))

    def test_waveform_mass_no_units(self):
        """Test that the correct waveform is generated if no units are supplied."""
        mass1 = 40 * u.solMass
        mass2 = 20 * u.solMass

        time_domain_1 = self.model.time_domain(parameters={"mass_1": mass1,
                                                           "mass_2": mass2})
        time_domain_2 = self.model.time_domain(parameters={"mass_1": mass1.to(u.kilogram).value,
                                                           "mass_2": mass2.to(u.kilogram).value})

        self.assertTrue(time_domain_1['plus'].data == time_domain_2['plus'].data)
