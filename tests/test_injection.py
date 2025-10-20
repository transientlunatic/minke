#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
test_injection
----------------------------------

Tests for `minke.injection` module.
"""

import unittest
import numpy as np
import astropy.units as u

from minke.injection import (
    calculate_network_snr_for_distance,
    find_distance_for_network_snr,
    make_injection,
    injection_parameters_add_units
)
from minke.models.lalsimulation import IMRPhenomXPHM
from minke.models.lalnoise import KNOWN_PSDS
from minke.detector import KNOWN_IFOS


class TestSNRCalculations(unittest.TestCase):
    """Tests for SNR calculation functions."""

    def setUp(self):
        """Set up test fixtures."""
        self.waveform_model = IMRPhenomXPHM()
        self.parameters = {
            'ra': 0,
            'dec': 0,
            'psi': 0,
            'theta_jn': 0,
            'phase': 0,
            'm1': 30 * u.solMass,
            'm2': 30 * u.solMass,
            'luminosity_distance': 100 * u.megaparsec
        }
        self.sample_rate = 4096
        self.duration = 4
        self.epoch = 0
        
        # Create detector and PSD objects
        self.detectors = [KNOWN_IFOS['AdvancedLIGOHanford'](), KNOWN_IFOS['AdvancedLIGOLivingston']()]
        self.psd_models = [
            KNOWN_PSDS['AdvancedLIGO'](),
            KNOWN_PSDS['AdvancedLIGO']()
        ]
        
        # Get times array
        psd_temp = self.psd_models[0].time_series(
            duration=self.duration,
            sample_rate=self.sample_rate,
            epoch=self.epoch
        )
        self.times = psd_temp.times

    def tearDown(self):
        """Tear down test fixtures."""
        pass

    def test_calculate_network_snr_for_distance(self):
        """Test that network SNR calculation returns a positive value."""
        distance = 100  # Mpc
        snr = calculate_network_snr_for_distance(
            distance,
            self.waveform_model,
            self.parameters,
            self.detectors,
            self.psd_models,
            self.times
        )
        
        self.assertIsInstance(snr, (float, np.floating))
        self.assertGreater(snr, 0)
        
    def test_network_snr_decreases_with_distance(self):
        """Test that SNR decreases as distance increases."""
        distance_near = 50  # Mpc
        distance_far = 200  # Mpc
        
        snr_near = calculate_network_snr_for_distance(
            distance_near,
            self.waveform_model,
            self.parameters,
            self.detectors,
            self.psd_models,
            self.times
        )
        
        snr_far = calculate_network_snr_for_distance(
            distance_far,
            self.waveform_model,
            self.parameters,
            self.detectors,
            self.psd_models,
            self.times
        )
        
        self.assertGreater(snr_near, snr_far)
        
    def test_network_snr_scales_inversely_with_distance(self):
        """Test that SNR approximately scales as 1/distance."""
        distance1 = 100  # Mpc
        distance2 = 200  # Mpc
        
        snr1 = calculate_network_snr_for_distance(
            distance1,
            self.waveform_model,
            self.parameters,
            self.detectors,
            self.psd_models,
            self.times
        )
        
        snr2 = calculate_network_snr_for_distance(
            distance2,
            self.waveform_model,
            self.parameters,
            self.detectors,
            self.psd_models,
            self.times
        )
        
        # SNR should be approximately halved when distance is doubled
        ratio = snr1 / snr2
        self.assertAlmostEqual(ratio, 2.0, delta=0.1)
        
    def test_find_distance_for_network_snr(self):
        """Test finding distance for a target SNR."""
        target_snr = 10.0
        
        distance = find_distance_for_network_snr(
            target_snr,
            self.waveform_model,
            self.parameters,
            self.detectors,
            self.psd_models,
            self.times
        )
        
        self.assertIsInstance(distance, u.Quantity)
        self.assertEqual(distance.unit, u.megaparsec)
        self.assertGreater(distance.value, 10)
        self.assertLess(distance.value, 10000)
        
        # Verify that the found distance gives the target SNR
        actual_snr = calculate_network_snr_for_distance(
            distance.value,
            self.waveform_model,
            self.parameters,
            self.detectors,
            self.psd_models,
            self.times
        )
        
        self.assertAlmostEqual(actual_snr, target_snr, delta=0.2)
        
    def test_find_distance_for_high_snr(self):
        """Test finding distance for a high SNR value."""
        target_snr = 100.0
        
        distance = find_distance_for_network_snr(
            target_snr,
            self.waveform_model,
            self.parameters,
            self.detectors,
            self.psd_models,
            self.times
        )
        
        # High SNR should require close distance
        self.assertLess(distance.value, 100)
        
    def test_find_distance_for_low_snr(self):
        """Test finding distance for a low SNR value."""
        target_snr = 5.0
        
        distance = find_distance_for_network_snr(
            target_snr,
            self.waveform_model,
            self.parameters,
            self.detectors,
            self.psd_models,
            self.times
        )
        
        # Low SNR should allow larger distance
        self.assertGreater(distance.value, 50)


class TestMakeInjection(unittest.TestCase):
    """Tests for the make_injection function."""

    def setUp(self):
        """Set up test fixtures."""
        self.basic_parameters = {
            'm1': 30,
            'm2': 30,
        }
        self.detectors = {
            'AdvancedLIGOHanford': 'AdvancedLIGO',
            'AdvancedLIGOLivingston': 'AdvancedLIGO'
        }
        self.duration = 4
        self.sample_rate = 4096
        self.epoch = 0

    def tearDown(self):
        """Tear down test fixtures."""
        pass

    def test_make_injection_with_distance(self):
        """Test making an injection with specified luminosity distance."""
        parameters = self.basic_parameters.copy()
        parameters['luminosity_distance'] = 100
        
        injections = make_injection(
            waveform=IMRPhenomXPHM,
            injection_parameters=parameters,
            detectors=self.detectors,
            duration=self.duration,
            sample_rate=self.sample_rate,
            epoch=self.epoch,
            channel="TestInjection"
        )
        
        self.assertEqual(len(injections), 2)
        self.assertIn('H1', injections)
        self.assertIn('L1', injections)
        
        # Check that injections have the correct length
        for det_name, injection in injections.items():
            expected_length = self.duration * self.sample_rate
            self.assertEqual(len(injection.data), expected_length)

    def test_make_injection_with_snr(self):
        """Test making an injection with specified target SNR."""
        parameters = self.basic_parameters.copy()
        parameters['snr'] = 20.0
        
        injections = make_injection(
            waveform=IMRPhenomXPHM,
            injection_parameters=parameters,
            detectors=self.detectors,
            duration=self.duration,
            sample_rate=self.sample_rate,
            epoch=self.epoch,
            channel="TestInjection"
        )
        
        self.assertEqual(len(injections), 2)
        self.assertIn('H1', injections)
        self.assertIn('L1', injections)
        
        # The function should have calculated a luminosity distance
        # We can't easily verify the exact SNR without recalculating,
        # but we can check that injections were created successfully
        for det_name, injection in injections.items():
            expected_length = self.duration * self.sample_rate
            self.assertEqual(len(injection.data), expected_length)
            
    def test_make_injection_with_custom_times(self):
        """Test making an injection with custom time array."""
        parameters = self.basic_parameters.copy()
        parameters['luminosity_distance'] = 100
        
        times = np.linspace(0, self.duration, self.duration * self.sample_rate)
        
        injections = make_injection(
            waveform=IMRPhenomXPHM,
            injection_parameters=parameters,
            detectors=self.detectors,
            times=times,
            channel="TestInjection"
        )
        
        self.assertEqual(len(injections), 2)
        
    def test_make_injection_single_detector(self):
        """Test making an injection for a single detector."""
        parameters = self.basic_parameters.copy()
        parameters['luminosity_distance'] = 100
        
        single_detector = {'AdvancedLIGOHanford': 'AdvancedLIGO'}
        
        injections = make_injection(
            waveform=IMRPhenomXPHM,
            injection_parameters=parameters,
            detectors=single_detector,
            duration=self.duration,
            sample_rate=self.sample_rate,
            epoch=self.epoch,
            channel="TestInjection"
        )
        
        self.assertEqual(len(injections), 1)
        self.assertIn('H1', injections)


class TestInjectionParametersAddUnits(unittest.TestCase):
    """Tests for the injection_parameters_add_units function."""

    def setUp(self):
        """Set up test fixtures."""
        pass

    def tearDown(self):
        """Tear down test fixtures."""
        pass

    def test_add_units_to_masses(self):
        """Test adding units to mass parameters."""
        parameters = {
            'm1': 30,
            'm2': 25
        }
        
        result = injection_parameters_add_units(parameters)
        
        self.assertIsInstance(result['m1'], u.Quantity)
        self.assertIsInstance(result['m2'], u.Quantity)
        self.assertEqual(result['m1'].unit, u.solMass)
        self.assertEqual(result['m2'].unit, u.solMass)
        self.assertEqual(result['m1'].value, 30)
        self.assertEqual(result['m2'].value, 25)

    def test_add_units_to_distance(self):
        """Test adding units to luminosity distance."""
        parameters = {
            'luminosity_distance': 100
        }
        
        result = injection_parameters_add_units(parameters)
        
        self.assertIsInstance(result['luminosity_distance'], u.Quantity)
        self.assertEqual(result['luminosity_distance'].unit, u.megaparsec)
        self.assertEqual(result['luminosity_distance'].value, 100)

    def test_preserve_existing_units(self):
        """Test that existing units are preserved."""
        parameters = {
            'm1': 30 * u.solMass,
            'luminosity_distance': 100 * u.megaparsec
        }
        
        result = injection_parameters_add_units(parameters)
        
        self.assertEqual(result['m1'], 30 * u.solMass)
        self.assertEqual(result['luminosity_distance'], 100 * u.megaparsec)

    def test_parameters_without_units(self):
        """Test parameters that don't need units."""
        parameters = {
            'ra': 0,
            'dec': 0,
            'psi': 0.5
        }
        
        result = injection_parameters_add_units(parameters)
        
        # These parameters should not have units added
        self.assertNotIsInstance(result['ra'], u.Quantity)
        self.assertNotIsInstance(result['dec'], u.Quantity)
        self.assertNotIsInstance(result['psi'], u.Quantity)

    def test_mixed_parameters(self):
        """Test a mix of parameters with and without units."""
        parameters = {
            'm1': 30,
            'm2': 25,
            'luminosity_distance': 100,
            'ra': 1.5,
            'dec': -0.3,
            'phase': 0
        }
        
        result = injection_parameters_add_units(parameters)
        
        # Check that correct parameters have units
        self.assertIsInstance(result['m1'], u.Quantity)
        self.assertIsInstance(result['m2'], u.Quantity)
        self.assertIsInstance(result['luminosity_distance'], u.Quantity)
        
        # Check that other parameters don't have units
        self.assertNotIsInstance(result['ra'], u.Quantity)
        self.assertNotIsInstance(result['dec'], u.Quantity)
        self.assertNotIsInstance(result['phase'], u.Quantity)


if __name__ == '__main__':
    unittest.main()
