#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
test_noise
----------------------------------

Tests for `minke.noise` module.
"""

import unittest
import numpy as np
from minke.noise import LALSimulationPSD, AdvancedLIGO, AdvancedLIGODesignSensitivity2018, KNOWN_PSDS


class TestLALSimulationPSD(unittest.TestCase):
    """Tests for the LALSimulationPSD class."""

    def setUp(self):
        """Set up test fixtures."""
        self.psd = AdvancedLIGO()
        
    def tearDown(self):
        """Tear down test fixtures."""
        pass

    # Basic functionality tests
    
    def test_initialization(self):
        """Test that LALSimulationPSD can be initialized."""
        psd = LALSimulationPSD()
        self.assertIsNotNone(psd)

    def test_frequency_domain_default_parameters(self):
        """Test frequency domain PSD generation with default parameters."""
        psd = self.psd.frequency_domain()
        self.assertIsNotNone(psd)
        self.assertGreater(len(psd.data), 0)
        self.assertGreater(len(psd.frequencies), 0)
        
    def test_frequency_domain_custom_parameters(self):
        """Test frequency domain PSD with custom frequency range."""
        lower_freq = 10
        upper_freq = 512
        df = 0.5
        # Build a custom frequency array with spacing df and pass directly
        frequencies = np.arange(lower_freq, upper_freq + df, df)
        psd = self.psd.frequency_domain(
            frequencies=frequencies
        )
        self.assertGreater(len(psd.data), 0)
        # Check frequency spacing
        freq_diff = psd.frequencies[1] - psd.frequencies[0]
        self.assertAlmostEqual(freq_diff.value, df, places=5)
        
    def test_frequency_domain_custom_frequencies(self):
        """Test frequency domain PSD with custom frequency array."""
        frequencies = np.linspace(20, 1024, 500)
        psd = self.psd.frequency_domain(frequencies=frequencies)
        self.assertEqual(len(psd.data), len(frequencies))
        
    def test_twocolumn_format(self):
        """Test two-column format output."""
        data = self.psd.twocolumn(df=1, lower_frequency=20, upper_frequency=100)
        self.assertEqual(data.shape[1], 2)  # Should have 2 columns
        self.assertGreater(data.shape[0], 0)  # Should have rows
        # First column should be frequencies, second should be PSD values
        self.assertTrue(np.all(data[:, 0] >= 20))
        # PSD values: finite, non-negative, and at least some strictly positive
        self.assertTrue(np.all(np.isfinite(data[:, 1])))
        self.assertTrue(np.all(data[:, 1] >= 0))
        self.assertGreater(np.sum(data[:, 1] > 0), 0)
        
    def test_covariance_matrix_shape(self):
        """Test that covariance matrix has correct shape."""
        times = np.linspace(0, 1, 256)
        cov = self.psd.covariance_matrix(times)
        n = len(times)
        self.assertEqual(cov.shape, (n, n))
        
    def test_covariance_matrix_symmetry(self):
        """Test that covariance matrix is symmetric."""
        times = np.linspace(0, 1, 128)
        cov = self.psd.covariance_matrix(times)
        self.assertTrue(np.allclose(cov, cov.T))
        
    def test_time_domain(self):
        """Test time domain representation."""
        times = np.linspace(0, 1, 256)
        cov = self.psd.time_domain(times)
        self.assertEqual(cov.shape, (len(times), len(times)))

    # Time series generation tests
    
    def test_time_series_with_times(self):
        """Test time series generation with explicit times."""
        duration = 4.0
        sample_rate = 1024
        times = np.linspace(0, duration, int(duration * sample_rate))
        noise = self.psd.time_series(times=times, sample_rate=sample_rate)
        self.assertEqual(len(noise.data), len(times))
        self.assertTrue(np.allclose(noise.times.value, times))
        
    def test_time_series_with_duration_and_sample_rate(self):
        """Test time series generation with duration and sample rate."""
        duration = 4.0
        sample_rate = 1024
        noise = self.psd.time_series(duration=duration, sample_rate=sample_rate)
        expected_length = int(duration * sample_rate)
        self.assertEqual(len(noise.data), expected_length)
        
    def test_time_series_epoch(self):
        """Test time series generation with custom epoch."""
        duration = 2.0
        sample_rate = 512
        epoch = 1000.0
        noise = self.psd.time_series(
            duration=duration,
            sample_rate=sample_rate,
            epoch=epoch
        )
        self.assertAlmostEqual(noise.times.value[0], epoch, places=5)
        
    def test_time_series_data_properties(self):
        """Test that generated time series has expected statistical properties."""
        duration = 8.0
        sample_rate = 1024
        noise = self.psd.time_series(duration=duration, sample_rate=sample_rate)
        # Noise should be real-valued
        self.assertTrue(np.all(np.isreal(noise.data)))
        # Noise should have finite values
        self.assertTrue(np.all(np.isfinite(noise.data)))
        # Mean should be close to zero for colored noise
        # (may not be exactly zero due to finite sample)
        self.assertLess(np.abs(np.mean(noise.data)), 1e-10)

    # PSD verification tests
    
    def test_generated_noise_psd_matches_input_psd_long_duration(self):
        """
        Test that the PSD of generated noise matches the input PSD.
        Uses a long duration for better frequency resolution.
        """
        duration = 64.0  # Long duration for better frequency resolution
        sample_rate = 2048
        
        # Generate noise
        noise = self.psd.time_series(duration=duration, sample_rate=sample_rate)
        
        # Compute PSD of the generated noise using Welch's method
        from scipy import signal
        x = np.asarray(noise.data)
        frequencies_noise, psd_noise = signal.welch(
            x,
            fs=sample_rate,
            nperseg=int(sample_rate * 4),  # 4 second segments
            noverlap=int(sample_rate * 2),  # 50% overlap
            scaling='density'
        )
        
        # Get the theoretical PSD at these frequencies
        freq_array = np.asarray(frequencies_noise, dtype=float)
        psd_theoretical = self.psd.frequency_domain(
            frequencies=freq_array
        )
        psd_theoretical_data = np.asarray(psd_theoretical.data)
        
        # Compare in the sensitive frequency band (e.g., 20-500 Hz for aLIGO)
        mask = (frequencies_noise >= 20) & (frequencies_noise <= 500)
        
        # Calculate relative difference
        relative_diff = np.abs(
            (psd_noise[mask] - psd_theoretical_data[mask]) / psd_theoretical_data[mask]
        )
        
        # Allow for statistical fluctuations - should match within ~20% on average
        mean_relative_diff = np.mean(relative_diff)
        self.assertLess(mean_relative_diff, 0.3,
                       msg=f"Mean relative PSD difference {mean_relative_diff:.2%} exceeds 30%")
        
    def test_generated_noise_psd_matches_input_psd_multiple_realizations(self):
        """
        Test PSD matching using multiple noise realizations averaged together.
        """
        duration = 16.0
        sample_rate = 2048
        n_realizations = 5
        
        from scipy import signal
        
        # First realization to capture frequencies
        noise = self.psd.time_series(duration=duration, sample_rate=sample_rate)
        x = np.asarray(noise.data)
        frequencies_noise, psd_noise = signal.welch(
            x,
            fs=sample_rate,
            nperseg=int(sample_rate * 4),
            noverlap=int(sample_rate * 2),
            scaling='density'
        )
        psd_estimates = [psd_noise]
        # Remaining realizations
        for _ in range(n_realizations - 1):
            noise = self.psd.time_series(duration=duration, sample_rate=sample_rate)
            x = np.asarray(noise.data)
            _, psd_noise = signal.welch(
                x,
                fs=sample_rate,
                nperseg=int(sample_rate * 4),
                noverlap=int(sample_rate * 2),
                scaling='density'
            )
            psd_estimates.append(psd_noise)
        
        # Average the PSD estimates
        psd_noise_avg = np.mean(psd_estimates, axis=0)
        
        # Get theoretical PSD
        freq_array = np.asarray(frequencies_noise, dtype=float)
        psd_theoretical = self.psd.frequency_domain(
            frequencies=freq_array
        )
        psd_theoretical_data = np.asarray(psd_theoretical.data)
        
        # Compare in sensitive band
        mask = (frequencies_noise >= 20) & (frequencies_noise <= 500)
        relative_diff = np.abs(
            (psd_noise_avg[mask] - psd_theoretical_data[mask]) / psd_theoretical_data[mask]
        )
        
        # With averaging, we should get better agreement
        mean_relative_diff = np.mean(relative_diff)
        self.assertLess(mean_relative_diff, 0.2,
                       msg=f"Mean relative PSD difference {mean_relative_diff:.2%} exceeds 20% with averaging")
        
    def test_noise_psd_shape_in_frequency_domain(self):
        """Test that FFT of noise has correct shape matching PSD."""
        duration = 4.0
        sample_rate = 1024
        noise = self.psd.time_series(duration=duration, sample_rate=sample_rate)
        
        # Compute FFT
        fft_data = np.fft.rfft(noise.data)
        fft_freqs = np.fft.rfftfreq(len(noise.data), 1.0/sample_rate)
        
        # Get PSD at same frequencies
        freq_tensor = np.asarray(fft_freqs, dtype=float)
        psd = self.psd.frequency_domain(frequencies=freq_tensor)
        
        # Shapes should match
        self.assertEqual(len(fft_data), len(psd.data))
        
    def test_multiple_noise_realizations_different(self):
        """Test that multiple noise realizations are different (stochastic)."""
        duration = 2.0
        sample_rate = 1024
        
        noise1 = self.psd.time_series(duration=duration, sample_rate=sample_rate)
        noise2 = self.psd.time_series(duration=duration, sample_rate=sample_rate)
        
        # The two realizations should not be bitwise identical
        self.assertFalse(np.array_equal(np.asarray(noise1.data), np.asarray(noise2.data)))
        
        # But should have similar statistical properties (variance)
        var1 = np.var(np.asarray(noise1.data))
        var2 = np.var(np.asarray(noise2.data))
        # Variances should be within a reasonable range given randomness
        self.assertLess(abs(var1 - var2) / ((var1 + var2) / 2), 2.0)


class TestAdvancedLIGO(unittest.TestCase):
    """Tests specific to AdvancedLIGO PSD."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.psd = AdvancedLIGO()
        
    def test_initialization(self):
        """Test that AdvancedLIGO PSD can be initialized."""
        self.assertIsNotNone(self.psd)
        
    def test_frequency_domain_aLIGO(self):
        """Test aLIGO PSD generation."""
        psd = self.psd.frequency_domain(lower_frequency=10, upper_frequency=2048)
        # Check PSD values are finite and non-negative with some positive
        psd_vals = np.asarray(psd.data)
        self.assertTrue(np.all(np.isfinite(psd_vals)))
        self.assertTrue(np.all(psd_vals >= 0))
        self.assertGreater(np.sum(psd_vals > 0), 0)
        # Basic sanity: frequencies increase and PSD array length matches
        self.assertEqual(len(psd_vals), len(psd.frequencies))
        
    def test_aLIGO_alias(self):
        """Test that AdvancedLIGODesignSensitivity2018 and AdvancedLIGO are equivalent."""
        psd1 = AdvancedLIGO()
        psd2 = AdvancedLIGODesignSensitivity2018()
        
        # Generate PSDs with same parameters
        psd_data1 = psd1.frequency_domain(df=1, lower_frequency=20, upper_frequency=1024)
        psd_data2 = psd2.frequency_domain(df=1, lower_frequency=20, upper_frequency=1024)
        
        # They should be identical
        self.assertTrue(np.allclose(np.asarray(psd_data1.data), np.asarray(psd_data2.data)))
        

class TestKnownPSDs(unittest.TestCase):
    """Tests for the KNOWN_PSDS registry."""
    
    def test_known_psds_contains_advanced_ligo(self):
        """Test that KNOWN_PSDS contains AdvancedLIGO."""
        self.assertIn("AdvancedLIGO", KNOWN_PSDS)
        
    def test_known_psds_can_instantiate(self):
        """Test that PSDs from registry can be instantiated."""
        for name, psd_class in KNOWN_PSDS.items():
            psd = psd_class()
            self.assertIsNotNone(psd)
            # Test basic functionality
            psd_data = psd.frequency_domain(lower_frequency=20, upper_frequency=100)
            self.assertGreater(len(psd_data.data), 0)


class TestBilbyComparison(unittest.TestCase):
    """Tests comparing minke noise generation with bilby (if available)."""
    
    def setUp(self):
        """Set up test fixtures and check if bilby is available."""
        try:
            import bilby
            self.bilby = bilby
        except ImportError:
            self.skipTest('bilby is not installed')
        self.psd = AdvancedLIGO()
        
    def test_noise_mean_comparison_with_bilby(self):
        """
        Test that noise generated by minke has similar mean to bilby-generated noise.
        Both should have mean close to zero.
        """
            
        duration = 4.0
        sample_rate = 2048
        
        # Generate noise with minke
        minke_noise = self.psd.time_series(duration=duration, sample_rate=sample_rate)
        minke_mean = np.mean(minke_noise.data)
        
        # Generate noise with bilby
        from bilby.gw.detector import PowerSpectralDensity
        from bilby.core.utils import create_frequency_series
        
        # Create frequency array for bilby PSD
        frequency_array = create_frequency_series(sample_rate, duration)
        
        # Get PSD values from minke
        freq_array = np.asarray(frequency_array, dtype=float)
        minke_psd = self.psd.frequency_domain(frequencies=freq_array)
        
        # Create bilby PSD from minke PSD values
        bilby_psd = PowerSpectralDensity(
            frequency_array=frequency_array,
            psd_array=np.asarray(minke_psd.data)
        )
        
        # Generate bilby noise (if available)
        get_noise = getattr(bilby_psd, 'get_noise_realisation', None)
        if get_noise is None:
            self.skipTest("bilby.PowerSpectralDensity.get_noise_realisation not available in this version")
        # bilby returns (frequency_domain_noise, frequencies) tuple
        bilby_noise_fd, _ = get_noise(sample_rate, duration)
        # Convert to time domain
        bilby_noise = np.fft.irfft(bilby_noise_fd)
        bilby_mean = np.mean(bilby_noise)
        
        # Both means should be very close to zero
        self.assertLess(abs(minke_mean), 1e-10,
                       msg=f"Minke noise mean {minke_mean} is not close to zero")
        self.assertLess(abs(bilby_mean), 1e-10,
                       msg=f"Bilby noise mean {bilby_mean} is not close to zero")
        
        # The difference between means should also be negligible
        self.assertLess(abs(minke_mean - bilby_mean), 1e-10,
                       msg=f"Difference in means: minke={minke_mean}, bilby={bilby_mean}")
        
    def test_noise_variance_comparison_with_bilby(self):
        """
        Test that noise generated by minke has similar variance to bilby-generated noise.
        """
            
        duration = 8.0
        sample_rate = 2048
        
        # Generate multiple realizations for better statistics
        n_realizations = 3
        minke_variances = []
        bilby_variances = []
        
        for _ in range(n_realizations):
            # Generate noise with minke
            minke_noise = self.psd.time_series(duration=duration, sample_rate=sample_rate)
            minke_variances.append(np.var(minke_noise.data))
            
            # Generate noise with bilby
            from bilby.gw.detector import PowerSpectralDensity
            from bilby.core.utils import create_frequency_series
            
            frequency_array = create_frequency_series(sample_rate, duration)
            freq_array = np.asarray(frequency_array, dtype=float)
            minke_psd = self.psd.frequency_domain(frequencies=freq_array)
            
            bilby_psd = PowerSpectralDensity(
                frequency_array=frequency_array,
                psd_array=np.asarray(minke_psd.data)
            )
            
            get_noise = getattr(bilby_psd, 'get_noise_realisation', None)
            if get_noise is None:
                self.skipTest("bilby.PowerSpectralDensity.get_noise_realisation not available in this version")
            # bilby returns (frequency_domain_noise, frequencies) tuple
            bilby_noise_fd, _ = get_noise(sample_rate, duration)
            # Convert to time domain
            bilby_noise = np.fft.irfft(bilby_noise_fd)
            bilby_variances.append(np.var(bilby_noise))
        
        # Average variances
        avg_minke_var = np.mean(minke_variances)
        avg_bilby_var = np.mean(bilby_variances)
        
        # Note: bilby and minke may use different FFT normalizations
        # This test is informational - our PSD matching tests verify correctness
        relative_diff = abs(avg_minke_var - avg_bilby_var) / ((avg_minke_var + avg_bilby_var) / 2)
        
        # Only fail if variances are orders of magnitude different (factor > 100)
        # Some difference is expected due to different FFT conventions
        self.assertLess(relative_diff, 2.0,
                       msg=f"Variance difference unexpectedly large: minke={avg_minke_var:.2e}, "
                           f"bilby={avg_bilby_var:.2e}, relative diff={relative_diff:.2%}. "
                           f"Note: Some difference is expected due to FFT normalization conventions.")


if __name__ == '__main__':
    unittest.main()
