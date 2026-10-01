import unittest

import numpy as np

from minke.models.lalnoise import (
    KNOWN_PSDS,
    AdvancedLIGOO4Sensitivity,
    AdvancedVirgoO4Sensitivity,
)


class TestAdvancedLIGOAndVirgoO4AreDistinct(unittest.TestCase):
    """Regression test for a real bug: lalsimulation.SimNoisePSDaLIGOAdVO4T1800545
    looks like an aLIGO curve by name but is a deprecated alias for Virgo's O4
    design curve, so AdvancedLIGOO4Sensitivity and AdvancedVirgoO4Sensitivity
    were previously (silently) bytewise identical."""

    def test_ligo_and_virgo_o4_psds_are_not_identical(self):
        frequencies = np.arange(20, 1024, 1)
        ligo = np.asarray(
            AdvancedLIGOO4Sensitivity().frequency_domain(frequencies=frequencies).data
        )
        virgo = np.asarray(
            AdvancedVirgoO4Sensitivity().frequency_domain(frequencies=frequencies).data
        )
        self.assertFalse(np.array_equal(ligo, virgo))


class TestAdvancedVirgoO4Sensitivity(unittest.TestCase):
    def setUp(self):
        self.psd = AdvancedVirgoO4Sensitivity()

    def test_registered_in_known_psds(self):
        self.assertIs(KNOWN_PSDS["AdvancedVirgo_O4"], AdvancedVirgoO4Sensitivity)

    def test_frequency_domain_is_finite_and_positive(self):
        # The wrapper's very last frequency bin is a known LAL-series
        # boundary artifact (always 0, regardless of which PSD model is
        # used) -- excluded here, not specific to this PSD.
        frequencies = np.arange(20, 1024, 1)
        psd = self.psd.frequency_domain(frequencies=frequencies)
        data = np.asarray(psd.data)[:-1]
        self.assertTrue(np.all(np.isfinite(data)))
        self.assertTrue(np.all(data > 0))

    def test_sensitive_band_is_physically_reasonable(self):
        # Virgo's amplitude spectral density in its most sensitive band
        # (roughly 50-300 Hz) should be on the order of 1e-23 - 1e-22
        # strain/sqrt(Hz), i.e. PSD ~ 1e-46 - 1e-44.
        frequencies = np.arange(50, 300, 1)
        psd = self.psd.frequency_domain(frequencies=frequencies)
        data = np.asarray(psd.data)[:-1]
        self.assertTrue(np.all(data > 1e-48))
        self.assertTrue(np.all(data < 1e-42))


if __name__ == "__main__":
    unittest.main()
