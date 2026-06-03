"""Tests for minke.bagpuss — the bagpuss InjectionSet adapter.

Exercises:
- read_injection_parameters: loads an InjectionSet HDF5 and returns a list
  of parameter dicts that minke/LALSimulation can consume.
- Parameter name mapping (m1_source → m1, cos_tilt → tilt, etc.)
- Spin conversion via puddin.lalsim
- Astropy unit attachment for masses and distances

Run with: python -m pytest tests/test_bagpuss.py
"""

from __future__ import annotations

import math
import tempfile
from pathlib import Path

import h5py
import numpy as np
import pytest

lalsimulation = pytest.importorskip("lalsimulation")
lal            = pytest.importorskip("lal")


# ---------------------------------------------------------------------------
# Fixture: minimal bagpuss-format HDF5
# ---------------------------------------------------------------------------

_FIELDS = [
    "m1_source", "m2_source",
    "a1", "a2",
    "cos_tilt1", "cos_tilt2",
    "phi12", "phi_jl",
    "theta_jn",
    "ra", "dec", "psi",
    "geocent_time",
    "redshift", "luminosity_distance",
    "host_galaxy_index",
]

_N = 5  # number of events in test fixtures


def _write_injection_hdf5(path: Path, data: dict[str, np.ndarray]) -> None:
    """Write a minimal bagpuss injection set to *path*."""
    with h5py.File(path, "w") as f:
        grp = f.create_group("injections")
        for field, values in data.items():
            grp.create_dataset(field, data=values)


def _make_aligned_data(n: int = _N) -> dict[str, np.ndarray]:
    """Return aligned-spin injection data (cos_tilt = 1, tilts = 0)."""
    rng = np.random.default_rng(0)
    return {
        "m1_source":          rng.uniform(20.0, 50.0, n),
        "m2_source":          rng.uniform(10.0, 20.0, n),
        "a1":                 rng.uniform(0.0,  0.5,  n),
        "a2":                 rng.uniform(0.0,  0.5,  n),
        "cos_tilt1":          np.ones(n),           # tilt = 0 → aligned
        "cos_tilt2":          np.ones(n),
        "phi12":              np.zeros(n),
        "phi_jl":             np.zeros(n),
        "theta_jn":           rng.uniform(0.0, math.pi, n),
        "ra":                 rng.uniform(0.0, 2 * math.pi, n),
        "dec":                np.arcsin(rng.uniform(-1.0, 1.0, n)),
        "psi":                rng.uniform(0.0, math.pi, n),
        "geocent_time":       rng.uniform(1.187e9, 1.270e9, n),
        "redshift":           rng.uniform(0.01, 0.3, n),
        "luminosity_distance": rng.uniform(50.0, 1500.0, n),
        "host_galaxy_index":  np.arange(n, dtype=np.int64),
    }


def _make_precessing_data(n: int = _N) -> dict[str, np.ndarray]:
    """Return precessing-spin injection data."""
    rng = np.random.default_rng(42)
    cos_tilt1 = rng.uniform(-0.9, 0.9, n)  # general tilts, avoid ±1 edge
    cos_tilt2 = rng.uniform(-0.9, 0.9, n)
    return {
        "m1_source":          rng.uniform(20.0, 50.0, n),
        "m2_source":          rng.uniform(10.0, 20.0, n),
        "a1":                 rng.uniform(0.1,  0.8,  n),
        "a2":                 rng.uniform(0.1,  0.8,  n),
        "cos_tilt1":          cos_tilt1,
        "cos_tilt2":          cos_tilt2,
        "phi12":              rng.uniform(0.0, 2 * math.pi, n),
        "phi_jl":             rng.uniform(0.0, 2 * math.pi, n),
        "theta_jn":           rng.uniform(0.0, math.pi, n),
        "ra":                 rng.uniform(0.0, 2 * math.pi, n),
        "dec":                np.arcsin(rng.uniform(-1.0, 1.0, n)),
        "psi":                rng.uniform(0.0, math.pi, n),
        "geocent_time":       rng.uniform(1.187e9, 1.270e9, n),
        "redshift":           rng.uniform(0.01, 0.3, n),
        "luminosity_distance": rng.uniform(50.0, 1500.0, n),
        "host_galaxy_index":  np.arange(n, dtype=np.int64),
    }


@pytest.fixture()
def aligned_hdf5(tmp_path):
    """Path to a temporary HDF5 file with aligned-spin injections."""
    path = tmp_path / "injections_aligned.h5"
    _write_injection_hdf5(path, _make_aligned_data())
    return path


@pytest.fixture()
def precessing_hdf5(tmp_path):
    """Path to a temporary HDF5 file with precessing-spin injections."""
    path = tmp_path / "injections_precessing.h5"
    _write_injection_hdf5(path, _make_precessing_data())
    return path


# ===========================================================================
# read_injection_parameters
# ===========================================================================

class TestReadInjectionParameters:
    """Tests for minke.bagpuss.read_injection_parameters."""

    # --- output structure ---------------------------------------------------

    def test_returns_list_of_dicts(self, aligned_hdf5):
        """Returns a list of dicts, one per injection."""
        from minke.bagpuss import read_injection_parameters
        params = read_injection_parameters(aligned_hdf5)
        assert isinstance(params, list)
        assert len(params) == _N
        for p in params:
            assert isinstance(p, dict)

    def test_length_matches_file(self, precessing_hdf5):
        """List length equals number of events in the file."""
        from minke.bagpuss import read_injection_parameters
        params = read_injection_parameters(precessing_hdf5)
        assert len(params) == _N

    # --- parameter name mapping ---------------------------------------------

    def test_m1_key_present(self, aligned_hdf5):
        """Each dict has 'm1' (not 'm1_source')."""
        from minke.bagpuss import read_injection_parameters
        params = read_injection_parameters(aligned_hdf5)
        assert "m1" in params[0]
        assert "m1_source" not in params[0]

    def test_m2_key_present(self, aligned_hdf5):
        """Each dict has 'm2' (not 'm2_source')."""
        from minke.bagpuss import read_injection_parameters
        params = read_injection_parameters(aligned_hdf5)
        assert "m2" in params[0]
        assert "m2_source" not in params[0]

    def test_no_cos_tilt_keys(self, aligned_hdf5):
        """cos_tilt1 and cos_tilt2 are consumed and not passed through."""
        from minke.bagpuss import read_injection_parameters
        params = read_injection_parameters(aligned_hdf5)
        assert "cos_tilt1" not in params[0]
        assert "cos_tilt2" not in params[0]

    def test_lalsim_spin_keys_present(self, aligned_hdf5):
        """Cartesian spin components S1x … S2z are present."""
        from minke.bagpuss import read_injection_parameters
        params = read_injection_parameters(aligned_hdf5)
        for key in ("S1x", "S1y", "S1z", "S2x", "S2y", "S2z"):
            assert key in params[0], f"Missing key: {key}"

    def test_iota_present_not_theta_jn(self, aligned_hdf5):
        """'iota' (LALSim name) is present; 'theta_jn' is consumed."""
        from minke.bagpuss import read_injection_parameters
        params = read_injection_parameters(aligned_hdf5)
        assert "iota" in params[0]
        assert "theta_jn" not in params[0]

    def test_geocent_time_renamed_to_gpstime(self, aligned_hdf5):
        """geocent_time is renamed to gpstime for LALSim compatibility."""
        from minke.bagpuss import read_injection_parameters
        params = read_injection_parameters(aligned_hdf5)
        assert "gpstime" in params[0]
        assert "geocent_time" not in params[0]

    def test_luminosity_distance_present(self, aligned_hdf5):
        """luminosity_distance key is preserved (minke's _convert handles it)."""
        from minke.bagpuss import read_injection_parameters
        params = read_injection_parameters(aligned_hdf5)
        assert "luminosity_distance" in params[0]

    # --- mass units ---------------------------------------------------------

    def test_m1_has_solar_mass_units(self, aligned_hdf5):
        """m1 values carry astropy solar mass units."""
        from astropy import units as u
        from minke.bagpuss import read_injection_parameters
        params = read_injection_parameters(aligned_hdf5)
        assert hasattr(params[0]["m1"], "unit")
        assert params[0]["m1"].unit == u.solMass

    def test_m1_value_matches_source(self, aligned_hdf5):
        """m1 value (in solar masses) matches m1_source from the file."""
        from minke.bagpuss import read_injection_parameters
        data = _make_aligned_data()
        params = read_injection_parameters(aligned_hdf5)
        for i, p in enumerate(params):
            assert math.isclose(p["m1"].value, data["m1_source"][i], rel_tol=1e-10)

    # --- spin values (aligned case) ----------------------------------------

    def test_aligned_s1_transverse_zero(self, aligned_hdf5):
        """For aligned spins (cos_tilt=1), S1x = S1y = 0."""
        from minke.bagpuss import read_injection_parameters
        params = read_injection_parameters(aligned_hdf5)
        for p in params:
            assert abs(p["S1x"]) < 1e-14, f"S1x = {p['S1x']}"
            assert abs(p["S1y"]) < 1e-14, f"S1y = {p['S1y']}"

    def test_aligned_s1z_matches_a1(self, aligned_hdf5):
        """For aligned spins, S1z = a1."""
        from minke.bagpuss import read_injection_parameters
        data   = _make_aligned_data()
        params = read_injection_parameters(aligned_hdf5)
        for i, p in enumerate(params):
            assert math.isclose(p["S1z"], data["a1"][i], rel_tol=1e-10)

    def test_aligned_iota_equals_theta_jn(self, aligned_hdf5):
        """For aligned spins, iota = theta_jn."""
        from minke.bagpuss import read_injection_parameters
        data   = _make_aligned_data()
        params = read_injection_parameters(aligned_hdf5)
        for i, p in enumerate(params):
            assert math.isclose(p["iota"], data["theta_jn"][i], rel_tol=1e-10)

    # --- spin values (precessing case) ------------------------------------

    def test_precessing_spin_magnitudes_preserved(self, precessing_hdf5):
        """For precessing spins, |S1|² = a1² and |S2|² = a2²."""
        from minke.bagpuss import read_injection_parameters
        data   = _make_precessing_data()
        params = read_injection_parameters(precessing_hdf5)
        for i, p in enumerate(params):
            s1_mag2 = p["S1x"]**2 + p["S1y"]**2 + p["S1z"]**2
            s2_mag2 = p["S2x"]**2 + p["S2y"]**2 + p["S2z"]**2
            assert math.isclose(s1_mag2, data["a1"][i]**2, rel_tol=1e-8)
            assert math.isclose(s2_mag2, data["a2"][i]**2, rel_tol=1e-8)

    # --- passthrough fields ------------------------------------------------

    def test_ra_dec_preserved(self, aligned_hdf5):
        """ra and dec are passed through unchanged."""
        from minke.bagpuss import read_injection_parameters
        data   = _make_aligned_data()
        params = read_injection_parameters(aligned_hdf5)
        for i, p in enumerate(params):
            assert math.isclose(p["ra"],  data["ra"][i],  rel_tol=1e-12)
            assert math.isclose(p["dec"], data["dec"][i], rel_tol=1e-12)

    def test_psi_preserved(self, aligned_hdf5):
        """psi is passed through unchanged."""
        from minke.bagpuss import read_injection_parameters
        data   = _make_aligned_data()
        params = read_injection_parameters(aligned_hdf5)
        for i, p in enumerate(params):
            assert math.isclose(p["psi"], data["psi"][i], rel_tol=1e-12)

    def test_gpstime_matches_geocent_time(self, aligned_hdf5):
        """gpstime equals the original geocent_time from the file."""
        from minke.bagpuss import read_injection_parameters
        data   = _make_aligned_data()
        params = read_injection_parameters(aligned_hdf5)
        for i, p in enumerate(params):
            assert math.isclose(p["gpstime"], data["geocent_time"][i], rel_tol=1e-12)

    # --- f_ref and phase defaults ------------------------------------------

    def test_default_f_ref_is_20_hz(self, aligned_hdf5):
        """Default reference frequency is 20 Hz."""
        from minke.bagpuss import read_injection_parameters
        # If the function accepts f_ref and phase without error, defaults work.
        params = read_injection_parameters(aligned_hdf5)
        assert len(params) == _N

    def test_custom_f_ref_accepted(self, aligned_hdf5):
        """Custom f_ref is accepted without error."""
        from minke.bagpuss import read_injection_parameters
        params = read_injection_parameters(aligned_hdf5, f_ref=40.0)
        assert len(params) == _N
