"""Tests for asimov blueprint generation from bagpuss injection sets.

Covers minke.bagpuss.make_blueprint and minke.bagpuss.write_blueprints.

Run with: python -m pytest tests/test_bagpuss_blueprints.py
"""

from __future__ import annotations

import math
import io
from pathlib import Path

import h5py
import numpy as np
import pytest
import yaml

lal = pytest.importorskip("lal")
lalsimulation = pytest.importorskip("lalsimulation")

MSUN_KG = lal.MSUN_SI
_N = 4


# ---------------------------------------------------------------------------
# Fixtures (reuse the helpers from test_bagpuss.py)
# ---------------------------------------------------------------------------

def _write_hdf5(path: Path, n: int = _N) -> dict:
    rng = np.random.default_rng(7)
    data = {
        "m1_source":          np.array([30.0, 20.0, 40.0, 15.0]),
        "m2_source":          np.array([20.0, 15.0, 30.0, 10.0]),
        "a1":                 np.zeros(n),
        "a2":                 np.zeros(n),
        "cos_tilt1":          np.ones(n),
        "cos_tilt2":          np.ones(n),
        "phi12":              np.zeros(n),
        "phi_jl":             np.zeros(n),
        "theta_jn":           rng.uniform(0.0, math.pi, n),
        "ra":                 rng.uniform(0.0, 2 * math.pi, n),
        "dec":                np.arcsin(rng.uniform(-1.0, 1.0, n)),
        "psi":                rng.uniform(0.0, math.pi, n),
        "geocent_time":       np.array([1.2e9, 1.2e9 + 100, 1.2e9 + 200, 1.2e9 + 300]),
        "redshift":           rng.uniform(0.01, 0.3, n),
        "luminosity_distance": rng.uniform(100.0, 1000.0, n),
        "host_galaxy_index":  np.arange(n, dtype=np.int64),
    }
    with h5py.File(path, "w") as f:
        grp = f.create_group("injections")
        for k, v in data.items():
            grp.create_dataset(k, data=v)
    return data


@pytest.fixture()
def injection_hdf5(tmp_path):
    path = tmp_path / "injections.h5"
    _write_hdf5(path)
    return path


@pytest.fixture()
def raw_data():
    return _write_hdf5(Path("/dev/null"))  # just the dict, no file needed


def _chirp_mass(m1_msun: float, m2_msun: float) -> float:
    """Analytic chirp mass in solar masses."""
    return (m1_msun * m2_msun) ** 0.6 / (m1_msun + m2_msun) ** 0.2


# ===========================================================================
# make_blueprint
# ===========================================================================

class TestMakeBlueprint:
    """Tests for minke.bagpuss.make_blueprint."""

    def _make_param(self, m1=30.0, m2=20.0, gpstime=1.2e9):
        """Return a minimal injection param dict."""
        from minke.bagpuss import read_injection_parameters
        # Build and read a single-event HDF5
        import tempfile, h5py, numpy as np
        rng = np.random.default_rng(0)
        with tempfile.NamedTemporaryFile(suffix=".h5") as f:
            with h5py.File(f.name, "w") as hf:
                grp = hf.create_group("injections")
                grp.create_dataset("m1_source",          data=np.array([m1]))
                grp.create_dataset("m2_source",          data=np.array([m2]))
                grp.create_dataset("a1",                 data=np.zeros(1))
                grp.create_dataset("a2",                 data=np.zeros(1))
                grp.create_dataset("cos_tilt1",          data=np.ones(1))
                grp.create_dataset("cos_tilt2",          data=np.ones(1))
                grp.create_dataset("phi12",              data=np.zeros(1))
                grp.create_dataset("phi_jl",             data=np.zeros(1))
                grp.create_dataset("theta_jn",           data=np.array([0.4]))
                grp.create_dataset("ra",                 data=np.array([1.0]))
                grp.create_dataset("dec",                data=np.array([0.5]))
                grp.create_dataset("psi",                data=np.array([0.3]))
                grp.create_dataset("geocent_time",       data=np.array([gpstime]))
                grp.create_dataset("redshift",           data=np.array([0.1]))
                grp.create_dataset("luminosity_distance", data=np.array([500.0]))
                grp.create_dataset("host_galaxy_index",  data=np.array([0]))
            return read_injection_parameters(f.name)[0]

    # --- structure ----------------------------------------------------------

    def test_returns_dict(self):
        """make_blueprint returns a dict."""
        from minke.bagpuss import make_blueprint
        bp = make_blueprint(self._make_param())
        assert isinstance(bp, dict)

    def test_kind_is_event(self):
        """Blueprint has kind='event'."""
        from minke.bagpuss import make_blueprint
        bp = make_blueprint(self._make_param())
        assert bp["kind"] == "event"

    def test_event_time_matches_gpstime(self):
        """event time equals the injection gpstime."""
        from minke.bagpuss import make_blueprint
        gpstime = 1_187_008_882.0
        bp = make_blueprint(self._make_param(gpstime=gpstime))
        assert math.isclose(bp["event time"], gpstime, rel_tol=1e-10)

    def test_name_set_from_gpstime_by_default(self):
        """Default name encodes the GPS time."""
        from minke.bagpuss import make_blueprint
        bp = make_blueprint(self._make_param(gpstime=1_200_000_000.0))
        assert "name" in bp
        assert "1200000000" in bp["name"]

    def test_custom_name_used(self):
        """Explicit name is used verbatim."""
        from minke.bagpuss import make_blueprint
        bp = make_blueprint(self._make_param(), name="GW_TEST_001")
        assert bp["name"] == "GW_TEST_001"

    # --- chirp mass prior ---------------------------------------------------

    def test_chirp_mass_prior_present(self):
        """priors.chirp mass has minimum and maximum keys."""
        from minke.bagpuss import make_blueprint
        bp = make_blueprint(self._make_param())
        assert "priors" in bp
        assert "chirp mass" in bp["priors"]
        assert "minimum" in bp["priors"]["chirp mass"]
        assert "maximum" in bp["priors"]["chirp mass"]

    def test_chirp_mass_prior_contains_true_value(self):
        """The true chirp mass lies strictly inside the prior range."""
        from minke.bagpuss import make_blueprint
        m1, m2 = 30.0, 20.0
        mc_true = _chirp_mass(m1, m2)
        bp = make_blueprint(self._make_param(m1=m1, m2=m2))
        mc_min = bp["priors"]["chirp mass"]["minimum"]
        mc_max = bp["priors"]["chirp mass"]["maximum"]
        assert mc_min < mc_true < mc_max, (
            f"mc_true={mc_true:.3f} not inside [{mc_min:.3f}, {mc_max:.3f}]"
        )

    def test_chirp_mass_prior_minimum_positive(self):
        """Minimum chirp mass prior is positive."""
        from minke.bagpuss import make_blueprint
        bp = make_blueprint(self._make_param(m1=5.0, m2=3.0))
        assert bp["priors"]["chirp mass"]["minimum"] > 0.0

    def test_chirp_mass_prior_ordered(self):
        """minimum < maximum."""
        from minke.bagpuss import make_blueprint
        bp = make_blueprint(self._make_param())
        assert (
            bp["priors"]["chirp mass"]["minimum"]
            < bp["priors"]["chirp mass"]["maximum"]
        )

    def test_custom_margin_widens_prior(self):
        """A larger margin produces a wider chirp mass prior."""
        from minke.bagpuss import make_blueprint
        p = self._make_param()
        bp_narrow = make_blueprint(p, chirp_mass_margin=0.2)
        bp_wide   = make_blueprint(p, chirp_mass_margin=0.8)
        width_narrow = (bp_narrow["priors"]["chirp mass"]["maximum"]
                        - bp_narrow["priors"]["chirp mass"]["minimum"])
        width_wide   = (bp_wide["priors"]["chirp mass"]["maximum"]
                        - bp_wide["priors"]["chirp mass"]["minimum"])
        assert width_wide > width_narrow


# ===========================================================================
# write_blueprints
# ===========================================================================

class TestWriteBlueprints:
    """Tests for minke.bagpuss.write_blueprints."""

    def test_writes_file(self, injection_hdf5, tmp_path):
        """write_blueprints creates the output file."""
        from minke.bagpuss import read_injection_parameters, write_blueprints
        params = read_injection_parameters(injection_hdf5)
        out = tmp_path / "blueprints.yaml"
        write_blueprints(params, out)
        assert out.exists()

    def test_output_is_valid_yaml(self, injection_hdf5, tmp_path):
        """Output file is valid multi-document YAML."""
        from minke.bagpuss import read_injection_parameters, write_blueprints
        params = read_injection_parameters(injection_hdf5)
        out = tmp_path / "blueprints.yaml"
        write_blueprints(params, out)
        docs = list(yaml.safe_load_all(out.read_text()))
        assert len(docs) == _N

    def test_each_document_is_event_blueprint(self, injection_hdf5, tmp_path):
        """Every YAML document has kind=event."""
        from minke.bagpuss import read_injection_parameters, write_blueprints
        params = read_injection_parameters(injection_hdf5)
        out = tmp_path / "blueprints.yaml"
        write_blueprints(params, out)
        for doc in yaml.safe_load_all(out.read_text()):
            assert doc["kind"] == "event"

    def test_event_times_match_injections(self, injection_hdf5, tmp_path):
        """Event times in the YAML match the injection gpstimes."""
        from minke.bagpuss import read_injection_parameters, write_blueprints
        params = read_injection_parameters(injection_hdf5)
        out = tmp_path / "blueprints.yaml"
        write_blueprints(params, out)
        docs = list(yaml.safe_load_all(out.read_text()))
        for doc, p in zip(docs, params):
            assert math.isclose(doc["event time"], p["gpstime"], rel_tol=1e-10)

    def test_accepts_string_path(self, injection_hdf5, tmp_path):
        """write_blueprints accepts a plain string path."""
        from minke.bagpuss import read_injection_parameters, write_blueprints
        params = read_injection_parameters(injection_hdf5)
        out = str(tmp_path / "blueprints.yaml")
        write_blueprints(params, out)
        assert Path(out).exists()
