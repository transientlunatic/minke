import numpy as np
import pytest
from scipy.signal import welch

from minke.models.lalnoise import KNOWN_PSDS
from minke.types import Waveform, WaveformDict
from minke.detector import KNOWN_IFOS


@pytest.mark.parametrize("name", ["AdvancedLIGO", "AdvancedLIGO_O4"])
@pytest.mark.parametrize("duration", [32, 128])
def test_noise_matches_psd(name, duration):
    """Generated noise must reproduce its PSD, independent of segment length."""
    fs = 2048
    model = KNOWN_PSDS[name]()
    np.random.seed(1)
    data = np.asarray(model.time_series(duration=duration, sample_rate=fs).data)
    f, p = welch(data, fs=fs, nperseg=4 * fs, average="median")
    ref = np.asarray(
        model.frequency_domain(frequencies=f[f >= 25], lower_frequency=f[f >= 25][0]).data
    )
    ratio = np.median(p[f >= 25] / ref)
    # median averaging of chi-squared bins is biased low by ~ln2/... ~ 6%
    assert ratio == pytest.approx(1.0, abs=0.12)


def test_project_skips_inclination_factors_for_inclined_waveforms():
    """Waveforms already evaluated at an inclination get only the antenna response."""
    n, dt = 64, 1 / 4096.0
    hp = Waveform(data=np.ones(n), dt=dt, t0=1000.0)
    hx = Waveform(data=np.zeros(n), dt=dt, t0=1000.0)
    params = {"ra": 1.0, "dec": 0.3, "psi": 0.4, "phase": 0.0, "theta_jn": 0.0}
    det = KNOWN_IFOS["AdvancedLIGOHanford"]()
    plain = WaveformDict(parameters=dict(params), plus=hp, cross=hx).project(det)
    inclined = WaveformDict(parameters={**params, "inclination": 1.0}, plus=hp, cross=hx).project(det)
    # theta_jn=0 gives the (1 + cos^2) = 2 prefactor; the inclined path gives 1.
    assert np.allclose(np.asarray(plain.data), 2 * np.asarray(inclined.data))


def test_project_inclination_only_and_phase():
    """`inclination` alone suffices, and `phase` still mixes the polarisations."""
    n, dt = 64, 1 / 4096.0
    hp = Waveform(data=np.ones(n), dt=dt, t0=1000.0)
    hx = Waveform(data=np.zeros(n), dt=dt, t0=1000.0)
    params = {"ra": 1.0, "dec": 0.3, "psi": 0.4, "phase": 0.0, "inclination": 1.0}
    det = KNOWN_IFOS["AdvancedLIGOHanford"]()
    zero = np.asarray(WaveformDict(parameters=dict(params), plus=hp, cross=hx).project(det).data)
    rotated = np.asarray(
        WaveformDict(parameters={**params, "phase": 0.7}, plus=hp, cross=hx).project(det).data
    )
    assert not np.allclose(zero, rotated)


def test_noise_odd_length_keeps_last_bin():
    """Odd-length series have no Nyquist bin, so the top bin must not be zeroed."""
    fs, N = 2048, 4097
    model = KNOWN_PSDS["AdvancedLIGO"]()
    np.random.seed(1)
    data = np.asarray(model.time_series(times=np.arange(N) / fs).data)
    assert len(data) == N and np.all(np.isfinite(data))
    assert abs(np.fft.rfft(data)[-1]) > 0
