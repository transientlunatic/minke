"""Adapter for reading bagpuss injection sets into minke.

Bagpuss (https://github.com/transientlunatic/bagpuss) generates simulated
gravitational-wave injection sets stored as HDF5 files.  Each file contains
one dataset per parameter under the ``/injections`` group, using the bilby
parameter naming convention (``m1_source``, ``cos_tilt1``, ``phi12``, …).

This module bridges the two packages:

1. Reads the HDF5 file produced by ``bagpuss.injection.InjectionSet.to_hdf5``.
2. Converts the bilby-style spin parameterisation to the Cartesian spin
   components that LALSimulation expects, using
   :func:`puddin.lalsim.spins_to_lalsim`.
3. Returns a list of parameter dicts, one per injection, that can be passed
   directly to :func:`minke.injection.make_injection`.

Parameter name mapping
----------------------

+--------------------+---------------------------+-----------------------------------+
| bagpuss name       | minke/LALSim name         | Notes                             |
+====================+===========================+===================================+
| ``m1_source``      | ``m1``                    | With ``astropy.units.solMass``    |
+--------------------+---------------------------+-----------------------------------+
| ``m2_source``      | ``m2``                    | With ``astropy.units.solMass``    |
+--------------------+---------------------------+-----------------------------------+
| ``cos_tilt1/2``    | consumed                  | Used to compute tilt angles       |
+--------------------+---------------------------+-----------------------------------+
| ``a1``, ``a2``     | consumed                  | Used in spin conversion           |
+--------------------+---------------------------+-----------------------------------+
| ``theta_jn``       | consumed → ``iota``       | Transformed by ``spins_to_lalsim``|
+--------------------+---------------------------+-----------------------------------+
| ``phi_jl``         | consumed                  | Used in spin conversion           |
+--------------------+---------------------------+-----------------------------------+
| ``phi12``          | consumed                  | Used in spin conversion           |
+--------------------+---------------------------+-----------------------------------+
| *spin transform*   | ``S1x``, ``S1y``, ``S1z`` | L-frame Cartesian components      |
|                    | ``S2x``, ``S2y``, ``S2z`` |                                   |
+--------------------+---------------------------+-----------------------------------+
| *spin transform*   | ``iota``                  | L-frame inclination (radians)     |
+--------------------+---------------------------+-----------------------------------+
| ``geocent_time``   | ``gpstime``               | Geocentric GPS merger time (s)    |
+--------------------+---------------------------+-----------------------------------+
| ``luminosity_distance`` | ``luminosity_distance`` | With ``astropy.units.Mpc``    |
+--------------------+---------------------------+-----------------------------------+
| ``ra``, ``dec``,   | unchanged                 | Passed through as plain floats    |
| ``psi``, …         |                           |                                   |
+--------------------+---------------------------+-----------------------------------+

Examples
--------
Read an injection set and inject the first event into simulated noise:

.. code-block:: python

    from minke.bagpuss import read_injection_parameters
    from minke.injection import make_injection

    params = read_injection_parameters("injections.h5")

    injections = make_injection(
        injection_parameters=params[0],
        detectors={"AdvancedLIGOHanford": "aLIGOZeroDetHighPower"},
        duration=4,
        sample_rate=4096,
        epoch=params[0]["gpstime"],
    )
"""

from __future__ import annotations

import numpy as np

import math

try:
    import yaml
except ImportError as exc:
    raise ImportError(
        "minke.bagpuss requires PyYAML.  Install it with 'pip install pyyaml'."
    ) from exc

try:
    import h5py
except ImportError as exc:
    raise ImportError(
        "minke.bagpuss requires h5py.  Install it with 'pip install h5py'."
    ) from exc

try:
    from astropy import units as u
except ImportError as exc:
    raise ImportError(
        "minke.bagpuss requires astropy.  Install it with 'pip install astropy'."
    ) from exc

try:
    import lal as _lal
except ImportError as exc:
    raise ImportError(
        "minke.bagpuss requires lal.  "
        "Install it via 'conda install -c conda-forge lalsuite'."
    ) from exc

from puddin import lalsim as _puddin_lalsim

__all__ = ["read_injection_parameters", "make_blueprint", "write_blueprints"]


def read_injection_parameters(
    path: str,
    f_ref: float = 20.0,
    phase: float = 0.0,
) -> list[dict]:
    r"""Read a bagpuss ``InjectionSet`` HDF5 file and return minke-ready dicts.

    Loads all events from *path*, performs the bilby-to-LALSimulation spin
    frame transformation via :func:`puddin.lalsim.spins_to_lalsim`, and
    returns one parameter dict per injection.

    Parameters
    ----------
    path : str or path-like
        Path to an HDF5 file produced by
        ``bagpuss.injection.InjectionSet.to_hdf5``.  The file must contain
        an ``/injections`` group with the standard bagpuss datasets.
    f_ref : float, optional
        Gravitational-wave reference frequency in Hz at which the spin
        components are defined.  Default 20 Hz (typical for O3-era analyses).
    phase : float, optional
        Reference orbital phase in radians.  Default 0.

    Returns
    -------
    list[dict]
        One dict per injection.  Keys are ready for direct use with
        :func:`minke.injection.make_injection`.  See the module docstring for
        the full parameter name mapping.

    Notes
    -----
    The spin transformation calls
    ``lalsimulation.SimInspiralTransformPrecessingNewInitialConditions``
    for precessing spins and uses an analytic fast path for aligned/anti-aligned
    spins; see :func:`puddin.lalsim.spins_to_lalsim` for details.

    The returned ``iota`` is the inclination of the *orbital* angular momentum
    :math:`\mathbf{L}` relative to the line of sight — the quantity used by
    LALSimulation waveform generators.  This differs from the input
    ``theta_jn``, which is the angle between the *total* angular momentum
    :math:`\mathbf{J}` and the line of sight.

    Examples
    --------
    >>> params = read_injection_parameters("injections.h5", f_ref=20.0)
    >>> len(params)          # number of injections in the file
    500
    >>> params[0].keys()
    dict_keys(['m1', 'm2', 'S1x', 'S1y', 'S1z', 'S2x', 'S2y', 'S2z',
               'iota', 'luminosity_distance', 'ra', 'dec', 'psi',
               'gpstime', 'redshift'])
    """
    # ── 1. Load raw arrays from HDF5 ─────────────────────────────────────────
    with h5py.File(path, "r") as f:
        grp = f["injections"]
        raw: dict[str, np.ndarray] = {key: grp[key][()] for key in grp.keys()}

    n = len(raw["m1_source"])

    # ── 2. Convert cos_tilt → tilt angle (radians) ───────────────────────────
    tilt1 = np.arccos(raw["cos_tilt1"])
    tilt2 = np.arccos(raw["cos_tilt2"])

    # ── 3. Masses in SI (kg) for the spin frame transformation ───────────────
    # Use lal.MSUN_SI for consistency with the rest of the minke/LALSim stack.
    m1_kg = raw["m1_source"] * _lal.MSUN_SI
    m2_kg = raw["m2_source"] * _lal.MSUN_SI

    # ── 4. Spin frame transformation: J-frame → L-frame Cartesian ────────────
    iota, s1x, s1y, s1z, s2x, s2y, s2z = _puddin_lalsim.spins_to_lalsim(
        theta_jn=raw["theta_jn"],
        phi_jl=raw["phi_jl"],
        tilt1=tilt1,
        tilt2=tilt2,
        phi12=raw["phi12"],
        a1=raw["a1"],
        a2=raw["a2"],
        m1=m1_kg,
        m2=m2_kg,
        f_ref=np.full(n, float(f_ref)),
        phase=np.full(n, float(phase)),
    )

    # ── 5. Build per-event parameter dicts ───────────────────────────────────
    params = []
    for i in range(n):
        params.append({
            # Masses with units so minke's unit-conversion layer is happy
            "m1": float(raw["m1_source"][i]) * u.solMass,
            "m2": float(raw["m2_source"][i]) * u.solMass,
            # Cartesian spin components in the L-frame (dimensionless)
            "S1x": float(s1x[i]),
            "S1y": float(s1y[i]),
            "S1z": float(s1z[i]),
            "S2x": float(s2x[i]),
            "S2y": float(s2y[i]),
            "S2z": float(s2z[i]),
            # L-frame inclination (replaces theta_jn which was consumed above)
            "iota": float(iota[i]),
            # Distance with units
            "luminosity_distance": float(raw["luminosity_distance"][i]) * u.Mpc,
            # Sky position and polarisation (plain radians)
            "ra":  float(raw["ra"][i]),
            "dec": float(raw["dec"][i]),
            "psi": float(raw["psi"][i]),
            # Merger time: renamed geocent_time → gpstime (LALSim convention)
            "gpstime": float(raw["geocent_time"][i]),
            # Cosmological
            "redshift": float(raw["redshift"][i]),
        })

    return params


def make_blueprint(
    param: dict,
    name: str | None = None,
    chirp_mass_margin: float = 0.5,
    frame_files: dict | None = None,
    network_snr: float | None = None,
) -> dict:
    """Generate a minimal asimov event blueprint from a single injection.

    Produces a ``kind: event`` blueprint containing the GPS event time and a
    broad chirp-mass prior centred on the true injected value.  All other
    fields (data channels, interferometers, likelihood settings) are left for
    the analyst to fill in via asimov's own configuration system.

    Parameters
    ----------
    param : dict
        A single-event parameter dict as returned by
        :func:`read_injection_parameters`.
    name : str or None, optional
        Event name to embed in the blueprint.  Defaults to
        ``inj_{gpstime:.3f}`` if not given.
    chirp_mass_margin : float, optional
        Fractional margin applied symmetrically around the true chirp mass
        to define the prior range:

        .. code-block:: text

            minimum = Mc / (1 + margin)
            maximum = Mc * (1 + margin)

        Default 0.5 gives a range of [Mc/1.5, Mc*1.5] — broad enough to
        encompass typical measurement uncertainty at moderate SNR.
    frame_files : dict or None, optional
        Per-detector frame metadata as returned by
        :func:`~minke.injection.make_injection` (its second return value),
        keyed by detector abbreviation with ``"path"``, ``"channel"`` and
        ``"snr"`` entries.  When given, the blueprint gains
        ``interferometers``, ``data`` (channels and frame paths) and
        ``injected.snr`` fields describing where to find the injection data.
    network_snr : float or None, optional
        Combined optimal SNR across detectors, as returned by
        :func:`~minke.injection.make_injection` (its third return value).
        Recorded under ``injected.network snr`` when *frame_files* is also
        given.

    Returns
    -------
    dict
        Asimov event blueprint ready for serialisation with
        :func:`write_blueprints` or ``yaml.dump``.

    Examples
    --------
    >>> from minke.bagpuss import read_injection_parameters, make_blueprint
    >>> params = read_injection_parameters("injections.h5")
    >>> bp = make_blueprint(params[0])
    >>> bp["kind"]
    'event'
    >>> bp["priors"]["chirp mass"]
    {'minimum': ..., 'maximum': ...}
    """
    m1_msun = param["m1"].value  # astropy Quantity → float in solar masses
    m2_msun = param["m2"].value
    gpstime  = param["gpstime"]

    # Chirp mass in solar masses — cast to plain float for YAML serialisation
    mc = float((m1_msun * m2_msun) ** 0.6 / (m1_msun + m2_msun) ** 0.2)

    event_name = name if name is not None else f"inj_{gpstime:.3f}"

    bp: dict = {
        "kind":       "event",
        "name":       event_name,
        "event time": float(gpstime),
        "priors": {
            "chirp mass": {
                "minimum": round(mc / (1.0 + chirp_mass_margin), 6),
                "maximum": round(mc * (1.0 + chirp_mass_margin), 6),
            },
        },
    }

    if frame_files:
        ifos = sorted(frame_files.keys())
        bp["interferometers"] = ifos
        bp["data"] = {
            "channels":    {ifo: frame_files[ifo]["channel"] for ifo in ifos},
            "data files": {ifo: [frame_files[ifo]["path"]]  for ifo in ifos},
        }
        injected: dict = {
            "snr": {ifo: frame_files[ifo]["snr"] for ifo in ifos if "snr" in frame_files[ifo]},
        }
        if network_snr is not None:
            injected["network snr"] = network_snr
        if injected:
            bp["injected"] = injected

    return bp


def write_blueprints(
    params: list[dict],
    path: str,
    name_prefix: str | None = None,
    chirp_mass_margin: float = 0.5,
    frame_files: list[dict] | None = None,
    network_snrs: list[float] | None = None,
) -> None:
    """Write asimov event blueprints for a list of injections to a YAML file.

    Each injection produces one YAML document separated by ``---``, following
    the multi-document format expected by ``asimov apply``.

    Parameters
    ----------
    params : list[dict]
        List of injection parameter dicts as returned by
        :func:`read_injection_parameters`.
    path : str or path-like
        Output file path.  An existing file is overwritten.
    name_prefix : str or None, optional
        If given, events are named ``{name_prefix}_{i:04d}`` where *i* is the
        zero-based index.  If *None* (default), names are derived from the GPS
        time: ``inj_{gpstime:.3f}``.
    chirp_mass_margin : float, optional
        Passed to :func:`make_blueprint`.  Default 0.5.
    frame_files : list[dict] or None, optional
        Per-injection frame metadata, one entry per element of *params*, each
        as returned by :func:`~minke.injection.make_injection` (its second
        return value).  Passed through to :func:`make_blueprint`.
    network_snrs : list[float] or None, optional
        Per-injection combined optimal SNR, one entry per element of
        *params*, as returned by :func:`~minke.injection.make_injection`
        (its third return value).  Passed through to :func:`make_blueprint`.

    Examples
    --------
    >>> from minke.bagpuss import read_injection_parameters, write_blueprints
    >>> params = read_injection_parameters("injections.h5")
    >>> write_blueprints(params, "blueprints.yaml")
    """
    blueprints = []
    for i, p in enumerate(params):
        name = f"{name_prefix}_{i:04d}" if name_prefix is not None else None
        ff = frame_files[i] if frame_files is not None else None
        snr = network_snrs[i] if network_snrs is not None else None
        blueprints.append(make_blueprint(p, name=name, chirp_mass_margin=chirp_mass_margin, frame_files=ff, network_snr=snr))

    with open(path, "w") as f:
        yaml.safe_dump_all(blueprints, f, default_flow_style=False, sort_keys=False)
