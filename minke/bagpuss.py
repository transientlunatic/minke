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

__all__ = ["read_injection_parameters"]


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
