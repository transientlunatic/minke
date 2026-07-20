Injections from a bagpuss catalogue
=====================================

`Bagpuss <https://github.com/transientlunatic/bagpuss>`_ is a companion
package that simulates galaxy catalogues and binary black hole (BBH)
populations for cosmological inference testing.  It produces
``InjectionSet`` HDF5 files containing fully-specified events — masses,
spins, sky positions, redshifts, and merger times.

This tutorial shows how to take a bagpuss injection set, turn each
event into a real detector strain frame using Minke, and generate
`Asimov <https://asimov.docs.ligo.org/>`_ event blueprints so that
parameter estimation can be launched automatically once the frames are
ready.

Prerequisites
-------------

Both packages must be installed, along with ``puddin`` (for the spin
frame transformation) and LALSuite:

.. code-block:: bash

    pip install bagpuss puddin
    conda install -c conda-forge lalsuite

Generating a bagpuss injection set
------------------------------------

If you do not already have an HDF5 file, the following snippet generates
one from scratch.  See the `bagpuss documentation
<https://bagpuss.readthedocs.io>`_ for a full explanation of each stage.

.. code-block:: python

    import numpy as np
    from astropy.cosmology import Planck18

    from bagpuss.catalogue import MagnitudeLimitedSurvey
    from bagpuss.injection import DistanceThreshold, create_injection_set
    from bagpuss.luminosity import SchechterLuminosityModel
    from bagpuss.population import (
        IsotropicSpinDistribution,
        PopulationModel,
        PowerLawPlusPeakMassDistribution,
    )
    from bagpuss.universe import PointProcess, Universe

    rng = np.random.default_rng(42)

    # Stages 1–3: simulated galaxy catalogue
    universe = Universe(
        cosmology=Planck18,
        structure=PointProcess(z_max=0.3),
        luminosity=SchechterLuminosityModel(
            phi_star=1.61e-2, m_star=-19.66, alpha=-1.16,
            m_min=-25.0, m_max=-14.0,
        ),
    )
    galaxies = universe.sample(5_000, rng=rng)
    catalogue = MagnitudeLimitedSurvey(m_lim=19.5).apply(galaxies, Planck18)

    # Stage 4: BBH population (O3 best-fit hyperparameters)
    population = PopulationModel(
        mass=PowerLawPlusPeakMassDistribution(
            alpha=3.5, beta_q=1.4, m_min=5.0, m_max=87.0,
            lambda_peak=0.03, mu_m=34.0, sigma_m=3.6, delta_m=4.8,
        ),
        spin=IsotropicSpinDistribution(),
    )

    # Stage 5: injection set, keeping only events within 800 Mpc
    injections = create_injection_set(
        catalogue=catalogue,
        population=population,
        cosmology=Planck18,
        n_draw=200,
        detectable=DistanceThreshold(d_max=800.0),
        rng=rng,
    )

    injections.to_hdf5("injections.h5")
    print(f"Saved {len(injections)} injections")

Reading the injection set into Minke
--------------------------------------

:func:`minke.bagpuss.read_injection_parameters` loads the HDF5 file and
returns a list of parameter dicts — one per event — that are ready for
:func:`minke.injection.make_injection`.

The function handles the spin parameterisation conversion automatically:
bagpuss stores spins as tilt angles and magnitudes (bilby convention),
while LALSimulation expects Cartesian components in the orbital frame.
The conversion is performed via
:func:`puddin.lalsim.spins_to_lalsim`.

.. code-block:: python

    from minke.bagpuss import read_injection_parameters

    params = read_injection_parameters("injections.h5", f_ref=20.0)

    print(f"Loaded {len(params)} injections")
    print("First event keys:", list(params[0].keys()))
    # ['m1', 'm2', 'S1x', 'S1y', 'S1z', 'S2x', 'S2y', 'S2z',
    #  'iota', 'luminosity_distance', 'ra', 'dec', 'psi', 'gpstime', 'redshift']

The ``f_ref`` argument sets the gravitational-wave reference frequency
(in Hz) at which the spin components are defined.  20 Hz is the standard
choice for O3-era analyses.

Injecting a single event
-------------------------

Pass any element of the list directly to :func:`~minke.injection.make_injection`:

.. code-block:: python

    from minke.injection import make_injection

    detectors = {
        "AdvancedLIGOHanford":    "aLIGOZeroDetHighPower",
        "AdvancedLIGOLivingston": "aLIGOZeroDetHighPower",
    }

    event = params[0]

    injections, frame_files, network_snr = make_injection(
        injection_parameters=event,
        detectors=detectors,
        duration=8,
        sample_rate=4096,
        epoch=event["gpstime"] - 4,   # 4 s before merger
        framefile="injection",
    )

    # injections["H1"] and injections["L1"] are GWPy TimeSeries objects.
    # The .gwf files are written alongside a cache file automatically.
    # frame_files maps each detector to its frame path, channel name and
    # optimal SNR; network_snr is the combined optimal SNR across detectors.

This produces ``H1_injection.gwf`` and ``L1_injection.gwf`` containing
the signal injected into coloured Gaussian noise.

Using a different waveform approximant
---------------------------------------

By default :func:`~minke.injection.make_injection` uses
:class:`~minke.models.lalsimulation.IMRPhenomXPHM`.  Any LALSimulation
approximant can be selected by name using
:func:`~minke.models.lalsimulation.get_approximant`:

.. code-block:: python

    from minke.injection import make_injection
    from minke.models.lalsimulation import get_approximant

    injections, frame_files, network_snr = make_injection(
        waveform=get_approximant("SEOBNRv4PHM"),
        injection_parameters=params[0],
        detectors={"AdvancedLIGOHanford": "aLIGOZeroDetHighPower"},
        duration=8,
        sample_rate=4096,
        epoch=params[0]["gpstime"] - 4,
    )

Any approximant name accepted by
``lalsimulation.GetApproximantFromString`` is valid — including
``NRSur7dq4``, ``IMRPhenomXO4a``, ``SEOBNRv5PHM``, and so on.

Generating asimov blueprints
-----------------------------

Once the frames exist, PE can be launched via
`Asimov <https://asimov.docs.ligo.org/>`_.  Asimov reads *event blueprints*
— YAML documents describing the event time and analysis priors — and uses
them to configure and submit PE jobs automatically.

:func:`~minke.bagpuss.write_blueprints` produces a multi-document YAML file
from the injection parameter list, with one ``kind: event`` document per
injection.  Each document contains:

* ``event time`` — the geocentric GPS merger time taken directly from the injection.
* ``priors.chirp mass`` — a broad prior range centred on the true injected
  chirp mass.  The default margin is ±50 %, i.e. the range
  ``[Mc / 1.5, Mc × 1.5]``, which is wide enough to accommodate typical
  measurement uncertainty without wasting sampler iterations.

.. code-block:: python

    from minke.bagpuss import read_injection_parameters, write_blueprints

    params = read_injection_parameters("injections.h5")
    write_blueprints(params, "blueprints.yaml")

The resulting file looks like this (two events shown):

.. code-block:: yaml

    kind: event
    name: inj_1187008882.000
    event time: 1187008882.0
    priors:
      chirp mass:
        minimum: 18.145807
        maximum: 40.828065
    ---
    kind: event
    name: inj_1187009000.000
    event time: 1187009000.0
    priors:
      chirp mass:
        minimum: 5.663247
        maximum: 12.742305

Each document can be passed directly to ``asimov apply``:

.. code-block:: bash

    asimov apply --file blueprints.yaml

Asimov will create one event per document and queue them for PE once the
corresponding injection frames are available.

Customising event names and prior width
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

By default each event is named after its GPS time (``inj_<gpstime>``).
Pass ``name_prefix`` to get sequentially numbered names instead, which is
cleaner for large catalogues:

.. code-block:: python

    write_blueprints(params, "blueprints.yaml", name_prefix="bbh_study")
    # → bbh_study_0000, bbh_study_0001, …

To tighten or loosen the chirp mass prior, adjust ``chirp_mass_margin``.
A value of ``0.3`` gives a [Mc/1.3, Mc×1.3] range; ``1.0`` gives [Mc/2,
Mc×2]:

.. code-block:: python

    write_blueprints(params, "blueprints.yaml", chirp_mass_margin=0.3)

For a single event, :func:`~minke.bagpuss.make_blueprint` returns the
blueprint dict directly without writing to disk:

.. code-block:: python

    from minke.bagpuss import make_blueprint

    bp = make_blueprint(params[0], name="GW_test_001")
    print(bp["priors"]["chirp mass"])
    # {'minimum': 18.145807, 'maximum': 40.828065}

Recording frame locations and SNR in the blueprint
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you have already generated frames with
:func:`~minke.injection.make_injection`, pass its ``frame_files`` and
``network_snr`` outputs through to :func:`~minke.bagpuss.make_blueprint` (or
:func:`~minke.bagpuss.write_blueprints`) so that Asimov knows where to find
the data and does not need to search for it itself:

.. code-block:: python

    injections, frame_files, network_snr = make_injection(
        injection_parameters=event,
        detectors=detectors,
        duration=8,
        sample_rate=4096,
        epoch=event["gpstime"] - 4,
        framefile="injection",
    )

    bp = make_blueprint(event, frame_files=frame_files, network_snr=network_snr)

This adds ``interferometers``, ``data.channels``, ``data.data files`` and
``injected.snr`` / ``injected.network snr`` fields to the blueprint,
pointing directly at the frame(s) just written.

Full end-to-end example
------------------------

The following snippet runs the complete pipeline: generate injection frames
for every event, then write the asimov blueprints — including the frame
locations and SNRs — ready for PE submission.

.. code-block:: python

    import os
    from minke.bagpuss import read_injection_parameters, write_blueprints
    from minke.injection import make_injection

    detectors = {
        "AdvancedLIGOHanford":    "aLIGOZeroDetHighPower",
        "AdvancedLIGOLivingston": "aLIGOZeroDetHighPower",
    }

    params = read_injection_parameters("injections.h5")

    os.makedirs("frames", exist_ok=True)

    all_frame_files = []
    all_network_snrs = []
    for i, event in enumerate(params):
        _, frame_files, network_snr = make_injection(
            injection_parameters=event,
            detectors=detectors,
            duration=8,
            sample_rate=4096,
            epoch=event["gpstime"] - 4,
            framefile=f"frames/event_{i:04d}",
        )
        all_frame_files.append(frame_files)
        all_network_snrs.append(network_snr)

    # Write blueprints for all events in one go, including frame locations and SNRs
    write_blueprints(
        params, "blueprints.yaml", name_prefix="bbh_study",
        frame_files=all_frame_files, network_snrs=all_network_snrs,
    )

    print(f"Injected {len(params)} events.")
    print("Run 'asimov apply --file blueprints.yaml' to queue PE.")

Parameter reference
--------------------

The table below lists every key in the dicts returned by
:func:`~minke.bagpuss.read_injection_parameters`.

+-------------------------+-------------------------------------+----------------------------------+
| Key                     | Type / units                        | Description                      |
+=========================+=====================================+==================================+
| ``m1``, ``m2``          | ``astropy.Quantity`` (``solMass``)  | Source-frame component masses    |
+-------------------------+-------------------------------------+----------------------------------+
| ``S1x``, ``S1y``,       | ``float`` (dimensionless)           | L-frame Cartesian spin           |
| ``S1z``, ``S2x``,       |                                     | components; \|S\| = *a*          |
| ``S2y``, ``S2z``        |                                     |                                  |
+-------------------------+-------------------------------------+----------------------------------+
| ``iota``                | ``float`` (radians)                 | Inclination of **L** to l.o.s.   |
+-------------------------+-------------------------------------+----------------------------------+
| ``luminosity_distance`` | ``astropy.Quantity`` (``Mpc``)      | Luminosity distance              |
+-------------------------+-------------------------------------+----------------------------------+
| ``ra``, ``dec``         | ``float`` (radians)                 | Sky position of host galaxy      |
+-------------------------+-------------------------------------+----------------------------------+
| ``psi``                 | ``float`` (radians)                 | Polarisation angle               |
+-------------------------+-------------------------------------+----------------------------------+
| ``gpstime``             | ``float`` (GPS seconds)             | Geocentric merger time           |
+-------------------------+-------------------------------------+----------------------------------+
| ``redshift``            | ``float`` (dimensionless)           | Host galaxy redshift             |
+-------------------------+-------------------------------------+----------------------------------+

.. note::

    The spin components ``S1x`` … ``S2z`` and ``iota`` are derived from the
    bagpuss parameters ``a1``, ``a2``, ``cos_tilt1``, ``cos_tilt2``,
    ``phi12``, ``phi_jl``, and ``theta_jn`` via the LALSimulation frame
    transformation.  The original bilby-convention parameters are consumed
    and do not appear in the output dicts.

API reference
-------------

.. autofunction:: minke.bagpuss.read_injection_parameters

.. autofunction:: minke.bagpuss.make_blueprint

.. autofunction:: minke.bagpuss.write_blueprints
