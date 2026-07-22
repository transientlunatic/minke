Minke Changelog
===============

Please accept my apologies for the rattiness of this CHANGELOG; this is an old project and it didn't have the most organised of starts.

2.2.0
=====

This is a minor feature release which adds support for generating injections from ``bagpuss`` catalogues and automatically producing asimov event blueprints from them, and fixes an intermittent frame-writing failure. It also fixes a bug which caused calculated SNRs to be incorrect.

Major New Features
-------------------

**Bagpuss catalogue adapter**
  A new ``minke.bagpuss`` module reads HDF5 injection sets produced by `bagpuss <https://github.com/transientlunatic/bagpuss>`_, converts the bilby-style spin parameterisation (``cos_tilt1``/``cos_tilt2``, ``a1``/``a2``, ``theta_jn``, ``phi_jl``, ``phi12``) to the Cartesian spin components expected by LALSimulation, and returns parameter dictionaries ready for ``minke.injection.make_injection``.

**Asimov event blueprints from injections**
  ``make_blueprint`` and ``write_blueprints`` generate minimal ``kind: event`` asimov YAML blueprints from a bagpuss injection, including a chirp-mass prior centred on the true injected value and, once available, the injection's frame-file paths and network SNR.

**Runtime-selectable LALSimulation approximants**
  ``LALSimulationApproximant`` gained a ``get_approximant(name)`` factory and an optional ``approximant`` constructor argument, so an approximant can be chosen by name at runtime instead of requiring a dedicated subclass.

Changes
-------

**GWF epoch-precision fix**
  gwpy's ``to_lal()`` computes a timeseries' epoch via a GPS-float to string to GPS round trip, which can lose the last significant digit(s) around GPS ~1.26e9s. When that pushed the series epoch earlier than the frame epoch, LALFrame would intermittently reject the write with "Series start time is earlier than frame start time". Frame writing in ``minke.injection`` now computes the epoch the same way LALFrame does, avoiding the round trip and the intermittent failure.

**SNR calculation fix**
  Corrected the FFT normalisation and frequency-bin handling in ``calculate_network_snr_for_distance`` and ``make_injection`` (``df = sample_rate / N``, one-sided PSD factor of ``2 * df``). SNRs calculated by previous releases were incorrect.

**``make_injection`` now writes frame caches and reports SNR**
  ``make_injection`` writes a matching ``.cache`` file alongside each detector's frame file.

**New PSD**
  Added ``AdvancedLIGOO4Sensitivity`` (``lalsimulation.SimNoisePSDaLIGOAdVO4T1800545``) to ``KNOWN_PSDS`` in ``minke.models.lalnoise``.

**Asimov pipeline fixes**
  Added an ``htcondor2``/``htcondor`` import fallback, switched job submission to use ``schedd.submit``, and fixed the results cache-file glob to look under ``<rundir>/cache``.

**Documentation**
  Added a tutorial walking through generating injections from a bagpuss catalogue (``docs/tutorial-bagpuss.rst``); doctest and tutorial examples updated for gwpy's ``.value`` attribute (previously ``.data``).

**Minor**
  ``minke/sources.py`` now imports ``random`` from ``numpy`` rather than the deprecated ``scipy.random``.

Breaking changes
-----------------

**``make_injection`` return value**
  ``make_injection`` now returns a ``(injections, frame_files, network_snr)`` tuple instead of just ``injections``. Callers relying on the previous single-value return will need to be updated.

Pull requests
-------------
This release contains the following PRs:
+ `github#21 <https://github.com/transientlunatic/minke/pull/21>`_ Bagpuss updates

2.1.2
=====

This is a bug-fix release and does not introduce any backwards-incompatible changes.

Breaking changes
----------------
This release is not believed to introduce any breaking changes.

Pull requests
-------------
This release contains the following PRs:
+ `github#19 <https://github.com/transientlunatic/minke/pull/19>`_ Fix job submission process in the Asimov class to use the ``schedd.submit`` method.


2.1.1
=====

This is a bug-fix release and does not introduce any backwards-incompatible changes.

Breaking changes
----------------
This release is not believed to introduce any breaking changes.

Pull requests
-------------
This release contains the following PRs:
+ `github#17 <https://github.com/transientlunatic/minke/pull/17>`_ Update the htcondor bindings to allow htcondor2 bindings.

2.1.0
=====

This is a minor feature release which introduces significant new functionality for SNR-based injections, updates the noise generation capabilities, and migrates to modern LIGO infrastructure dependencies.

Major New Features
------------------

**SNR-Based Injection Functionality**
  Minke now supports creating injections based on target signal-to-noise ratio (SNR) rather than just physical parameters. This includes functions to calculate network SNR for a given luminosity distance and to find the distance that produces a target network SNR. The ``make_injection`` function has been updated to support SNR-based injection specifications.

**Enhanced Noise Generation**
  The noise generation module has been substantially refactored to improve PSD calculation, support dynamic array library selection (including optional PyTorch support), and provide better control over noise generation parameters. Comprehensive unit tests have been added to ensure reliability.

Changes
-------

**LIGO Infrastructure Migration**
  Updated to use ``igwn-ligolw`` instead of the older ``python-ligo-lw`` package, aligning with current LIGO infrastructure standards. This change is handled automatically by pip during installation.

**Documentation Improvements**
  Enhanced documentation with new doctest examples, expanded noise module documentation with detailed usage instructions, and added tutorials for using Minke with Asimov workflows and generating injections with colored noise.

**Asimov Interface Updates**
  Improved the Asimov interface with better pretty printing support and various bug fixes to enhance integration with the Asimov automation framework.

Breaking Changes
----------------

**Dependency Update**
  The migration from ``python-ligo-lw`` to ``igwn-ligolw`` requires users to have ``igwn-ligolw`` installed in their environment. This is handled automatically by pip when installing or upgrading minke, but users with pinned environments may need to update their dependency specifications.

Pull Requests
-------------

This release contains the following PRs:

+ `github#14 <https://github.com/transientlunatic/minke/pull/14>`_ Asimov fixes
+ `github#13 <https://github.com/transientlunatic/minke/pull/13>`_ Add SNR-based injection functionality and corresponding tests
+ `github#12 <https://github.com/transientlunatic/minke/pull/12>`_ Refactor imports to use igwn_ligolw and update dependencies in pyproject.toml
+ `github#11 <https://github.com/transientlunatic/minke/pull/11>`_ Improve the documentation
+ `github#10 <https://github.com/transientlunatic/minke/pull/10>`_ Update the noise generation and add tests
+ `github#9 <https://github.com/transientlunatic/minke/pull/9>`_ Bump pypa/gh-action-pypi-publish in CI workflow
+ `github#8 <https://github.com/transientlunatic/minke/pull/8>`_ Merge v2-preview branch with SNR calculation
+ `github#7 <https://github.com/transientlunatic/minke/pull/7>`_ Update the asimov interface

2.0.1
=====

This is a bug-fix release and does not introduce any backwards-incompatible changes.

Breaking changes
----------------

This release is not believed to introduce any breaking changes.

Pull requests
-------------

This release contains the following PRs:

+ `github#5 <https://github.com/transientlunatic/minke/pull/5>`_ Minor bug fixes.

2.0.0
=====

Version 2.0.0 is a major feature version, and represents the start of efforts to modernise the codebase.
We have added initial support for running minke using the asimov automation tool, and some initial support for interaction via a commandline interface.
We have started to refactor the package to work more closely with the astropy and gwpy packages in order to support useful features such as physical units for quantities.
Additionally, we have added support for a wider variety of waveform types than was previously possible in minke, and we now provide initial support for making injections of compact binary (CBC) waveforms.

1.1.9 (2021-05-04)
==================

This is a bug-fix release and does not introduce any backwards-incompatible changes.

Changes
-------

+ Corrected the behaviour of numerical relativity waveforms in geometric units by incorporating the extraction radius, and fixed a related bug in the numerical relativity class.
+ Updated the documentation for NR hyperbolic waveforms.

1.1.8 (2020-11-10)
==================

This is a maintenance release focused on the project's build and documentation infrastructure.

Changes
-------

+ Migrated documentation and package building to GitHub Actions.
+ Fixed the destination used when publishing built documentation.
+ Fixed a small bug in the ``Hyperbolic`` waveform class.

1.1.7 (2020-08-24)
==================

Changes
-------

+ Added ``seed`` as an argument to the noise timeseries generator, allowing reproducible noise realisations.
+ Updated the ``Sources`` object to allow selection of a specific epoch.
+ Removed the Python 2.7 wheel build; minke is now Python 3 only.
+ Documentation updates.

1.1.6 (2019-06-22)
==================

This release provides provisional Python 3.6 and 3.7 support.

Changes
-------

+ Began the move to a fully GitLab-CI-oriented build and test workflow.
+ Replaced ad-hoc ``print`` statements with proper logging calls.
+ Made various changes required for Python 3.6/3.7 compatibility.
+ Tests which depend on numerical-relativity data files are now skipped when those files are unavailable.

1.1.5 (2019-06-22)
==================

Changes
-------

+ Added provisional support for numerical relativity waveform files containing hyperbolic encounters.
+ Removed the dependency on ``glue.segments``.

1.1.4 (2019-05-07)
==================

Changes
-------

+ Fixed a bug in frame (GWF) file creation for supernova injections, where the output directory was not always created correctly.

1.1.3 (2019-03-27)
==================

Changes
-------

+ Wheels are now only built and pushed to PyPI when a tag is pushed, rather than on every commit.
+ Minor documentation fixes.

1.1.2 (2019-03-25)
==================

This release contains a number of bug fixes to improve hardware injection support, and continues the removal of the deprecated ``pylal`` dependency.

Changes
-------

+ Corrected the handling of inclination for injected waveforms.
+ Fixed generation of the times list for all loaded XML tables (previously only produced correctly for burst tables).
+ Migrated frame (GWF) production from ``pylal`` to ``lalframe``.
+ Fixed a potential memory leak in supernova waveform tail production.
+ Ensured the same random seed is used to generate a waveform for every detector in a network, so that the same physical signal is injected coherently.
+ Added test coverage for white-noise-burst (WNB) injections.

1.1.1 (2018-03-27)
==================

Changes
-------

+ Updated the package metadata classifiers in ``setup.py`` to allow upload to PyPI.

1.1.0 "Luce Bay" (2018-03-27)
==============================

This is a major feature release, and the first to include a proper automated test suite.

Major New Features
------------------

+ Added full support for hardware injection production, including output of hardware-injection ASCII files.
+ Added the Yakunin supernova waveform family, with accompanying tests.
+ Added ringdown waveform support, including generation of ringdown MDCs from XML tables.
+ Added support for string cusp waveforms.
+ Added support for all lscsoft XML table formats, and support for editing XML files directly.
+ Added experimental support for arbitrary distance-rescaled injections (ADI) and distance rescaling for supernova waveforms, including handling of the memory effect for 3D supernova waveforms.

Changes
-------

+ Added a proper automated test suite to the repository, run via GitLab CI.
+ Migrated XML handling away from ``pylal`` and towards ``glue``, removing further ``pylal`` dependencies.
+ Fixed generation of GraveEn log files so that they comply with LIGO-T040020-01, including a fix to the projected right ascension recorded in the log file.
+ Fixed a bug which caused hardware injections to be all-NaN for ad-hoc waveforms.
+ Added CI-driven documentation builds (GitLab Pages), coverage reporting, and automatic wheel building/PyPI deployment from tagged commits.
+ Added an experimental Singularity container definition, and an experimental Mattermost notification hook.
+ Added the ability to change the size of source plots.

1.0.1 (2017-01-23)
==================

This release adds provisional support for hardware injection production.

Changes
-------

+ Added the ``HWInj`` set object and the ability to output a hardware injection ASCII file.
+ Added the ability to produce multiple ``burst_dist`` values in a single run.
+ Fixed a bug in the read-out of GPS times for ``sim_burst`` tables.

1.0.0 "The Anniversary Update" (2016-09-14)
=============================================

This release brought minke up to readiness for ER10 and O2.

Major New Features
------------------

+ Full support for producing supernova injections.
+ Full support for producing LALSimulation burst waveforms.
+ An extensible interface for specifying parameter distributions.

0.1.0 (2016-03-14)
==================

First release of minke on PyPI.
