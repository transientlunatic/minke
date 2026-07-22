=================================
Minke: Making Gravitational Waves
=================================

.. image:: https://zenodo.org/badge/53331163.svg
   :target: https://zenodo.org/badge/latestdoi/53331163

.. image:: https://img.shields.io/pypi/v/minke.svg
   :target: https://pypi.python.org/pypi/minke
   :alt: PyPI version

.. image:: https://github.com/transientlunatic/minke/actions/workflows/docs.yml/badge.svg
   :target: https://github.com/transientlunatic/minke/actions/workflows/docs.yml
   :alt: Documentation build status

.. image:: https://code.daniel-williams.co.uk/minke/_images/minke.png
   :alt: Project Minke Logo


Minke is a Python package to create simulated gravitational-wave signals for a number of different sources, including compact binaries, supernovae, and other burst and ringdown morphologies, and to inject them into detector data or hardware-injection pipelines.

* Free software: ISC license
* Source code: https://github.com/transientlunatic/minke
* Documentation: https://code.daniel-williams.co.uk/minke/

Features
--------

* Produces compact binary coalescence (CBC) waveforms using LALSimulation approximants, selectable by name at runtime
* Produces burst MDCs with Gaussian, SineGaussian, and White Noise Burst ad-hoc waveforms
* Produces ringdown waveforms for MDCs
* Produces numerical relativity burst MDCs for supernovae, including hyperbolic encounters and long-duration searches
* Generates injections directly from `bagpuss <https://github.com/transientlunatic/bagpuss>`_ catalogues
* Specifies injections either by physical parameters or by a target network signal-to-noise ratio (SNR)
* Produces coloured detector noise for a range of known detector PSDs
* Produces GWF frame files, matching frame caches, and GravEn-format log files for MDCs
* Produces hardware-injection ready data files
* Produces SimBurstTable XML files for MDCs
* Integrates with `asimov <https://asimov.docs.ligo.org/asimov/>`_ for automated analysis pipelines, including generating asimov event blueprints from injections
