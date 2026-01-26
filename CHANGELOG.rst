Minke Changelog
===============

Please accept my apologies for the rattiness of this CHANGELOG; this is an old project and it didn't have the most organised of starts.

2.1.1
=====

This is a bug-fix release and does not introduce any backwards-incompatible changes.

Breaking changes
----------------
This release is not believed to introduce any breaking changes.

Pull requests
-------------
This release contains the following PRs:
+ `github#17 <https://github.com/transientlunatic/minke/pull/17>`_ Fix the interface with htcondor to allow htcondor2 bindings.

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

+ `github!5<https://github.com/transientlunatic/minke/pull/5>`_ Minor bug fixes.

2.0.0
=====

Version 2.0.0 is a major feature version, and represents the start of efforts to modernise the codebase.
We have added initial support for running minke using the asimov automation tool, and some initial support for interaction via a commandline interface.
We have started to refactor the package to work more closely with the astropy and gwpy packages in order to support useful features such as physical units for quantities.
Additionally, we have added support for a wider variety of waveform types than was previously possible in minke, and we now provide initial support for making injections of compact binary (CBC) waveforms.

