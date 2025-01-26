Minke Changelog
===============

Please accept my apologies for the rattiness of this CHANGELOG; this is an old project and it didn't have the most organised of starts.

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

