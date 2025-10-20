.. minke documentation master file, created by
   sphinx-quickstart on Tue Jul  9 22:26:36 2013.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


Minke
=====

Minke is a package which aims to make it easier to produce gravitational waveforms in python.
It works with waveform models in the ``lalsimulation`` library maintained by the LIGO-Virgo-Kagra collaborations, and makes them available in a format which is easy to work with while using modern python packages including ``astropy`` and ``gwpy``.

You can use Minke as part of your own code if you need to generate a waveform, or you can use it as part of a larger workflow thanks to its integration with ``asimov``.
It can produce waveform injections for sensitivity and recovery studies, and it can produce large data sets for training machine learning algorithms.

The Minke framework aims to be flexible, and easily extended to include new waveform morphologies, and signal distributions, so if you find something lacking feel free to make a pull request!

User guide
==========

.. toctree::
   :maxdepth: 2
   :caption: User guide
	     
   installation
   quickstart
   other-software

Supported waveforms
===================

.. toctree::
   :maxdepth: 3
   :caption: Waveforms

   burst
   cbc
   ringdowns
   numerical

Noise
=====

.. toctree::
   :maxdepth: 3
   :caption: Noise

   noise
   
Parameter distributions
=======================

.. toctree::
   :maxdepth: 3
   :caption: Parameter distributions

   distributions
   
Data Challenge Construction
===========================
.. toctree::
   :maxdepth: 3
   :caption: Data Challenges

   software-injections
   hardware-injections

   
Tutorials
=========
   
.. toctree::
   :maxdepth: 2
   :caption: Tutorials
	     
   tutorial
   tutorial-editing
   tutorial-SNR-threshold
   tutorial-supernova
   tutorial-hardware
   tutorial-noise-psd
   tutorial-asimov

Developer guide
===============

Minke is under continuous development, and more help is always needed! 

.. toctree::
   :maxdepth: 2
   :caption: Development guide

   contributing
   doctest-examples

Credits and Information
=======================

.. toctree::
   :maxdepth: 2
	   
   authors
   history

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

