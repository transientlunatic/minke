Notes on using Numerical Relativity Data
========================================

Minke is capable of producing signals which are based on waveforms
calculated by numerical relativity simulations, for example, supernova
and ADI waveforms.

.. note::
   In order to allow the production of these waveforms the standard LALSuite definition of a SimBurstTable has been altered on the `minke` branch of LALSuite to contain an extra column to contain a filepath to a sidecar file containing the raw waveform.

   This is *not* required to make injections directly using python, and is only required if injections should be generated with a SimBurstTable intermediary.

.. note::
      Minke is not distributed with numerical waveform data, however LSC members may access all of the data required for the production of any of the supported waveforms, for example, by cloning the internal `minke-supernova` repository, and linking to those files when producing the xml files. This repository can be cloned by running::
	
  $ git clone https://git.ligo.org/daniel-williams/minke-supernova.git

  If you're using Minke to generate signal sets, but aren't a member of the LSC then you can download the data for each individual signal set.

Encounter waveforms
-------------------

Minke is capable of producing injections using the parabolic encounter waveforms of Bae et al (https://arxiv.org/abs/1701.01548).

For this purpose Minke provides the ``minke.sources.Hyperbolic`` class, which can be used in conjunction with the numerical relativity strain files.
For example:

.. code-block:: python

		source = minke.sources.Hyperbolic(datafile="h_m1_L1.5_l2m2_r300.dat",
                         total_mass=100,
                         extraction=300,
                         distance=0.1,                         
                         time=1126630000,)

Here the extraction radius (``extraction``) must be specified in geometric units, but the ``distance`` parameter should be specified in SI units.
  
Supernova signals
-----------------
  
Dimmelmeier+08
~~~~~~~~~~~~~~

The rotational core-collapse simulations detailed in Dimmelmeier et
al. (2008) are available from
http://mpa.iwww.mpg.de/177514/Gravitational-Waveform-Catalog.


Ott+13
~~~~~~

The 3D hydrodynamic 3-species neutrino heating waveforms from Ott et
al. (2013) are available from https://stellarcollapse.org/ottetal2013.

