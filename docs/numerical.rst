Notes on using Numerical Relativity Data
========================================

Minke is capable of producing signals which are based on waveforms
calculated by numerical relativity simulations, for example, supernova
and ADI waveforms.

In order to allow the production of these waveforms the standard
LALSuite definition of a SimBurstTable has been altered on the `minke`
branch of LALSuite to contain an extra column to contain a filepath to
a sidecar file containing the raw waveform.

Minke is not distributed with numerical waveform data, however LSC
members may access all of the data required for the production of any
of the supported waveforms by cloning the internal `minke-supernova`
repository, and linking to those files when producing the xml
files. This repository can be cloned by running::
  $ git clone https://git.ligo.org/daniel-williams/minke-supernova.git

If you're using Minke to generate signal sets, but aren't a member of
the LSC then you can download the data for each individual signal set.

Dimmelmeier+08
--------------

The rotational core-collapse simulations detailed in Dimmelmeier et
al. (2008) are available from
http://mpa.iwww.mpg.de/177514/Gravitational-Waveform-Catalog.


Ott+13
------

The 3D hydrodynamic 3-species neutrino heating waveforms from Ott et
al. (2013) are available from https://stellarcollapse.org/ottetal2013.

