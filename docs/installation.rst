.. highlight:: shell

============
Installation
============

Installation of the Minke package is a little more involved than the
majority of Python packages, due to its dependence upon LALSuite, and,
unforunately, a development version of LALSuite.

We've included a script to install the necessary LAL software into a
virtual environment alongside the Minke code, which should make your
life slightly easier, and in the future we hope to find a way to
bundle this installation into the Minke install script.

To install lalsuite and minke you should run::
  $ ./setup_environment.sh
in the root of the project, which will attempt to download and build the 
latest version of lalsuite, using the correct branch to run the minke code.
You can then install the Minke python package by running::
  $ python setup.py install

=======================
Use on the LDG Clusters
=======================

Daniel Williams currently maintains a virtual environment on the LLO
cluster of the LIGO Data Grid for general use of Minke. You can
activate the virtual environment, which gives access to both the
correct branch of LALSuite and to the Minke Python package, by running
the following command in BASH::
  $ source ~daniel.williams/.virtualenvs/lal-minke/bin/activate

Any code which is run in this virtual environment should link
automatically to the correct python packages to produce signals using
Minke.
