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

Installing in a virtual environment
===================================

We need to setup a virtual environment first. If you're using
virtualenvwrapper this is straightforward ::
  
  $ mkvirtualenv <name for venv>

First install lalsuite from the minke branch on git ::
  
  $ git clone https://github.com/lscsoft/lalsuite.git
  $ git checkout minke
  $ ./00boot
  $ configure --prefix=$VIRTUALENV
  $ make && make install

Then clone the minke repository, and install the code and dependencies::
  
  $ git clone https://daniel-williams@git.ligo.org/daniel-williams/minke.git
  $ pip install numpy scipy matplotlib pandas
  $ cd minke
  $ python setup.py install

  
System-wide installation
========================

We do not recommend system-wide installation of minke, however, it is possible.

First install lalsuite from the minke branch on git ::
  
  $ git clone https://github.com/lscsoft/lalsuite.git
  $ git checkout minke
  $ ./00boot
  $ configure
  $ make
  $ sudo make install

Then clone the minke repository, and install the code ::
  
  $ git clone https://daniel-williams@git.ligo.org/daniel-williams/minke.git
  $ sudo pip install numpy scipy matplotlib pandas
  $ cd minke
  $ sudo python setup.py install


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


Using Minke in docker
=====================

Minke is now available in a docker container, allowing you to run
Minke on any system which has docker installed. The container is
available on both dockerhub and git.ligo.org's registry.

To get the minke image from dockerhub just run::

  docker pull lpmn/minke

The process for fetching the image from gitlab is a little more
involved, as you need to be logged in to the registry first. ::

  docker login git.ligo.org
  docker pull git.ligo.org/registry/daniel.williams/minke


