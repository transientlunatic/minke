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




At the command line::

    $ easy_install minki

Or, if you have virtualenvwrapper installed::

    $ mkvirtualenv minki
    $ pip install minki
