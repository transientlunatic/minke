Making Software MDCs for Supernova Searches
===========================================

Minke is currently capable of producing software injections for mock
data challenges for a large number of waveforms, including all of the
waveforms supported by LALSimulation, and also numerical waveforms,
such as supernova waveforms.

In this tutorial we will produce a set of supernova injections using a
set of pre-rotated waveforms.

We'll start by importing minke::
  >>> from minke import mdctools, distribution, sources

Then we can create the MDC Set. To do that we need to tell Minke which
interferometers took part in the run. For O1 this was LIGO Livingston,
and LIGO Hanford (4km).::
  >>> mdcset = mdctools.MDCSet(['L1', 'H1'])

Next, we need to set up the distribution which we'll use to select the
injection times. For the O1 injection set we inject signals at a fixed
cadence throughout the run, but each injection point can have a small
amount of "jitter", so the signals could be injected slightly off the
fixed locations. This is simply a random offset up to 20 seconds
either side of the fixed location, in the case that `jitter = 20`.
Note that the start and end times are the gps times which bounded the
observing run. ::
  >>> times = distribution.even_time(start = 1126620016, stop = 1136995216, rate = 630720, jitter = 20)

We could also have injected them at random times drawn from a uniform
distribution, for example, if we wanted to make 1000 injections over
the O1 period.::
  >>> times = distribution.uniform_time(start =  1126620016, stop = 1136995216, number = 1000)

We also need a distribution of theta and phi values for the
orientation of the supernova. Since we're using pre-computed rotations
we need to use the `supernova_angle` distribution, which by default
assumes that there are 100 evenly spaced combinations. ::
  >>> angles = distribution.supernova_angle(len(times))

Making a single injection is simple. All of the injection waveforms
are located in Minke's `source` module. We can create an Ott+13
waveform, for example, with the line ::
  >>> sn = sources.Ott2013(theta = 0 , phi = 0, time=1126620016, 
  ...                      filepath="/home/daniel.williams/repositories/snsearch/RawWaveforms/", 
  ...			   family="R1E1CA")

We can add this injection to our MDC set using standard python syntax.::
  >>> mdcset + wnb

Of course, we need injections for the entire run, so we can set this
up in a Python loop.::
  >>> for time, (theta, phi) in zip(times, angles):
  ...    sn = sources.Ott2013(theta = 0 , phi = 0, time=1126620016, 
  ...                         filepath="/home/daniel.williams/repositories/snsearch/RawWaveforms/", 
  ...	   		      family="R1E1CA")             
  ...    mdcset + sn

Now that we have the MDC set (and it might take a while, especially
for white noise burst sets), we can produce the various data products
that we need for MDC analyses.

The XML file which defines the injection set can be produced using the
save_xml method of the mdcset.::
  >>> mdcset.save_xml('ott13.xml.gz')
