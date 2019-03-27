Making Burst MDCs
=================

Minke is currently capable of producing software injections for mock
data challenges for a large number of waveforms, including all of the
waveforms supported by LALSimulation, and also numerical waveforms,
such as supernova waveforms.

Let's get started with a straight-forward example; we'll make an MDC
set of white noise bursts for the O1 run.


We'll start by importing minke:
  >>> from minke import mdctools, distribution, sources

Then we can create the MDC Set. To do that we need to tell Minke which
interferometers took part in the run. For O1 this was LIGO Livingston,
and LIGO Hanford (4km).
  >>> mdcset = mdctools.MDCSet(['L1', 'H1'])

Next, we need to set up the distribution which we'll use to select the
injection times. For the O1 injection set we injected signals at a
fixed cadence throughout the run, but each injection point had a small
amount of "jitter", so the signals could be injected slightly off the
fixed locations.
  >>> times = distribution.even_time(start = 1126620016, stop = 1136995216, rate = 630720, jitter = 20)

We could also have injected them at random times
drawn from a uniform distribution, for example, if we wanted to make 1000 injections over the O1 period.
  >>> times = distribution.uniform_time(start =  1126620016, stop = 1136995216, numer = 1000)

Now that we have a distribution on the times, it's time to make the
MDCs. In O1, however, we wanted a set of injections which had varying
strain amplitudes, or hrss values. So we can define another
distribution for those.
  >>> hrss_values = distribution.log_uniform(5e-23, 1e-20, len(times))

Making a single injection is simple. All of the injection waveforms
are located in Minke's source module. For our white noise burst we can
make the injection with a single line. 
  >>> wnb = sources.WhiteNoiseBurst(duration=0.1, bandwidth=10, frequency=1000,
  ...                               hrss=1e-23, time=1126630000, seed=3)

We can add this injection to our MDC set using standard python syntax.
  >>> mdcset + wnb

Of course, we need injections for the entire run, so we can set this
up in a Python loop.
  >>> for hrss, time in zip(hrss_values, times):
  ...    wnb = sources.WhiteNoiseBurst(duration=0.1, bandwidth=10, frequency=1000, 
  ...                                  hrss=hrss, time=time, seed=3)
  ...    mdcset + wnb

Now that we have the MDC set (and it might take a while, especially
for white noise burst sets), we can produce the various data products
that we need for MDC analyses.

The XML file which defines the injection set can be produced using the
save_xml method of the mdcset.
  >>> mdcset.save_xml('wnb1000b10tau0d1.xml.gz')

We may also want to produce frame files. Minke can do this, but we'll
need to have another piece of information: the start times and
durations of the data frames which we're going to inject the waveforms
into. They should come as an ASCII file of the start times and the
durations, along with the list of interferometers which contributed to
that dataframe. Let's make a frameset.
::
   >>> o1 = mdctools.FrameSet('frame_list.dat')

And now we can make the frames for the MDC set:
::
   >>> mdc_folder = "frames"
   >>> for o1frame in o1.frames:
   ...    o1frame.generate_gwf(mdcset, mdc_folder, 'SCIENCE')

The last argument of the generate_gwf method determines the channel
which the injections are added to in the frame file.

Finally, we can make the GravEn log file:
::
   >>> o1.full_logfile(mdcset, 'frames/logfile.txt')

Let's have a look at the full script:
::
   from minke import mdctools, distribution, sources

   mdcset = mdctools.MDCSet(['L1', 'H1'])

   times = distribution.even_time(start = 1126620016, stop = 1136995216, rate = 630720, jitter = 20)
   hrss_values = distribution.log_uniform(5e-23, 1e-20, len(times))

   for hrss, time in zip(hrss_values, times):
      wnb = sources.WhiteNoiseBurst(duration=0.1, bandwidth=10, frequency=1000,
                                    hrss=hrss, time=time, seed=3)
      mdcset + wnb

   mdcset.save_xml('wnb1000b10tau0d1.xml.gz')

   o1 = mdctools.FrameSet('frame_list.dat')

   mdc_folder = "frames"
   for o1frame in o1.frames:
      o1frame.generate_gwf(mdcset, mdc_folder, 'SCIENCE')

   o1.full_logfile(mdcset, 'frames/logfile.txt')
