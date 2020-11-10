Signal sets with an SNR cut-off
================================

While Minke's primary purpose is to produce MDC sets for analysis, it
is also capable of producing sets of simulation signals for more
general purposes.

In this tutorial we'll make a set of SineGaussian signals, and impose
an SNR threshold on the waveforms which are included in the set, with
regards to the PSD from the detector.

We'll start by importing minke::
  
  >>> from minke import mdctools, distribution, sources, noise

Then we can create the signal Set. To do that we need to tell Minke which
interferometers are being simulated. We'll assume we're simulating the O1 configuration, which was LIGO Livingston,
and LIGO Hanford (4km).::
  
  >>> mdcset = mdctools.MDCSet(['L1', 'H1'])

Next we define a distribution for the times of the injections. Here we'll produce 1000 signals over the entire O1 run, made randomly according to a uniform distribution over times. ::
  
  >>> times = distribution.uniform_time(tstart=1126620016, tstop=1136995216, number=1000)

We can also define a distribution over strains.::
  
  >>> hrss_values = distribution.log_uniform(5e-23, 1e-20, len(times))

In order to perform thresholding of the injections according to their
SNR we need to define a PSD. We can do that by loading in a file. Here
we'll use the O1 semi-analytic PSD from LIGO-P1200087. ::

  >>> o1psd = noise.PSD("LIGO-P1200087-v18-AdV_EARLY_HIGH.txt")

Making a single injection is simple. All of the injection waveforms
are located in Minke's source module. For our Sine Gaussian we can
make the injection with a single line. ::
  
  >>> sg = sources.SineGaussian(q = 10, frequency=10, hrss=1e-15, time=1126630000, polarisation="linear")

We now calculate the SNRs, both the network and the individual detector SNRs. ::
  
  >>> network, snrs = o1psd.SNR(sg, ['L1', 'H1'])

We can add this injection to our MDC set iff its network SNR > 8, using standard python syntax.::
  
  >>> if network > 8 :
  ...    mdcset + sg

Of course, we need injections for the entire run, so we can set this
up in a Python loop.::
  
  >>> for hrss, time in zip(hrss_values, times):
  ...    sg = sources.WhiteNoiseBurst(q=10, frequency=1000, 
  ...                                  hrss=hrss, time=time)
  ...    network, snrs = o1psd.SNR(sg, ['L1', 'H1'])
  ...    if network > 8 :
  ...       mdcset + sg

Now that we have the MDC set (and it might take a while, especially
for white noise burst sets), we can produce the various data products
that we need for MDC analyses.

The XML file which defines the injection set can be produced using the
save_xml method of the mdcset.::
  
  >>> mdcset.save_xml('signal_set.xml.gz')
