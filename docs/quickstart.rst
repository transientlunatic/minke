A Quick Start Guide for Minke
=============================

This tutorial will quickly walk you through the steps required to produce a waveform using one of Minke's supported waveform generators, and is a good place to start if you're encountering the package for the first time.

First you'll need to import minke, you can do this by running ::

  import minke

We'll also need to have access to the ``astropy.units`` package, which we'll import as ``u``. ::

  import astropy.units as u

This will allow us to explicitly state the units which we're using to retrieve the waveform, for example ``10 * u.solMass`` indicates a quantity of 10 solar masses.

Next we'll import and set up our waveform generator ::

  from minke.models.cbc import IMRPhenomPv2
  model = IMRPhenomPv2()

The waveform generator will allow us to produce representations of the waveform in both the time and frequency domains.
At the moment let's concentrate on producing a waveform in the time domain, for which we'll need to specify a few parameters. ::

  parameters = {"mass_ratio": 0.7, "total_mass": 100*u.solMass, "luminosity_distance": 10*u.megaparsec}

Then to generate the waveform we use the `time_domain` method on the generator ::

  data = model.time_domain(parameters)

Alternatively you can specify the parameters as keyword arguments ::

  data = model.time_domain(mass_ratio=0.7, total_mass=100*u.solMass, luminosity_distance: 10*u.megaparsec)

Either of these will return a ``WaveformDict`` which contains the two polarisations of the waveform, ``plus`` and ``cross``.

You can plot an individual polarisation ::

  f = data['plus'].plot()
  f.savefig("waveform.png")

The waveforms which Minke returns use ``gwpy`` timeseries to store the waveform, so it is possible to produce spectrograms, and various other useful plots which gwpy supports.
For example ::

  specgram = data['plus'].spectrogram(20, fftlength=8, overlap=4) ** (1/2.)
  specgram.plot(norm='log', vmin=1e-23, vmax=1e-19)


We can also produce burst-like waveforms using the ``minke.models.bursts`` module:

For example we can produce a sine-Gaussian with ::

  from minke.models.bursts import SineGaussian
  
  model = SineGaussian()

  parameters = {"centre_frequency": 20,
                "phase": 0,
                "eccentricity": 0,
                "q": 1.,
                "hrss": 1e-22,
                "duration": 2}
		
  data = model.time_domain(parameters)
  f = data['plus'].plot()


Normally we'll want to be able to produce the waveform as it is received by a detector.
For this we need to project the waveform onto the detector.
For example, if we want to project in for LIGO Hanford we need to specify the detector, the time, sky location, polarisation angle, the initial phase, and the inclination angle of the source. ::

  from minke.detector import AdvancedLIGOHanford
  
  detector = AdvancedLIGOHanford()

  projected = data.project(detector,
               time=1000,
	       ra=1, dec=0.5,
	       iota=0.4,
	       phi_0=0,
	       psi=0
	       )

  f = projected.plot()
  f.savefig("projected_waveform.png")

Creating Injections with Target SNR
------------------------------------

When creating injections for testing parameter estimation pipelines, you may want to specify a target signal-to-noise ratio (SNR) rather than a luminosity distance.
Minke can automatically calculate the required luminosity distance to achieve a specific network SNR.

To create an injection with a target SNR of 20 across multiple detectors::

  from minke.injection import make_injection

  detectors = {"AdvancedLIGOHanford": "AdvancedLIGO", 
               "AdvancedLIGOLivingston": "AdvancedLIGO"}

  parameters = {"m1": 30,
                "m2": 30,
                "snr": 20,  # Target network SNR instead of luminosity_distance
                "ra": 0,
                "dec": 0,
                "psi": 0,
                "theta_jn": 0,
                "phase": 0}

  injections = make_injection(detectors=detectors, 
                              injection_parameters=parameters, 
                              duration=4, 
                              sample_rate=4096, 
                              epoch=998)

The ``make_injection`` function will determine the luminosity distance needed to produce a network SNR of 20.
The network SNR is the quadrature sum of the individual detector SNRs, so each detector will have an SNR that depends on its antenna pattern and the source sky location.

You can access individual detector injections::

  h1_injection = injections['H1']
  l1_injection = injections['L1']

  # Plot the injection
  f = h1_injection.plot()
  f.savefig("h1_snr_injection.png")


