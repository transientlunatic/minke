Noise sources
=============

An important part of preparing injections with Minke is being able to add them to a known noise background.

Minke interfaces with the ``lalsimulation`` noise sources to simulate noise for a variety of detectors and design scenarios.

Overview
--------

Minke provides several ways to work with detector noise:

- Generate colored Gaussian noise from power spectral densities (PSDs)
- Use built-in LALSimulation PSDs for various detector configurations
- Work with PSDs in both time and frequency domains
- Create realistic noise backgrounds for injection studies

The noise module is built around the ``LALSimulationPSD`` class, which provides a flexible interface for generating noise colored according to detector sensitivity curves.

Basic Usage
-----------

Generating a Colored Noise Time Series
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The simplest way to generate colored noise is to use one of the built-in detector models:

.. plot::
   :include-source:

   from minke.noise import AdvancedLIGO

   # Create a noise model for Advanced LIGO
   noise = AdvancedLIGO()

   # Generate a time series of colored noise
   # duration: length in seconds
   # sample_rate: sampling frequency in Hz
   data = noise.time_series(duration=4, sample_rate=16384)

   # Plot the noise
   fig = data.plot()
   plt.tight_layout()

The ``time_series`` method generates zero-mean Gaussian noise with a power spectral density matching the Advanced LIGO design sensitivity.

Setting the Time Epoch
~~~~~~~~~~~~~~~~~~~~~~~

You can control the start time of the time series using the ``epoch`` parameter:

::

   from minke.noise import AdvancedLIGO

   noise = AdvancedLIGO()
   
   # Generate noise starting at GPS time 1000
   data = noise.time_series(duration=4, sample_rate=16384, epoch=1000)
   
   print(f"Start time: {data.times[0]}")  # Will be 1000
   print(f"End time: {data.times[-1]}")   # Will be ~1004

Working with PSDs
-----------------

Frequency Domain Representation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The PSD can be accessed in the frequency domain, which is useful for understanding the noise characteristics:

.. plot::
   :include-source:

   from minke.noise import AdvancedLIGO
   import numpy as np

   noise = AdvancedLIGO()

   # Generate PSD at specific frequencies
   frequencies = np.arange(10, 2048, 1)  # 10 to 2048 Hz
   psd = noise.frequency_domain(frequencies=frequencies)

   # Plot the PSD
   fig, ax = plt.subplots(1, 1, figsize=(10, 6))
   ax.loglog(psd.frequencies, np.sqrt(psd.data))
   ax.set_xlabel('Frequency [Hz]')
   ax.set_ylabel('Strain Noise [1/√Hz]')
   ax.set_title('Advanced LIGO Design Sensitivity')
   ax.grid(True, alpha=0.3)
   plt.tight_layout()

Two-Column Format
~~~~~~~~~~~~~~~~~~

For compatibility with other tools, you can export the PSD in two-column ASCII format:

::

   from minke.noise import AdvancedLIGO
   import numpy as np

   noise = AdvancedLIGO()
   
   # Get PSD in two-column format [frequency, strain]
   psd_data = noise.twocolumn(
       lower_frequency=10,
       upper_frequency=2048,
       df=1
   )
   
   # Save to file
   np.savetxt("aLIGO_psd.txt", psd_data)

Covariance Matrix
~~~~~~~~~~~~~~~~~

For some applications, you may need the time-domain covariance matrix:

::

   from minke.noise import AdvancedLIGO
   import numpy as np

   noise = AdvancedLIGO()
   
   # Define time array
   times = np.linspace(0, 1, 1024)  # 1 second at 1024 Hz
   
   # Get covariance matrix
   cov = noise.covariance_matrix(times)
   
   print(f"Covariance matrix shape: {cov.shape}")

Available Noise Models
----------------------

Built-in Models
~~~~~~~~~~~~~~~

Minke currently provides the following built-in noise models:

- ``AdvancedLIGO``: Advanced LIGO design sensitivity (alias for ``AdvancedLIGODesignSensitivity2018``)
- ``AdvancedLIGODesignSensitivity2018``: Advanced LIGO design sensitivity from T1800044

Example usage:

::

   from minke.noise import AdvancedLIGO
   
   # All of these are equivalent
   noise1 = AdvancedLIGO()
   
   from minke.noise import AdvancedLIGODesignSensitivity2018
   noise2 = AdvancedLIGODesignSensitivity2018()

Creating Custom PSD Models
~~~~~~~~~~~~~~~~~~~~~~~~~~~

You can create custom PSD models by subclassing ``LALSimulationPSD`` and providing a LALSimulation PSD function:

::

   from minke.noise import LALSimulationPSD
   import lalsimulation

   class MyCustomPSD(LALSimulationPSD):
       # Point to a LALSimulation PSD function
       psd_function = lalsimulation.SimNoisePSDaLIGODesignSensitivityT1800044

   # Use it
   custom_noise = MyCustomPSD()
   data = custom_noise.time_series(duration=4, sample_rate=4096)

Advanced Features
-----------------

Noise Normalization
~~~~~~~~~~~~~~~~~~~

The noise generation follows standard conventions for one-sided PSDs:

- Interior frequency bins use amplitude ``sqrt(PSD * sample_rate / 2)``
- The DC bin is set to zero (no mean offset)
- The Nyquist bin (when present) is real-valued with amplitude ``sqrt(PSD * sample_rate)``
- The output is normalized to account for the inverse FFT scaling

Multiple Realizations
~~~~~~~~~~~~~~~~~~~~~~

Each call to ``time_series`` generates a new random realization of the noise:

::

   from minke.noise import AdvancedLIGO

   noise = AdvancedLIGO()

   # Generate three independent noise realizations
   noise1 = noise.time_series(duration=4, sample_rate=4096)
   noise2 = noise.time_series(duration=4, sample_rate=4096)
   noise3 = noise.time_series(duration=4, sample_rate=4096)

   # Each will be different (different random draws)

Combining with Signals
~~~~~~~~~~~~~~~~~~~~~~

Noise is commonly used as a background for signal injections:

::

   from minke.noise import AdvancedLIGO
   from minke.models.cbc import IMRPhenomXPHM
   from minke.detector import AdvancedLIGOHanford
   import astropy.units as u

   # Generate noise
   noise = AdvancedLIGO()
   noise_ts = noise.time_series(duration=16, sample_rate=4096, epoch=1000)

   # Generate a signal
   model = IMRPhenomXPHM()
   parameters = {
       "m1": 36 * u.solMass,
       "m2": 29 * u.solMass,
       "luminosity_distance": 400 * u.megaparsec,
       "gpstime": 1008
   }
   waveform = model.time_domain(parameters, times=noise_ts.times)

   # Project onto detector
   detector = AdvancedLIGOHanford()
   signal = waveform.project(detector, ra=1.5, dec=-0.5, psi=0.8, iota=0.4, phi_0=0)

   # Add signal to noise
   injection = noise_ts + signal
   injection.channel = "H1:INJECTION"

   # Save to frame file
   injection.write("H1_injection.gwf", format="gwf")

See Also
--------

- :doc:`tutorial-noise-psd` - Detailed tutorial on using PSDs for injections
- :doc:`software-injections` - Creating injection frames
- :doc:`tutorial-asimov` - Using Minke in Asimov workflows   
