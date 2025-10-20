Generating Injections with Colored Noise
=========================================

This tutorial demonstrates how to generate gravitational wave injections into noise that is colored according to a specific power spectral density (PSD). This is useful for realistic simulations that account for detector noise characteristics.

Overview
--------

When creating realistic injection sets, it's important to add signals to noise that matches the characteristics of real detector data. Minke provides tools to:

1. Generate noise colored by detector PSDs
2. Create gravitational wave signals
3. Combine signals with noise
4. Save the results to frame files

Quick Start Examples
--------------------

Here are some quick examples you can try.

Creating a noise model::

   >>> from minke.noise import AdvancedLIGO
   >>> noise_model = AdvancedLIGO()
   >>> noise_model  # doctest: +ELLIPSIS
   <minke.noise.AdvancedLIGO object at 0x...>

Generating colored noise::

   >>> import numpy as np
   >>> from minke.noise import AdvancedLIGO
   >>> noise_model = AdvancedLIGO()
   >>> noise_data = noise_model.time_series(duration=1, sample_rate=1024)
   >>> len(noise_data.data)
   1024
   >>> noise_data.times[0]
   0.0
   >>> noise_data.times[-1]  # doctest: +ELLIPSIS
   0.999...

Getting the PSD in frequency domain::

   >>> from minke.noise import AdvancedLIGO
   >>> import numpy as np
   >>> noise_model = AdvancedLIGO()
   >>> frequencies = np.arange(10, 100, 10)
   >>> psd = noise_model.frequency_domain(frequencies=frequencies)
   >>> len(psd.data)
   9
   >>> psd.frequencies[0]
   10.0
   >>> all(psd.data > 0)  # All PSD values should be positive
   True

Checking noise properties::

   >>> from minke.noise import AdvancedLIGO
   >>> noise_model = AdvancedLIGO()
   >>> noise1 = noise_model.time_series(duration=1, sample_rate=1024)
   >>> noise2 = noise_model.time_series(duration=1, sample_rate=1024)
   >>> # Each realization is different (random)
   >>> import numpy as np
   >>> np.allclose(noise1.data, noise2.data)
   False
   >>> # But both have the same length
   >>> len(noise1.data) == len(noise2.data)
   True

Basic Noise Generation
-----------------------

First, let's generate colored noise using a built-in detector PSD. Minke interfaces with ``lalsimulation`` to provide noise models for various detectors.

.. plot::
   :include-source:

   from minke.noise import AdvancedLIGO

   # Create a noise model using the Advanced LIGO design sensitivity PSD
   noise_model = AdvancedLIGO()

   # Generate a time series of colored noise
   # duration: length of the time series in seconds
   # sample_rate: sampling frequency in Hz
   noise_data = noise_model.time_series(duration=16, sample_rate=4096)

   # Plot the noise time series
   fig = noise_data.plot()
   plt.tight_layout()

The ``time_series`` method generates Gaussian noise with a power spectral density matching the Advanced LIGO design sensitivity curve.

Understanding the PSD
----------------------

You can examine the power spectral density that's being used to color the noise:

.. plot::
   :include-source:

   from minke.noise import AdvancedLIGO
   import numpy as np

   noise_model = AdvancedLIGO()

   # Generate the PSD in the frequency domain
   frequencies = np.arange(10, 2048, 1)  # Hz
   psd = noise_model.frequency_domain(frequencies=frequencies)

   # Plot the PSD
   fig, ax = plt.subplots(1, 1, figsize=(10, 6))
   ax.loglog(psd.frequencies, np.sqrt(psd.data))
   ax.set_xlabel('Frequency [Hz]')
   ax.set_ylabel('Strain noise [1/√Hz]')
   ax.set_title('Advanced LIGO Design Sensitivity')
   ax.grid(True, alpha=0.3)
   plt.tight_layout()

Using a Custom PSD from a File
-------------------------------

In many cases, you may want to use a PSD from actual detector data rather than a theoretical model. While Minke's built-in ``LALSimulationPSD`` class is designed for LALSimulation PSDs, you can extend it to load custom PSDs.

Here's an approach for working with a custom two-column PSD file (frequency, strain):

.. code-block:: python

   import numpy as np
   from minke.noise import LALSimulationPSD
   from minke.types import PSD
   import scipy.interpolate

   class CustomPSD(LALSimulationPSD):
       """A PSD model loaded from a two-column ASCII file."""
       
       def __init__(self, filename):
           super().__init__()
           # Load the PSD data from file
           # Format: frequency [Hz], strain noise [1/√Hz]
           self.psd_data = np.loadtxt(filename)
           self.freq_data = self.psd_data[:, 0]
           self.strain_data = self.psd_data[:, 1]**2  # Convert ASD to PSD
           
       def frequency_domain(self, df=1, frequencies=None, 
                          lower_frequency=20, upper_frequency=1024,
                          mask_below=20):
           """Return the PSD at the requested frequencies."""
           if frequencies is None:
               frequencies = np.arange(lower_frequency, upper_frequency + df, df)
           
           # Interpolate the PSD data to the requested frequencies
           interpolator = scipy.interpolate.interp1d(
               self.freq_data, self.strain_data, 
               kind='linear', fill_value='extrapolate'
           )
           psd_values = interpolator(frequencies)
           
           return PSD(psd_values, frequencies=frequencies)

   # Use the custom PSD
   custom_noise = CustomPSD('my_detector_psd.txt')
   noise_data = custom_noise.time_series(duration=16, sample_rate=4096)

Creating an Injection into Colored Noise
-----------------------------------------

Now let's create a complete example that generates a gravitational wave signal and adds it to colored noise:

.. plot::
   :include-source:

   from minke.noise import AdvancedLIGO
   from minke.models.cbc import IMRPhenomXPHM
   from minke.detector import AdvancedLIGOHanford
   import astropy.units as u

   # Step 1: Generate colored noise
   noise_model = AdvancedLIGO()
   noise_ts = noise_model.time_series(duration=16, sample_rate=4096, epoch=1000)

   # Step 2: Create a binary black hole waveform
   waveform_model = IMRPhenomXPHM()
   
   parameters = {
       "m1": 36 * u.solMass,           # Primary mass
       "m2": 29 * u.solMass,           # Secondary mass
       "luminosity_distance": 400 * u.megaparsec,
       "ra": 1.5,                       # Right ascension
       "dec": -0.5,                     # Declination
       "theta_jn": 0.4,                 # Inclination angle
       "phase": 0,                      # Orbital phase
       "psi": 0.8,                      # Polarization angle
       "gpstime": 1008                  # Merger time (within our noise segment)
   }

   # Generate the waveform
   waveform = waveform_model.time_domain(parameters, times=noise_ts.times)

   # Step 3: Project the waveform onto the detector
   detector = AdvancedLIGOHanford()
   projected_signal = waveform.project(
       detector,
       ra=parameters['ra'],
       dec=parameters['dec'],
       psi=parameters['psi'],
       iota=parameters['theta_jn'],
       phi_0=parameters['phase']
   )

   # Step 4: Add the signal to the noise
   injection = noise_ts + projected_signal
   injection.channel = "H1:INJECTION"

   # Step 5: Visualize
   fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(12, 10))
   
   # Plot noise
   ax1.plot(noise_ts.times, noise_ts.data)
   ax1.set_ylabel('Strain')
   ax1.set_title('Colored Noise')
   ax1.grid(True, alpha=0.3)
   
   # Plot signal
   ax2.plot(projected_signal.times, projected_signal.data)
   ax2.set_ylabel('Strain')
   ax2.set_title('Gravitational Wave Signal')
   ax2.grid(True, alpha=0.3)
   
   # Plot injection (noise + signal)
   ax3.plot(injection.times, injection.data)
   ax3.set_xlabel('Time [s]')
   ax3.set_ylabel('Strain')
   ax3.set_title('Injection (Noise + Signal)')
   ax3.grid(True, alpha=0.3)
   
   plt.tight_layout()

Saving Injections to Frame Files
---------------------------------

Once you've created an injection, you can save it to a LIGO frame file for use in analysis pipelines:

.. code-block:: python

   # Save to a frame file
   injection.write("H1_injection_colored_noise.gwf", format="gwf")

Multiple Detector Injections
-----------------------------

For a realistic multi-detector analysis, you'll want to create injections for multiple detectors:

.. code-block:: python

   from minke.noise import AdvancedLIGO
   from minke.models.cbc import IMRPhenomXPHM
   from minke.detector import AdvancedLIGOHanford, AdvancedLIGOLivingston
   import astropy.units as u

   # Define detectors
   detectors = {
       'H1': AdvancedLIGOHanford(),
       'L1': AdvancedLIGOLivingston()
   }

   # Generate noise for each detector
   noise_model = AdvancedLIGO()
   injections = {}

   # Waveform parameters
   parameters = {
       "m1": 36 * u.solMass,
       "m2": 29 * u.solMass,
       "luminosity_distance": 400 * u.megaparsec,
       "ra": 1.5,
       "dec": -0.5,
       "theta_jn": 0.4,
       "phase": 0,
       "psi": 0.8,
       "gpstime": 1008
   }

   # Generate the waveform once (it's the same for all detectors)
   waveform_model = IMRPhenomXPHM()
   
   for det_name, detector in detectors.items():
       # Generate colored noise
       noise_ts = noise_model.time_series(duration=16, sample_rate=4096, epoch=1000)
       
       # Generate waveform at this detector
       waveform = waveform_model.time_domain(parameters, times=noise_ts.times)
       
       # Project onto this detector
       projected = waveform.project(
           detector,
           ra=parameters['ra'],
           dec=parameters['dec'],
           psi=parameters['psi'],
           iota=parameters['theta_jn'],
           phi_0=parameters['phase']
       )
       
       # Combine with noise
       injection = noise_ts + projected
       injection.channel = f"{det_name}:INJECTION"
       injections[det_name] = injection
       
       # Save to frame file
       injection.write(f"{det_name}_injection.gwf", format="gwf")

   print(f"Created injections for {len(injections)} detectors")

Calculating the Signal-to-Noise Ratio
--------------------------------------

You can calculate the expected signal-to-noise ratio (SNR) of your injection:

.. code-block:: python

   import numpy as np
   from minke.filters import inner_product

   # Generate the signal in frequency domain
   signal_f = np.fft.rfft(projected_signal.data) / sample_rate

   # Get the PSD at the appropriate frequencies
   N = len(noise_ts.data)
   df = sample_rate / N
   frequencies = np.arange(0, N // 2 + 1) * df
   psd = noise_model.frequency_domain(frequencies=frequencies)

   # Calculate optimal SNR using matched filtering
   optimal_snr = np.sqrt(inner_product(signal_f, signal_f, psd.data))
   print(f"Optimal SNR: {optimal_snr:.2f}")

Best Practices
--------------

1. **Choose appropriate duration and sample rate**: The duration should be long enough to contain the entire signal plus some padding. The sample rate should be at least twice the highest frequency of interest (Nyquist criterion).

2. **Match time grids**: Always use the same time grid for the noise and the waveform to ensure proper addition.

3. **Set the epoch correctly**: The ``epoch`` parameter sets the start time of your data. Make sure your signal's GPS time falls within the noise time series.

4. **PSD frequency range**: Ensure your PSD covers the frequency range of your signal. For binary black holes, this is typically 10-2048 Hz for Advanced LIGO.

5. **Random seeds**: When generating multiple noise realizations, you may want to set different random seeds to ensure independent noise samples.

Summary
-------

This tutorial covered:

- Generating noise colored by detector PSDs
- Understanding and visualizing PSDs
- Loading custom PSDs from files
- Creating gravitational wave injections into colored noise
- Working with multiple detectors
- Saving injections to frame files
- Calculating expected SNR values

These techniques are essential for creating realistic injection sets for detector characterization, algorithm development, and validation studies.
