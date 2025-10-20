Doctest Examples
================

This document contains testable code examples for the Minke documentation.
Run ``make doctest`` in the docs directory to verify these examples.

Noise Module Tests
------------------

Basic Noise Generation
~~~~~~~~~~~~~~~~~~~~~~~

Test creating a noise model::

   >>> from minke.noise import AdvancedLIGO
   >>> noise = AdvancedLIGO()
   >>> type(noise).__name__
   'AdvancedLIGO'

Test time series generation with different durations::

   >>> noise = AdvancedLIGO()
   >>> ts_1sec = noise.time_series(duration=1, sample_rate=1024)
   >>> ts_2sec = noise.time_series(duration=2, sample_rate=1024)
   >>> len(ts_1sec.data)
   1024
   >>> len(ts_2sec.data)
   2048

Test time series with epoch::

   >>> noise = AdvancedLIGO()
   >>> ts = noise.time_series(duration=1, sample_rate=1024, epoch=1000)
   >>> ts.times[0]
   1000.0
   >>> ts.times[1]  # doctest: +ELLIPSIS
   1000.000...

PSD Generation
~~~~~~~~~~~~~~

Test frequency domain PSD::

   >>> import numpy as np
   >>> from minke.noise import AdvancedLIGO
   >>> noise = AdvancedLIGO()
   >>> frequencies = np.array([10.0, 20.0, 50.0, 100.0])
   >>> psd = noise.frequency_domain(frequencies=frequencies)
   >>> len(psd.data) == len(frequencies)
   True
   >>> all(psd.frequencies == frequencies)
   True

Test PSD with frequency range::

   >>> noise = AdvancedLIGO()
   >>> psd = noise.frequency_domain(lower_frequency=20, upper_frequency=100, df=10)
   >>> len(psd.data)
   9
   >>> float(psd.frequencies[0])
   20.0
   >>> float(psd.frequencies[-1])
   100.0

Test that PSD values are positive::

   >>> import numpy as np
   >>> noise = AdvancedLIGO()
   >>> psd = noise.frequency_domain(lower_frequency=10, upper_frequency=200, df=1)
   >>> all(psd.data > 0)
   True
   >>> np.isfinite(psd.data).all()
   True

Two-Column Format
~~~~~~~~~~~~~~~~~

Test two-column PSD output::

   >>> from minke.noise import AdvancedLIGO
   >>> noise = AdvancedLIGO()
   >>> psd_array = noise.twocolumn(lower_frequency=10, upper_frequency=20, df=1)
   >>> psd_array.shape[1]  # Should have 2 columns
   2
   >>> psd_array.shape[0]  # Should have frequency samples
   11

Time Domain Representation
~~~~~~~~~~~~~~~~~~~~~~~~~~

Test covariance matrix generation::

   >>> import numpy as np
   >>> from minke.noise import AdvancedLIGO
   >>> noise = AdvancedLIGO()
   >>> times = np.linspace(0, 1, 128)
   >>> cov = noise.covariance_matrix(times)
   >>> cov.shape
   (128, 128)
   >>> # Covariance matrix should be symmetric
   >>> np.allclose(cov, cov.T)
   True

Detector Module Tests
---------------------

Test detector imports::

   >>> from minke.detector import AdvancedLIGOHanford, AdvancedLIGOLivingston
   >>> h1 = AdvancedLIGOHanford()
   >>> l1 = AdvancedLIGOLivingston()
   >>> h1.abbreviation
   'H1'
   >>> l1.abbreviation
   'L1'

Types Module Tests
------------------

Test TimeSeries creation::

   >>> import numpy as np
   >>> from minke.types import TimeSeries
   >>> data = np.random.randn(100)
   >>> times = np.linspace(0, 1, 100)
   >>> ts = TimeSeries(data=data, times=times)
   >>> len(ts.data)
   100
   >>> len(ts.times)
   100

Test PSD creation::

   >>> import numpy as np
   >>> from minke.types import PSD
   >>> frequencies = np.array([10.0, 20.0, 30.0])
   >>> psd_data = np.array([1e-46, 2e-46, 3e-46])
   >>> psd = PSD(psd_data, frequencies=frequencies)
   >>> len(psd.data)
   3
   >>> len(psd.frequencies)
   3

Waveform Model Tests
--------------------

Test CBC model imports::

   >>> from minke.models.cbc import IMRPhenomXPHM
   >>> model = IMRPhenomXPHM()
   >>> type(model).__name__
   'IMRPhenomXPHM'

Integration Tests
-----------------

Test noise generation with different sample rates::

   >>> from minke.noise import AdvancedLIGO
   >>> noise = AdvancedLIGO()
   >>> ts_4096 = noise.time_series(duration=1, sample_rate=4096)
   >>> ts_8192 = noise.time_series(duration=1, sample_rate=8192)
   >>> len(ts_4096.data)
   4096
   >>> len(ts_8192.data)
   8192

Test that noise has approximately zero mean::

   >>> import numpy as np
   >>> from minke.noise import AdvancedLIGO
   >>> noise = AdvancedLIGO()
   >>> ts = noise.time_series(duration=10, sample_rate=1024)
   >>> abs(np.mean(ts.data)) < 1e-10  # Mean should be very close to zero
   True

Test PSD consistency::

   >>> import numpy as np
   >>> from minke.noise import AdvancedLIGO
   >>> noise = AdvancedLIGO()
   >>> # Get PSD twice - should be the same (deterministic)
   >>> psd1 = noise.frequency_domain(lower_frequency=20, upper_frequency=100, df=1)
   >>> psd2 = noise.frequency_domain(lower_frequency=20, upper_frequency=100, df=1)
   >>> np.allclose(psd1.data, psd2.data)
   True

Edge Cases
----------

Test with very short duration::

   >>> from minke.noise import AdvancedLIGO
   >>> noise = AdvancedLIGO()
   >>> ts = noise.time_series(duration=0.1, sample_rate=1024)
   >>> len(ts.data)
   102

Test with single frequency::

   >>> import numpy as np
   >>> from minke.noise import AdvancedLIGO
   >>> noise = AdvancedLIGO()
   >>> psd = noise.frequency_domain(frequencies=np.array([100.0]))
   >>> len(psd.data)
   1
   >>> psd.data[0] > 0
   True
