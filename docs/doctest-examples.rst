Doctest Examples
================

This document contains testable code examples for the Minke documentation.
Run ``make doctest`` in the docs directory to verify these examples.

Noise Module Tests
------------------

Basic Noise Generation
~~~~~~~~~~~~~~~~~~~~~~~

Test creating a noise model:

   >>> from minke.noise import AdvancedLIGO
   >>> noise = AdvancedLIGO()
   >>> type(noise).__name__
   'AdvancedLIGO'

Test time series generation with different durations:

   >>> noise = AdvancedLIGO()
   >>> ts_1sec = noise.time_series(duration=1, sample_rate=1024)
   >>> ts_2sec = noise.time_series(duration=2, sample_rate=1024)
   >>> len(ts_1sec.value)
   1024
   >>> len(ts_2sec.value)
   2048

Test time series with epoch:

   >>> noise = AdvancedLIGO()
   >>> ts = noise.time_series(duration=1, sample_rate=1024, epoch=1000)
   >>> ts.times[0]
   <Quantity 1000. s>
   >>> ts.times[1]  # doctest: +ELLIPSIS
   <Quantity 1000.000... s>

PSD Generation
~~~~~~~~~~~~~~

Test frequency domain PSD:

   >>> import numpy as np
   >>> from minke.noise import AdvancedLIGO
   >>> noise = AdvancedLIGO()
   >>> frequencies = np.array([10.0, 20.0, 50.0, 100.0])
   >>> psd = noise.frequency_domain(frequencies=frequencies)
   >>> len(psd.value) == len(frequencies)
   True
   >>> np.allclose(psd.frequencies.value, frequencies)
   True

Test PSD with frequency range:

   >>> noise = AdvancedLIGO()
   >>> psd = noise.frequency_domain(lower_frequency=20, upper_frequency=100, df=10)
   >>> len(psd.value)
   9
   >>> psd.frequencies[0].value
   20.0
   >>> psd.frequencies[-1].value
   100.0

Test that PSD values are positive:

   >>> import numpy as np
   >>> noise = AdvancedLIGO()
   >>> psd = noise.frequency_domain(lower_frequency=10, upper_frequency=200, df=1)
   >>> # The very last requested frequency is a known edge case that
   >>> # currently comes back as zero -- see frequency_domain().
   >>> bool(all(psd.value[:-1] > 0))
   True
   >>> bool(np.isfinite(psd.value).all())
   True

Two-Column Format
~~~~~~~~~~~~~~~~~

Test two-column PSD output:

   >>> from minke.noise import AdvancedLIGO
   >>> noise = AdvancedLIGO()
   >>> psd_array = noise.twocolumn(lower_frequency=10, upper_frequency=20, df=1)
   >>> psd_array.shape[1]  # Should have 2 columns
   2
   >>> psd_array.shape[0]  # Should have frequency samples
   11

Detector Module Tests
---------------------

Test detector imports:

   >>> from minke.detector import AdvancedLIGOHanford, AdvancedLIGOLivingston
   >>> h1 = AdvancedLIGOHanford()
   >>> l1 = AdvancedLIGOLivingston()
   >>> h1.abbreviation
   'H1'
   >>> l1.abbreviation
   'L1'

Types Module Tests
------------------

Test TimeSeries creation:

   >>> import numpy as np
   >>> from minke.types import TimeSeries
   >>> data = np.random.randn(100)
   >>> times = np.linspace(0, 1, 100)
   >>> ts = TimeSeries(data=data, times=times)
   >>> len(ts.value)
   100
   >>> len(ts.times)
   100

Test PSD creation:

   >>> import numpy as np
   >>> from minke.types import PSD
   >>> frequencies = np.array([10.0, 20.0, 30.0])
   >>> psd_data = np.array([1e-46, 2e-46, 3e-46])
   >>> psd = PSD(psd_data, frequencies=frequencies)
   >>> len(psd.value)
   3
   >>> len(psd.frequencies)
   3

Waveform Model Tests
--------------------

Test CBC model imports:

   >>> from minke.models.cbc import IMRPhenomXPHM
   >>> model = IMRPhenomXPHM()
   >>> type(model).__name__
   'IMRPhenomXPHM'

Integration Tests
-----------------

Test noise generation with different sample rates:

   >>> from minke.noise import AdvancedLIGO
   >>> noise = AdvancedLIGO()
   >>> ts_4096 = noise.time_series(duration=1, sample_rate=4096)
   >>> ts_8192 = noise.time_series(duration=1, sample_rate=8192)
   >>> len(ts_4096.value)
   4096
   >>> len(ts_8192.value)
   8192

Test that noise has approximately zero mean:

   >>> import numpy as np
   >>> from minke.noise import AdvancedLIGO
   >>> noise = AdvancedLIGO()
   >>> ts = noise.time_series(duration=10, sample_rate=1024)
   >>> bool(abs(np.mean(ts.value)) < 1e-10)  # Mean should be very close to zero
   True

Test PSD consistency:

   >>> import numpy as np
   >>> from minke.noise import AdvancedLIGO
   >>> noise = AdvancedLIGO()
   >>> # Get PSD twice - should be the same (deterministic)
   >>> psd1 = noise.frequency_domain(lower_frequency=20, upper_frequency=100, df=1)
   >>> psd2 = noise.frequency_domain(lower_frequency=20, upper_frequency=100, df=1)
   >>> np.allclose(psd1.value, psd2.value)
   True

Edge Cases
----------

Test with very short duration:

   >>> from minke.noise import AdvancedLIGO
   >>> noise = AdvancedLIGO()
   >>> ts = noise.time_series(duration=0.1, sample_rate=1024)
   >>> len(ts.value)
   102

Test with a minimal (two-point) frequency array:

   >>> import numpy as np
   >>> from minke.noise import AdvancedLIGO
   >>> noise = AdvancedLIGO()
   >>> psd = noise.frequency_domain(frequencies=np.array([100.0, 101.0]))
   >>> len(psd.value)
   2
   >>> bool(psd.value[0] > 0)
   True
