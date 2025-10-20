Generic burst waveforms
=======================

Perhaps the simplest waveforms which Minke is able to produce are so-caled "generic" or "ad-hoc" burst waveforms.
These waveform families include gaussian bursts, sine-gaussian wavelets, and white noise bursts.

Gaussian bursts
---------------

Perhaps the simplest conceivable model of a burst of gravitational waves is one where energy is emitted across a broadband range of frequencies over a fixed period of time, with a smooth rise and decay in amplitude.
Such a source can be modelled as with a Gaussian function, and may be a suitable model for broadband sources, such as the core-bounce during a core-collapse abbr:sn.

In searches the model for such a signal is

.. math::
   h(t) = A \exp\left( - \frac{ (t - t_{0})^{2} }{ 2 \sigma^{2} } \right),

for a strain :math:`h` at time :math:`t`, with an amplitude :math:`A`, central time :math:`t_{0}` and duration :math:`\sigma`.

Minke supports Gaussian bursts using the `minke.sources.Gaussian` class.

.. autoclass:: minke.sources.Gaussian


Sine-Gaussian bursts
--------------------

In addition to searching for broadband, time-constrained bursts of gravitational wave energy, some sources are expected to produce gravitational waves which are in a confined range of frequencies, in addition to being released over a short time-span.

Such a source can be approximated by a sinusoidal signal which is enveloped by a Gaussian rise and decay in amplitude.

The model used in gls:ligo searches for such signals is: 

.. math::

   h(t) = A \exp \left[ \frac{ - 2(t - t_{0})^{2} \pi^{2} f^{2}}{Q^{2}} \right] \cos\left[ 2 \pi f (t - t_{0}) \right],
   
for a strain :math:`h` at time :math:`t`, with :math:`A` the amplitude of the signal, :math:`t_{0}` its central time, :math:`Q` the quality factor of the burst, and :math:`f` is frequency.


A SineGaussian burst can be produced with a short script such as this ::

  import minke
  import astropy.units as u

  from minke.models.bursts import SineGaussian
  from minke.detector import AdvancedLIGOHanford

  model = SineGaussian()

  parameters = {"centre_frequency": 20,
                "phase": 0,
                "eccentricity": 0,
                "q": 1.,
                "sample_rate": 4096 * u.Hertz,
                "gpstime": 998,
                "hrss": 1e-22,
                "duration": 2*u.second}

  data = model.time_domain(parameters)

  detector = AdvancedLIGOHanford()

  projected = data.project(detector,
                ra=1, dec=0.5,
                iota=0.4,
                phi_0=0,
                psi=0
                )

Visualizing the projected signal:

.. plot::
   :include-source:

   from minke.models.burst import SineGaussian
   from minke.detector import AdvancedLIGOHanford
   import astropy.units as u

   model = SineGaussian()
   parameters = {"centre_frequency": 20,
                 "phase": 0,
                 "eccentricity": 0,
                 "q": 1.,
                 "sample_rate": 4096 * u.Hertz,
                 "gpstime": 998,
                 "hrss": 1e-22,
                 "duration": 2*u.second}

   data = model.time_domain(parameters)
   detector = AdvancedLIGOHanford()
   projected = data.project(detector,
                 ra=1, dec=0.5,
                 iota=0.4,
                 phi_0=0,
                 psi=0
                 )

   fig = projected.plot()
   plt.tight_layout()

  
Band-limited white noise bursts
-------------------------------

Astrophysical processes are unlikely to produce emission at a single frequency, or with a smooth evolution of amplitude, and so searches are normally expected to be sensitive to band-limited white noise bursts, which consist of band-limited uncorrelated noise within a Gaussian amplitude envelope.

.. autoclass:: minke.sources.WhiteNoiseBurst
