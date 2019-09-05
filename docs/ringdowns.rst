Ringdown Waveforms
==================

Minke is capable of producing signals which are based on waveforms
based on Binary Black Hole (BBH) ringdowns.

Ringdown-like signals, with a sudden rise, and exponential decay in amplitude are expected in the post-merger signal of CBC systems, and in some models of neutron star model excitation [2004PhRvD..70l4015B].

These take the form

.. math::
   h(t) = \exp (-t / \tau) \sin( 2 \pi f t)
   
for a strain :math:`h` at time :math:`t`, given a decay time :math:`\tau` and frequency :math:`f`.

BBHRingdown
-----------

The BBH Ringdown class in the sources module follows a similar style
to other source classes within `minke`. To produce a single ringdown
waveform:

>>> import minke.sources
>>> bbh = minke.sources.BBHRingdown(100.23, 1e-22, np.rad2deg(.1), 10., 0.97, 0.01, 1.0,45.)

then the waveforms themselves can be produced by running

>>> hp, hx, _, _ = bbh._generate()

The waveform will be generated at the LIGO sampling rate (16384 Hz) by
default, but this can be changed by specifying a `rate` keyword. The
waveforms are also generated for the (2,2) mode by default, but by
adding `l` and `m` keywords to the `_generate()` method this can also be changed.


.. autoclass:: minke.sources.BBHRingdown
