Hardware Injections
===================

In addition to being able to produce software injections using Minke, you can produce data in a format which is suitable for injection into the hardware systems of gravitational wave detector using the photometric calibrator (PCAL).

Production of hardware injections is supported for all of the waveform families supported by Minke.

In order to produce hardware injection 'frames', which are stored as ASCII-format text files, Minke provides an `HWFrameSet` class, which can be used analogously to the normal `FrameSet` class for producing GWF-format frame files.

.. autoclass:: minke.mdctools.HWFrameSet


Example
-------

::

   from minke import mdctools, distribution, sources
   mdcset = mdctools.MDCSet(['L1', 'H1', 'V1'])


   time = 1000 

   hrss = 1e-22 # This is for a harware injection, so this is a nominal strain; it can be rescaled by multiplying all of the strain values.


   sg = sources.SineGaussian(time=time, q=100, frequency=235, hrss=hrss, polarisation="elliptical")

   mdcset + sg

   mdcset.save_xml("sg_F235Q100H1E-22-elliptical-hwinj.xml.gz")

   inj = mdctools.HWFrameSet(ifos=["H1", "L1", "V1"])
   frames = inj.full_frameset(mdc=mdcset, directory="./")
