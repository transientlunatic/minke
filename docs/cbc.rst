Compact Binary Waveforms
========================

Compact binary coalescences are the sources which have produced all of the gravitational wave detections to the end of 2023.
Minke provides a way of producing waveforms from the various models which are supported by the ``lalsimulation`` library.


Binary Black Holes
------------------

IMRPhenomPv2
^^^^^^^^^^^^

IMRPhenomPv2 is a precessing, spinning black hole model.

.. plot::
   :include-source:

   import minke
   import astropy.units as u
   from minke.models.cbc import IMRPhenomPv2

   model = IMRPhenomPv2()
   parameters = {"mass_ratio": 0.7, 
                 "total_mass": 100*u.solMass, 
                 "luminosity_distance": 10*u.megaparsec}

   data = model.time_domain(parameters)

   fig = data['plus'].plot()
   plt.tight_layout()

	   
SEOBNRv3
^^^^^^^^

.. plot::
   :include-source:

   import minke
   import astropy.units as u
   from minke.models.cbc import SEOBNRv3

   model = SEOBNRv3()
   parameters = {"mass_ratio": 0.7, 
                 "total_mass": 100*u.solMass, 
                 "luminosity_distance": 10*u.megaparsec}

   data = model.time_domain(parameters)

   fig = data['plus'].plot()
   plt.tight_layout()
      
