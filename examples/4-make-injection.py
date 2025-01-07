from minke.injection import make_injection_zero_noise
from minke.noise import AdvancedLIGO
from minke.models.cbc import IMRPhenomPv2
from minke.detector import AdvancedLIGOHanford

import astropy.units as u


detectors = {"AdvancedLIGOHanford": "AdvancedLIGO"}

parameters = {"mass_ratio": 0.7,
              "total_mass": 100*u.solMass,
              "luminosity_distance": 100*u.megaparsec,
              "ra": 1,
              "dec": 0.5,
              "iota": 0.4,
              "phi_0": 0,
              "psi": 0,
              "gpstime": 1000}

injection = make_injection_zero_noise(detectors=detectors, injection_parameters=parameters, duration=4, sample_rate=16384, epoch=998)['H1']

f = injection.plot()
f.savefig("injection_function_zero.png")
