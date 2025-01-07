import minke
import astropy.units as u

from minke.models.cbc import IMRPhenomPv1
from minke.detector import AdvancedLIGOHanford

model = IMRPhenomPv1()

parameters = {"mass_ratio": 0.7, "total_mass": 100*u.solMass, "luminosity_distance": 10*u.megaparsec}

data = model.time_domain(parameters)

data = model.time_domain(mass_ratio=0.7, total_mass=100*u.solMass)

f = data['plus'].plot()
f.savefig("docs/images/cbc/waveform-imrphenompv1.png")

detector = AdvancedLIGOHanford()

projected = data.project(detector,
             ra=1, dec=0.5,
             iota=0.4,
             phi_0=0,
             psi=0
             )

f = projected.plot()
f.savefig("projected_waveform.png")
