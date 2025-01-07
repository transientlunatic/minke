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
              "gsptime": 998,
              "hrss": 1e-22,
              "duration": 2*u.second}

data = model.time_domain(parameters)
f = data['plus'].plot()
f.savefig("sinegaussian.png")

detector = AdvancedLIGOHanford()

projected = data.project(detector,
             ra=1, dec=0.5,
             iota=0.4,
             phi_0=0,
             psi=0
             )

f = projected.plot()
f.savefig("projected_sinegaussian.png")
