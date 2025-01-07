from minke.noise import AdvancedLIGO
import astropy.units as u

from minke.models.cbc import IMRPhenomPv2
from minke.detector import AdvancedLIGOHanford

noise = AdvancedLIGO()

noise_ts = noise.time_series(duration=4, sample_rate=16384, epoch=998)

model = IMRPhenomPv2()

parameters = {"mass_ratio": 0.7, "total_mass": 100*u.solMass, "luminosity_distance": 100*u.megaparsec, "gpstime": 1000}

data = model.time_domain(parameters, times=noise_ts.times)
f = data['plus'].plot()
f.savefig("waveform.png")

detector = AdvancedLIGOHanford()

projected = data.project(detector,
             ra=1, dec=0.5,
             iota=0.4,
             phi_0=0,
             psi=0
             )
injection = (noise_ts + projected)
f = injection.plot()
f.savefig("projected_injection.png")
injection.write("test_injection.gwf", format="gwf")
