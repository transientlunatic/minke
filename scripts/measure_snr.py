import minke
from minke import mdctools, sources, distribution, noise
import lalsimulation

sg = sources.SineGaussian(10, 1000, 1e-23, "linear", 1126620016)

o1psd = noise.PSD("/home/daniel/LIGO-P1200087-v18-AdV_EARLY_HIGH.txt")


print lalsimulation.MeasureSNR(sg._generate_for_detector(['L1']), o1psd.psd)
