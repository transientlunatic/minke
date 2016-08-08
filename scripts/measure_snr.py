import minke
from minke import mdctools, sources, distribution, noise
import lalsimulation

sg = sources.SineGaussian(10, 1000, 1e-22, "linear", 1126620016)

o1psd = noise.PSD("LIGO-P1200087-v18-AdV_EARLY_HIGH.txt")

print o1psd.SNR(sg, ['L1', 'H1'])
