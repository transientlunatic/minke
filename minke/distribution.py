from scipy import random
import numpy

def sky_params(net, time, ra, dec, psi, incl=numpy.pi):
    if not (0 < ra < 2 * numpy.pi) or not (0 < psi < 2 * numpy.pi) or not(-numpy.pi/2 < dec < numpy.pi/2):
        return 0.0

    fp, fx = 0, 0
    for ifo in net:
        tp, tx, _, _ = response(time, ra, dec, incl, psi, 'radians', ifo)
        fp += tp**2
        fx += fx**2
    fp, fx = math.sqrt(fp), math.sqrt(fx)
    return (fp + fx) / len(net)

def uniform_dec(num):
    """
    Declination distribution: uniform in sin(dec). num controls the number of draws.
    """
    return (numpy.pi / 2.) - numpy.arccos(2 * random.random_sample(num) - 1)

def uniform_theta(num):
    """
    Uniform in cos distribution. num controls the number of draws.
    """
    return numpy.arccos(2 * random.random_sample(num) - 1)

def uniform_phi(num):
    """
    Uniform in (0, 2pi) distribution. num controls the number of draws.
    """
    return random.random_sample(num) * 2 * numpy.pi

def uniform_interval(interval, num):
    """
    Uniform in an interval, with the interval being a 2-tuple. num controls the number of draws. See also:

    http://docs.scipy.org/doc/numpy/reference/generated/numpy.random.random_sample.html#numpy.random.random_sample
    """
    return (interval[1] - interval[0]) * random.random_sample(num) + interval[0]
