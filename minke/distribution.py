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

def supernova_angle(num, divisions = 10):
    """
    Draw from a discrete distribution of angles. 
    """
    theta = numpy.linspace(-1, 1, divisions)
    ctheta = numpy.arccos(theta)
    phi = numpy.linspace(0, numpy.pi*2, divisions)
    out_t =  numpy.random.choice(ctheta, num, replace = True)
    out_p =  numpy.random.choice(phi, num, replace = True)
    return zip(out_t, out_p)

def uniform_sky(number=1):
    """
    Get a set of (RA, declination, polarization) randomized appopriately to astrophysical sources isotropically distributed in the sky.
    """
    expnum = number
    ra = uniform_phi(expnum)
    dec = uniform_dec(expnum)
    pol = uniform_phi(expnum)
    return ra, dec, pol

def favorable_sky(net, time):
    """
    Wander through the skies, searching for a most favorable location --- draw extrinsic parameters as if the network antenna pattern magnitude were the PDF.
    """
    ndraws = len(time)
    ra_out, dec_out, psi_out = numpy.empty((3, len(time)))
    while ndraws > 0:
        ra = numpy.random.uniform(0, 2 * numpy.pi, 1000)
        dec = numpy.random.uniform(-numpy.pi / 2, numpy.pi / 2, 1000)
        psi = numpy.random.uniform(0, 2 * numpy.pi, 1000)
        rnd = numpy.random.uniform(0, 1, 1000)
        for r, d, p, n in zip(ra, dec, psi, rnd):
            if n < sky_params(net, time, r, d, p):
                ndraws -= 1
                ra_out[ndraws] = r
                dec_out[ndraws] = d
                psi_out[ndraws] = p
                break
    return ra_out, dec_out, psi_out

def uniform_interval(interval, num):
    """
    Return a number, or a list of numbers which are sampled
    from a uniform distribution.

    Parameters
    ----------
    interval : tuple (lower, upper)
       The interval which the uniform distribution covers.
    num : int
       The number of samples which should be drawn from the distribution.

    Returns
    -------
    sample : float or list
       A sample, or a list of samples.

    Notes
    -----
    
    http://docs.scipy.org/doc/numpy/reference/generated/numpy.random.random_sample.html#numpy.random.random_sample
    """
    return (interval[1] - interval[0]) * random.random_sample(num) + interval[0]

def uniform_time(tstart, tstop, number):
    """
    Get a set of randomized (integer) event times.
    """
    return random.randint(tstart, tstop, number) + random.rand(number)

def log_uniform(lower, upper, number):
    """
    Draw uniformly in the log of a predefined  range.

    Parameters
    ----------
    lower : float
       The lower hrss value (n.b. not the lower log(hrss) value)
    upper : float
       The upper hrss value (n.b. not the upper log(hrss) value)
    number : int
       The number of samples to be drawn from the distribution.

    Returns
    -------
    sample : float
       An sample value drawn from the log uniform distribution.
    """
    log10h = (numpy.log10(lower), numpy.log10(upper))
    return 10**uniform_interval(log10h, number)

def even_time(start, stop, rate, jitter=0):
    """
    Produce an evenly-distributed set of times which has some jitter.

    Parameters
    ----------
    start : int
       The gps start time of the set.
    stop : int
       The gps end time of the set.
    rate : int
       The rate of events per year.
    jitter : float
       The number of seconds around which the time can "jitter".
       This has the effect of adding a small random time on to each
       time returned in the distribution, up to a maximum of half 
       the jitter value, or taking away up to half the jitter value.

    Returns
    -------
    times : ndarray
       A numpy array of gps times where injections should be made.
    """
    expnum_exact = (stop-start) * rate / 365.0 / 24 / 3600
    interval = (stop - start) / expnum_exact
    time = numpy.array(map(lambda i: start + i*interval + jitter *(random.rand()-0.5), range(1,int(expnum_exact)+1)))
    return time

def burst_dist(minimum, maximum, size=1):
    """
    Generate an hrss drawn from the distribution
    \[ r + 50/r \]
    as desired by the Burst group for observing MDCs
    """
    #extent = numpy.log10(maximum) - numpy.log10(minimum)
    # Convert hrss to distance
    minimum, maximum = 1/minimum, 1/maximum
    #Set-up the distributions limits
    a = (minimum + 50.0/minimum)**2
    b = (maximum + 50.0/maximum)**2
    rnd = numpy.random.uniform(a,b, size)
    # Convert back to an hrss
    return 1.0/numpy.sqrt(rnd)
