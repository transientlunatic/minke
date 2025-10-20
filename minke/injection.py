"""
Code to create injections using the various models supported by heron.
"""

import logging
import os

import click

import numpy as np
from scipy.optimize import brentq

import astropy.units as u

import otter
import matplotlib.pyplot as plt

from .models.lalsimulation import SEOBNRv3, IMRPhenomPv2, IMRPhenomXPHM
from .models.lalnoise import KNOWN_PSDS
from .detector import KNOWN_IFOS
from .utils import load_yaml
from .filters import inner_product

logger = logging.getLogger("minke.injection")

def calculate_network_snr_for_distance(distance, waveform_model, parameters, detectors, psd_models, times):
    """Calculate the network SNR for a given luminosity distance."""
    params = parameters.copy()
    params['luminosity_distance'] = distance * u.megaparsec
    
    waveform = waveform_model.time_domain(params, times=times)
    
    network_snr_squared = 0.0
    sample_rate = 1.0 / (times[1] - times[0])
    
    for detector, psd_model in zip(detectors, psd_models):
        injection_data = waveform.project(detector)
        
        # Calculate SNR for this detector
        injection_data_f = np.fft.fft(injection_data.data, n=len(times)//2) / sample_rate
        
        N = len(times)
        df = 1. / sample_rate
        frequencies = np.arange(0, N // 2) * df
        psd_f = psd_model.frequency_domain(frequencies=frequencies)
        
        snr_squared = inner_product(injection_data_f, injection_data_f, np.array(psd_f.data))
        network_snr_squared += snr_squared
    
    return np.sqrt(network_snr_squared)


def find_distance_for_network_snr(target_snr, waveform_model, parameters, detectors, psd_models, times):
    """Find the luminosity distance that produces the target network SNR."""
    
    def snr_residual(distance):
        return calculate_network_snr_for_distance(distance, waveform_model, parameters, 
                                                  detectors, psd_models, times) - target_snr
    
    # Search between 10 Mpc and 10000 Mpc
    try:
        distance = brentq(snr_residual, 10, 10000, xtol=0.1)
        return distance * u.megaparsec
    except ValueError:
        snr_at_10 = calculate_network_snr_for_distance(10, waveform_model, parameters, 
                                                       detectors, psd_models, times)
        logger.error(f"Could not find distance for network SNR={target_snr}. Network SNR at 10 Mpc = {snr_at_10:.1f}")
        raise

def make_injection(
        waveform=IMRPhenomXPHM,
        injection_parameters={},
        times=None,
        epoch=None,
        detectors=None,
        framefile=None,
        duration=None,
        sample_rate=None,
        psdfile=None,
        channel=None,
):

    if channel is None:
        channel = "Injection"
    
    parameters = {"ra": 0, "dec": 0, "psi": 0, "theta_jn": 0, "phase": 0}
    parameters.update(injection_parameters)

    waveform_model = waveform()

    default_units = {"m1": u.solMass,
                     "m2": u.solMass,
                     "luminosity_distance": u.megaparsec,
                     }

    for p in parameters.keys():
        if p in default_units.keys() and not isinstance(parameters[p], u.Quantity):
            parameters[p] *= default_units[p]


    target_snr = injection_parameters.get("snr", None)
    if target_snr is not None:
        logger.info(f"Finding luminosity distance for target network SNR={target_snr}")
        
        # Prepare detector and PSD objects for all detectors
        detector_objects = []
        psd_objects = []
        
        for det_name, psd_name in detectors.items():
            detector_objects.append(KNOWN_IFOS[det_name]())
            psd_objects.append(KNOWN_PSDS[psd_name]())
        
        # Get times array
        if times is not None:
            times_array = times
        else:
            psd_temp = psd_objects[0].time_series(duration=duration, sample_rate=sample_rate, epoch=epoch)
            times_array = psd_temp.times
        
        parameters['luminosity_distance'] = find_distance_for_network_snr(
            target_snr, waveform_model, parameters, detector_objects, psd_objects, times_array
        )
        logger.info(f"Required luminosity distance for network SNR: {parameters['luminosity_distance']:.2f}")

    injections = {}
    detector_snrs = {}
    for detector, psd_model in detectors.items():
        logger.info(f"Making injection for {detector}")
        psd_model = KNOWN_PSDS[psd_model]()
        detector = KNOWN_IFOS[detector]()
        if times is not None:
            kwargs = {"times": times}
        else:
            kwargs = {"duration": duration, "sample_rate": sample_rate, "epoch": epoch}
        data = psd_model.time_series(**kwargs)

       
        print("data length", len(data))

        channel_n = f"{detector.abbreviation}:{channel}"
        
        waveform = waveform_model.time_domain(
            parameters,
            times=data.times
        )
        injection_data = waveform.project(detector)
        injection = data + injection_data
        injection.channel = channel_n

        injection_data_f = np.fft.fft(injection_data.data, n=len(data.times)//2)/sample_rate
        
        N = len(data.times)
        df = 1./sample_rate
        frequencies = np.arange(0, N // 2) * df
        print(len(frequencies))
        
        psd_f = psd_model.frequency_domain(frequencies = frequencies)
        det_snr = np.sqrt(inner_product(injection_data_f, injection_data_f, np.array(psd_f.data)))
        detector_snrs[detector.abbreviation] = det_snr
        print(f"Optimal SNR for {detector.abbreviation}: {det_snr:.2f}")
        
        print("length of injection", len(injection.data))
        print("duration of injection", (duration))

        print("duration in TS", data.times[-1] - data.times[0])
        print("duration by calc", data.dt * len(data.data))
        print("dt", data.dt)
        
        if psdfile:
            # Write the PSD file to an ascii file
            filename = f"{detector.abbreviation}_{psdfile}.dat"
            psd_model.to_file(filename, upper_frequency=(sample_rate/2)-1./duration,
                              lower_frequency=0, df=1./duration)

        
        if framefile:
            filename = f"{detector.abbreviation}_{framefile}.gwf"
            logger.info(f"Saving framefile to {filename}")
            injection.write(filename, format="gwf.lalframe")
        injections[detector.abbreviation] = injection

    network_snr = np.sqrt(sum(snr**2 for snr in detector_snrs.values()))
    logger.info(f"Network SNR: {network_snr:.2f}")
    print(f"Network SNR: {network_snr:.2f}")
    
    return injections


def make_injection_zero_noise(
        waveform=IMRPhenomPv2,
        injection_parameters={},
        times=None,
        epoch=None,
        detectors=None,
        framefile=None,
        duration=None,
        sample_rate=None,
        channel=None
):

    if channel is None:
        channel = "Injection"
    
    parameters = {"ra": 0, "dec": 0, "psi": 0, "theta_jn": 0, "phase": 0, 'gpstime': 4000}
    parameters.update(injection_parameters)

    waveform_model = waveform()

    if times is not None:
        kwargs = {"times": times}
    else:
        kwargs = {"duration": duration, "sample_rate": sample_rate, "epoch": epoch}


    injections = {}
    for detector, psd_model in detectors.items():
        detector = KNOWN_IFOS[detector]()
        channel = f"{detector.abbreviation}:{channel}"
        logger.info(f"Making injection for {detector} in channel {channel}")
        psd_model = KNOWN_PSDS[psd_model]()
        data = psd_model.time_series(**kwargs)
        waveform = waveform_model.time_domain(
            parameters,
            times=data.times
        )
        injection = waveform.project(detector)
        injection.channel = channel
        injections[detector.abbreviation] = injection
        if framefile:
            filename = f"{detector.abbreviation}_{framefile}.gwf"
            logger.info(f"Saving framefile to {filename}")
            injection.write(filename, format="gwf")
        injections[detector.abbreviation] = injection
        return injections
            

def injection_parameters_add_units(parameters):
    UNITS = {"luminosity_distance": u.megaparsec, "m1": u.solMass, "m2": u.solMass}

    for parameter, value in parameters.items():
        if not isinstance(value, u.Quantity) and parameter in UNITS:
            parameters[parameter] = value * UNITS[parameter]
    return parameters


@click.command()
@click.option("--settings")
def injection(settings):

    settings = load_yaml(settings)

    if "logging" in settings:

        level = settings.get("logging", {}).get("level", "warning")

        LOGGER_LEVELS = {
            "info": logging.INFO,
            "debug": logging.DEBUG,
            "warning": logging.WARNING,
        }

        logging.basicConfig(level=LOGGER_LEVELS[level])
    settings_ = settings
    settings = settings["injection"]
    parameters = injection_parameters_add_units(settings["parameters"])

    report = otter.Otter(os.path.join(settings_.get("pages directory", "."), "injection.html"),
                         title="Minke Injection"
                         )

    detector_dict = {
        settings["interferometers"][ifo]: settings["psds"][ifo]
        for ifo in settings["interferometers"]
    }

    
    with report:
        report + "# PSDs"

        for detector, psd in detector_dict.items():
            psd_o = KNOWN_PSDS[psd]()
            data = psd_o.frequency_domain()
            f, ax = plt.subplots(1,1, dpi=300)
            ax.plot(data.frequencies, np.array(data.data))
            report + f 
    
    with report:
        report + "# Injection frames"
        report + settings
    
    injections = make_injection(
        channel=settings['channel'],
        duration=settings['duration'],
        sample_rate=settings['sample_rate'],
        epoch=settings['epoch'],
        waveform=IMRPhenomXPHM,
        injection_parameters=parameters,
        detectors=detector_dict,
        framefile=settings['channel'],
        psdfile="psd",
    )
    data = injections

    for injection in injections.values():
        f, ax = plt.subplots(1,1, dpi=300)
        ax.plot(injection.times, np.array(injection.data))
        
        with report:
            report + f
