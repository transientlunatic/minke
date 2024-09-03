"""
Code to create injections using the various models supported by heron.
"""

import logging
import os

import click

import numpy as np

import astropy.units as u

import otter

from .models.lalsimulation import SEOBNRv3, IMRPhenomPv2
from .models.lalnoise import KNOWN_PSDS
from .detector import KNOWN_IFOS
from .utils import load_yaml

logger = logging.getLogger("minke.injection")


def make_injection(
        waveform=IMRPhenomPv2,
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

    injections = {}
    for detector, psd_model in detectors.items():
        logger.info(f"Making injection for {detector}")
        psd_model = KNOWN_PSDS[psd_model]()
        detector = KNOWN_IFOS[detector]()
        if times is not None:
            kwargs = {"times": times}
        else:
            kwargs = {"duration": duration, "sample_rate": sample_rate, "epoch": epoch}
        data = psd_model.time_series(**kwargs)

        channel_n = f"{detector.abbreviation}:{channel}"
        
        waveform = waveform_model.time_domain(
            parameters,
            times=data.times
        )
        injection = data + waveform.project(detector)
        injection.channel = channel_n

        print("length of injection", len(injection.data))
        print("duration of injection", (injection.times[-1] - injection.times[0]))
        
        if psdfile:
            # Write the PSD file to an ascii file
            filename = f"{detector.abbreviation}_{psdfile}.dat"
            psd_model.to_file(filename)

        
        if framefile:
            filename = f"{detector.abbreviation}_{framefile}.gwf"
            logger.info(f"Saving framefile to {filename}")
            injection.write(filename, format="gwf")
        injections[detector.abbreviation] = injection
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

    with report:
        report + "# Injection frames"
        report + settings
    
    detector_dict = {
        settings["interferometers"][ifo]: settings["psds"][ifo]
        for ifo in settings["interferometers"]
    }
    injections = make_injection(
        channel=settings['channel'],
        duration=settings['duration'],
        sample_rate=settings['sample_rate'],
        epoch=settings['epoch'],
        waveform=IMRPhenomPv2,
        injection_parameters=parameters,
        detectors=detector_dict,
        framefile=settings['channel'],
        psdfile="psd",
    )
    data = injections

    for injection in injections:
        f = injection.plot()
        with report:
            report + f
