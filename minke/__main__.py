import SEOBNRE
import pycbc.noise.reproduceable
import pycbc.psd.analytical
from pycbc.waveform import get_td_waveform
import numpy as np
import lal
import lalsimulation as lalsim
from matplotlib import pyplot as plt
from scipy import interpolate
from pycbc.detector import Detector
from math import ceil, floor

from asimov.utils import update

import os

import pycbc.noise
import pycbc.psd
from gwpy.timeseries import TimeSeries

import click
import yaml

import otter
import otter.bootstrap as bt
import otter.plot

def generate_interpolant(approximant="SEOBNRE", **kwargs):
    defaults = {'mass1': 10,
              'mass2': 10,
              'spin1z':0.0,
              'spin2z':0.0,
              'distance':400,
              'delta_t':1./16384,
              'eccentricity':0,
              'coa_phase':0,
              'f_lower':20,
              'mass_ratio': 1,
              'long_asc_nodes':0}
    defaults.update(kwargs)
    if "total_mass" in kwargs:
        defaults['mass1'] = (defaults['total_mass'])/(1+defaults['mass_ratio'])
        defaults['mass2'] = (defaults['mass_ratio'] * defaults['total_mass'])/(1+defaults['mass_ratio'])
        del(defaults['total_mass'], defaults['mass_ratio'])
    hp, hc = get_td_waveform(approximant=approximant,
                                     **defaults)
    return hp, hc

def create_frames(settings):

    defaults = {"injection": {"noise":
                              {"H1": pycbc.psd.aLIGOZeroDetHighPower,
                               "L1": pycbc.psd.aLIGOZeroDetHighPower,
                               "V1": pycbc.psd.Virgo}
                              }
                }
    settings = update(defaults, settings)
    
    if "report" in settings:
        webdir = settings['report']['location']
        os.makedirs(webdir, exist_ok=True)
        report = otter.Otter(
            f"{webdir}/index.html",
            author="Asimov",
            title="Injections Report",
            #theme_location=resource_filename("asimov.cli", "report-theme"),
            #config_file=os.path.join(".asimov", "asimov.conf"),
        )
    else:
        report = None
    
    event_name = settings['event']
    
    parameters = settings['injection']
    epoch = parameters.pop('epoch')
    
    approximant = settings['waveform']['approximant']
    sample_rate = settings['likelihood']['sample rate']
    
    hp, hc = generate_interpolant(approximant=approximant, sample_rate=sample_rate,
                                  **parameters)
    
    hp.start_time += epoch
    hc.start_time += epoch

    waveforms = {}
    detectors = settings['interferometers']
    for detector in detectors:
        waveforms[detector] = Detector(detector).project_wave(hp, hc,
                                                              settings.get("injection", {}).get("right ascension", None),
                                                              settings.get("injection", {}).get("declination", None),
                                                              settings.get("injection", {}).get("polarization", None))

        segment_length = settings.get("data", {}).get("segment length", None)
        # The color of the noise matches a PSD which you provide
        flow = settings.get('quality', {}).get('minimum frequency', {}).get(detector, None)
        delta_f = 1.0 / segment_length
        flen = int(sample_rate / delta_f) + 1
        psd = settings['injection']['noise'][detector](flen, delta_f, flow)

        delta_t = 1.0 / sample_rate

        strain_pad = pycbc.types.TimeSeries(
                pycbc.types.zeros(segment_length*sample_rate, dtype=waveforms[detector].dtype),
                delta_t=waveforms[detector].delta_t, copy=False, epoch=epoch - segment_length/2)

        tsamples = len(waveforms[detector].data)
        ts = pycbc.noise.noise_from_psd(sample_rate * segment_length, delta_t, psd, seed=100)
        ts.start_time += epoch - (segment_length/2)

        if len(strain_pad.sample_times) < len(waveforms[detector].sample_times):
            # We'll need to trim the ends off the waveform
            trim_left = sum(waveforms[detector].sample_times < strain_pad.sample_times[0])
            trim_right = -sum(waveforms[detector].sample_times > strain_pad.sample_times[-1]) + 1
            strain_pad = waveforms[detector][trim_left:trim_right]
        else:
            # We'll need to pad the waveform
            pad_left = floor(((sample_rate * segment_length) - len(waveforms[detector].sample_times))/2)
            pad_right = -ceil(((sample_rate * segment_length) - len(waveforms[detector].sample_times))/2)
            strain_pad[pad_left:pad_right] = waveforms[detector]

        gwpy_ts = TimeSeries(data=ts.data+strain_pad.data, times=strain_pad.sample_times, channel=f"{detector}:Injection", name=f"{detector}:Injection", dtype=np.float64)

        directory = os.path.join(settings['rundir'], "frames")
        os.makedirs(directory, exist_ok=True)
        gwpy_ts.write(os.path.join(directory, f"{detector}-injection.gwf"), format="gwf")

        if report:
            with report:
                report += f"#{detector}"

            with report:
                f, ax = plt.subplots(1, 1, dpi=300)
                ax.plot(ts)
                ax.plot(strain_pad)

                report += "## Waveform"
                report += f

            with report:
                specgram  = gwpy_ts.spectrogram2(fftlength=1/16., overlap=15/256.) ** (1/2.)
                specgram = specgram.crop_frequencies(20)
                plot = specgram.plot(norm='log', cmap='viridis', yscale='log')
                report += "## Spectrogram"
                report += otter.plot.Figure(plot)

        directory = os.path.join(settings['rundir'], "psds")
        os.makedirs(directory, exist_ok=True)
        psd.save(os.path.join(directory, f'{ifo}-psd.txt'))
                
        directory = os.path.join(settings['rundir'], "cache")
        os.makedirs(directory, exist_ok=True)
                
        with open(os.path.join(directory, f"{detector}.cache"), "w") as f:
           absolute = os.path.abspath(f"{detector}-injection.gwf")
           f.write(f"{detector[0]}\t{detector}_Injection\t{int(epoch-segment_length/2)}\t{segment_length}\tfile://localhost{absolute}\n")

@click.command()
@click.option('--settings', help='The injection settings')            
def main(settings):
    with open(settings, "r") as settings_file:
        settings = yaml.safe_load(settings_file)
    create_frames(settings)

if __name__ == '__main__':
    main()
