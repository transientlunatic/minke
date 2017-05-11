import sys,os,glob
import os.path
import otter
import re
import pandas

import matplotlib.pyplot as plt

import lalburst, lalsimulation, lalmetaio
from pylal.antenna import response
from logmake import *

#--------------------------------------#
# Startup information.
# ------
# This should be moved to a
# configuration file post haste
#--------------------------------------#

ifos = ['H1', 'L1']
runs = ['O1']
calibration = 'C01'

lasttime = 1136762880
#lasttime = 1137213440

import mdctools

o1 = mdctools.FrameSet('frame_list.dat')



inj_families = ['d08'] # The types of injections for which frames must be produced.

xml_folder = '/home/daniel.williams/data/mdc/O1/xml'    # The location for the XML files
mdc_folder = '/home/daniel.williams/data/mdc/O1/frames' # The place where MDC frames will be stored once they've been generated.

#report_location = '/home/daniel.williams/public_html/reports/mdc/o1'


# Look-up table to convert from the short to the long names for each injection family.
# inj_families_names = {'ga' : 'Gaussian',
#                       'sg' : 'SineGaussian',
#                       'wnb': 'BTLWNB',
#                       'd08': 'Dimmelmeier 08'
#                       }


import numpy

import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
font = {'size'   : 1}
matplotlib.rc('font', **font)

plt.style.use('ggplot')

from glue.ligolw import ligolw, utils, lsctables
lsctables.use_in(ligolw.LIGOLWContentHandler)

import lalburst, lalsimulation, lalmetaio

from pycbc.types import TimeSeries

def mkdir(path):
    sub_path = os.path.dirname(path)
    if not os.path.exists(sub_path):
        mkdir(sub_path)
    if not os.path.exists(path):
        os.mkdir(path)


report = otter.Otter(report_location+"/sn.html", 
	{'title':"O1 SN MDC",
	'author':'Daniel Williams',
	'subtitle':''})


report.write_page_header()
report.write_breadcrumb([('Burst MDCs', '../'),('O1', '#')])


frames = o1.frame_list

for family in inj_families:
    report.write_header(2,inj_families_names[family])
    report._write("<a href='{}'></a>".format(family))
    # Work on each family of injections sequentially.
    os.chdir(xml_folder)
    
    for injection in glob.glob("{}*".format(family)):
        # Each injection gets its own report
        subreport_loc="{}/{}.html".format(family,injection)
        report_inj = otter.Otter(report_location+"/"+subreport_loc, 
	                         {'title':"O1 Burst MDC",
	                          'author':'Daniel Williams',
	                          'subtitle':'{}'.format(injection)}
        )
        report_inj.write_page_header()
        report_inj.write_breadcrumb([('Burst MDCs', '../../'),('O1','../'),(inj_families_names[family], '../index.html#{}'.format(family)),(injection,'#')])

        print "Working on {}".format(injection)

        # Prepare the MDC Set
        mdcs = mdctools.MDCSet(['H1', 'L1'], injection, full=False)

        
        # Write a link into the main report for this injection family.
        report.write_row("<a href={}>{}</a>".format(subreport_loc, injection))


        # # Start the histograms of each of the file's parameters.
        report_inj.write_header(2,"Parameter histograms")
        
        fig = mdcs.plot_skymap()
        report_inj.write_plot(figure=fig)
        plt.close()

        fig = mdcs.plot_hist('hrss')
        report_inj.write_plot(figure=fig)
        plt.close()

        report_inj.write_warning('info', 'These injections can be found in {}'.format(mdc_folder+"/"+family+"/"+injection[:-16]+"/"))

        
        report_inj.write_footer()
        

report.write_footer()


