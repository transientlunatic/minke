"""
This is a script designed to read an XML file of injection signals, and 
produce the hardware injection ASCII files described by them. 

The script takes a single command-line argument, the filepath to the xml file.
"""

import matplotlib
matplotlib.use('agg')
from minke import mdctools, distribution, sources
import sys
mdcset = mdctools.MDCSet(['L1', 'H1'])

print "Constructing the frame set"
o1 = mdctools.HWFrameSet() #str(sys.argv[2]))

print "Loading the XML"
mdcset.load_xml(sys.argv[1], full=True)
print "XML Loaded"

mdc_folder = r"/home/daniel.williams/scripts/hwinj_test/"

print "Starting frame production"
for o1frame in o1.frames:
    o1frame.generate_pcal(mdcset, mdc_folder, 'SCIENCE')

    del o1frame

#o1.full_logfile(mdcset, '/home/daniel.williams/data/mdc/O2a/graven/{}-logfile.txt'.format(sys.argv[1]))
