Editing an MDC XML file
========================

As well as producing new MDC specifications, Minke is capable of
editing pre-existing SimBurstTable XML files in a Pythonic way.

In this tutorial we will load a pre-existing XML file, and change the
value of every row's "hrss" value, but the same approach can be used
to edit any of the values within the SimBurstTable.

First we load the required parts of minke::
  >>> from minke import mdctools

In order to load an XML file we need to specify an MDCSet; what we
specify in this step isn't terribly important if we're just editing
XML files, but is essential if we need to produce a frame file at some
point. So if we assume we'll be working with the two LIGO 4km detectors::
  >>> mdcset = mdctools.MDCSet(['L1', 'H1'])

We'll now load our XML file into the `mdcset` object. We'll assume
they're in a file, in the present working directory, called
"injections.xml.gz"::
  >>> mdcset.load("injections.xml.gz")

This loads each of the waveforms in the table into an array called
`mdcset.waveforms`, which can be indexed and iterated over just like
any Python list, and we can access individual rows a bit like this.::
  >>> mdcset.waveforms[1]

which gives us access to the underlying Swig SimBurst object from
LALSuite. We can then edit the values of this by calling them as
attributes, so, for example, to set the hrss,::
  >>> mdcset.waveforms[1].hrss = 1e-23

In order to perform a batch edit on all of the rows we just need a loop. For example, if we wanted to set every hrss to `1e-23`,::
  >>> for i in xrange(len(mdcset.waveforms)):
  >>>    mdcset.waveforms[i].hrss = 1e-23

We then need to save this edit as an xml file,::
  >>> mdcset.save_xml('edited_injections.xml.gz')

and we're then free to use this xml file in other pipelines·
