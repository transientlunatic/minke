Using Minke with Asimov Workflows
==================================

This tutorial demonstrates how to integrate Minke into an Asimov workflow for automated injection generation. 
`Asimov <https://asimov.docs.ligo.org/>`_ is a workflow management system designed for gravitational wave data analysis, and Minke provides a pipeline interface to generate injection frames as part of an Asimov-managed analysis.

What is Asimov?
---------------

Asimov is a workflow management and automation system used within the LIGO-Virgo-KAGRA collaboration. It helps manage complex analysis workflows by:

- Organizing analysis configurations in a git repository
- Automating job submission to computing clusters
- Tracking job status and outputs
- Managing dependencies between analysis steps
- Providing a web interface for monitoring

Minke's Asimov integration allows you to generate injection frames as part of a larger analysis workflow, making it easy to create consistent, reproducible injection sets.

Prerequisites
-------------

Before using Minke with Asimov, you'll need:

1. Asimov installed and configured
2. Minke installed in your analysis environment
3. Access to an HTCondor cluster (for job submission)
4. A git repository set up for Asimov event management

Installation
~~~~~~~~~~~~

.. code-block:: bash

   pip install asimov minke

You'll also need to configure Asimov. 
See the `Asimov documentation <https://asimov.docs.ligo.org/>`_ for detailed setup instructions.

Setting Up a Minke Production
------------------------------

In Asimov, a "production" is a unit of work that generates some output. 
To create injection frames with Minke, you'll define a Minke production in your event ledger.


Creating the Injection Blueprint
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Create a YAML configuration file in your event repository that specifies the injection parameters. This file should be named to match your production (e.g., ``injection-O1.yaml``).

.. code-block:: yaml

   # injection-O1.yaml
   injection:
     duration: 16
     sample_rate: 4096
     epoch: 1126259462  # Start time for O1
     channel: INJECTION
     
     parameters:
       # Binary black hole parameters
       m1: 36              # Solar masses
       m2: 29              # Solar masses
       luminosity_distance: 410  # Megaparsecs
       ra: 1.95            # Right ascension (radians)
       dec: -1.27          # Declination (radians)
       theta_jn: 0.4       # Inclination angle
       phase: 0            # Orbital phase
       psi: 0.82           # Polarization angle
       gpstime: 1126259478 # Merger time (must be within duration)
     
     waveform: IMRPhenomXPHM  # Waveform approximant
     
     interferometers:
       H1: AdvancedLIGOHanford
       L1: AdvancedLIGOLivingston
     
     psds:
       H1: AdvancedLIGO
       L1: AdvancedLIGO

Configuration Parameters
~~~~~~~~~~~~~~~~~~~~~~~~

Let's break down the configuration parameters:

**Injection Parameters**:

- ``duration``: Length of the data segment in seconds
- ``sample_rate``: Sampling frequency in Hz
- ``epoch``: GPS start time of the data segment
- ``channel``: Name of the channel in the output frame files

**Physical Parameters**:

- ``m1``, ``m2``: Component masses (in solar masses)
- ``luminosity_distance``: Distance to the source (in Megaparsecs)
- ``ra``, ``dec``: Sky position in equatorial coordinates (radians)
- ``theta_jn``: Inclination angle (radians)
- ``phase``: Orbital phase at coalescence
- ``psi``: Polarization angle (radians)
- ``gpstime``: GPS time of the merger (must fall within [epoch, epoch + duration])

**Waveform and Detectors**:

- ``waveform``: The waveform approximant to use (e.g., ``IMRPhenomXPHM``, ``SEOBNRv4``)
- ``interferometers``: Map of detector names to detector classes
- ``psds``: Map of detector names to PSD models

Updating the Event Ledger
~~~~~~~~~~~~~~~~~~~~~~~~~~

Add the production to your event in the Asimov ledger:

.. code-block:: console

   $ asimov apply -e GW150914_095045 --file injection-O1.yaml

Running the Production
----------------------

Using the Asimov Command Line
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Once your production is configured, you can manage it using the Asimov command-line interface:

.. code-block:: bash

   # Check the status of all productions
   asimov monitor
   
   # Apply changes to create/update the production
   asimov apply GW150914
   
   # Build and submit the job
   asimov build
   asimov submit
   
   # Monitor the job status
   asimov monitor GW150914

The Asimov build and submit process will:

1. Generate the necessary configuration files
2. Create an HTCondor submission file
3. Submit the job to the HTCondor scheduler
4. Track the job status


Understanding the Workflow
---------------------------

Behind the Scenes
~~~~~~~~~~~~~~~~~

When you submit a Minke production through Asimov:

1. **Job Creation**: Asimov creates an HTCondor job description using the ``minke.asimov.Asimov`` pipeline class

2. **Configuration**: The configuration from your event ledger is rendered into a YAML file using Minke's template

3. **Submission**: The job is submitted to HTCondor with appropriate resource requests and accounting information

4. **Execution**: HTCondor runs ``minke injection --settings <config>.yaml`` on a worker node

5. **Output Collection**: Generated frame files, cache files, and PSDs are collected and registered in the event metadata

6. **Completion**: Asimov marks the production as complete and makes outputs available to downstream productions

Output Files
~~~~~~~~~~~~

A successful Minke production generates several output files:

- **Frame files** (``*.gwf``): LIGO frame files containing the injection time series for each detector
- **Cache files** (``*.cache``): Lists of frame file locations for use with other tools
- **PSD files** (``*_psd.dat``): Two-column ASCII files with the PSDs used for noise generation

These outputs are automatically registered in the event metadata and can be used by subsequent analysis productions.

Advanced Usage
--------------

Creating Multiple Injections
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For studies requiring multiple injections with different parameters, you can create multiple productions:

.. code-block:: yaml

   productions:
     - name: injection-low-mass
       pipeline: minke
       status: ready
       meta:
         injection:
           parameters:
             m1: 10
             m2: 8
             luminosity_distance: 100
             # ... other parameters
     
     - name: injection-high-mass
       pipeline: minke
       status: ready
       meta:
         injection:
           parameters:
             m1: 80
             m2: 65
             luminosity_distance: 500
             # ... other parameters

Chaining with Analysis Productions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You can create downstream productions that depend on your injection:

.. code-block:: yaml

   productions:
     - name: injection-O1
       pipeline: minke
       status: ready
       meta:
         # ... injection config
     
     - name: pe-analysis
       pipeline: bilby
       status: wait
       needs:
         - injection-O1
       meta:
         data:
           channels:
             H1: H1:INJECTION
             L1: L1:INJECTION
         # ... other bilby config

Using Custom PSD Files
~~~~~~~~~~~~~~~~~~~~~~

If you want to use a PSD from a file instead of a LALSimulation model, you can reference it in your configuration. However, note that this requires extending Minke's PSD classes (see the noise-PSD tutorial).

.. code-block:: yaml

   psds:
     H1: /path/to/H1_psd.txt
     L1: /path/to/L1_psd.txt

Resource Requirements
~~~~~~~~~~~~~~~~~~~~~

You can adjust the computational resources requested for your job by modifying the pipeline's resource requests. This is done in the Asimov pipeline class, but you can also override in the event ledger:

.. code-block:: yaml

   meta:
     scheduler:
       request_memory: 2048  # MB
       request_disk: 4096    # MB
       accounting_group: ligo.dev.o4.cbc.pe.lalinference

Troubleshooting
---------------

Common Issues
~~~~~~~~~~~~~

**Job fails immediately**:

- Check that the ``gpstime`` parameter falls within ``[epoch, epoch + duration]``
- Verify that all required parameters are provided
- Ensure the waveform approximant name is correct

**Missing outputs**:

- Check the job error log (``<production-name>.err``)
- Verify that the output directory is writable
- Check that required software is available in the computing environment

**Incorrect sky position or detector response**:

- Verify that ``ra`` and ``dec`` are in radians, not degrees
- Check that ``psi`` (polarization) and ``theta_jn`` (inclination) are set correctly

Checking Logs
~~~~~~~~~~~~~

Asimov and HTCondor create several log files:

.. code-block:: bash

   # HTCondor logs
   cat /path/to/rundir/injection-O1.out  # Standard output
   cat /path/to/rundir/injection-O1.err  # Standard error
   cat /path/to/rundir/injection-O1.log  # HTCondor log
   
   # Asimov logs
   asimov log GW150914

Example: Complete Workflow
---------------------------

Here's a complete example workflow for creating injections and running a parameter estimation analysis:

.. code-block:: yaml

   # Event ledger entry
   name: Injection-Study-BBH
   repository: /home/username/studies/bbh-injections
   
   productions:
     # Generate injection frames
     - name: bbh-injection
       pipeline: minke
       status: ready
       comment: Binary black hole injection for PE study
       
       meta:
         injection:
           duration: 32
           sample_rate: 4096
           epoch: 1187008882
           channel: INJECTION
           
           parameters:
             m1: 50
             m2: 35
             luminosity_distance: 500
             ra: 2.5
             dec: 0.8
             theta_jn: 1.2
             phase: 0.5
             psi: 1.1
             gpstime: 1187008900
           
           interferometers:
             H1: AdvancedLIGOHanford
             L1: AdvancedLIGOLivingston
           
           psds:
             H1: AdvancedLIGO
             L1: AdvancedLIGO
         
         waveform:
           approximant: IMRPhenomXPHM
         
         scheduler:
           accounting_group: ligo.dev.o4.cbc.pe.lalinference
     
     # Run parameter estimation on the injection
     - name: bbh-pe
       pipeline: bilby
       status: wait
       needs:
         - bbh-injection
       comment: Parameter estimation on BBH injection
       
       meta:
         data:
           channels:
             H1: H1:INJECTION
             L1: L1:INJECTION
           frame_files:
             H1: '{bbh-injection:H1:frames}'
             L1: '{bbh-injection:L1:frames}'
         # ... bilby configuration

Then run:

.. code-block:: bash

   # Apply the configuration
   asimov apply Injection-Study-BBH
   
   # Build and submit the injection job
   asimov build Injection-Study-BBH bbh-injection
   asimov submit Injection-Study-BBH bbh-injection
   
   # Monitor progress
   asimov monitor Injection-Study-BBH
   
   # Once injection is complete, the PE job will automatically start
   # (because of the "needs" dependency)

Summary
-------

This tutorial covered:

- Setting up Minke as an Asimov pipeline
- Configuring injection productions in the event ledger
- Submitting and monitoring jobs through Asimov
- Understanding the workflow and outputs
- Chaining Minke with downstream analysis productions
- Troubleshooting common issues

Using Minke with Asimov provides a robust, reproducible way to generate injection frames as part of a larger gravitational wave data analysis workflow. This is particularly useful for:

- Systematic injection studies
- Algorithm validation
- Detector characterization
- Sensitivity estimates
- Parameter estimation validation

For more information, see:

- `Asimov Documentation <https://asimov.docs.ligo.org/>`_
- `Minke Source Code <https://github.com/transientlunatic/minke>`_
- :doc:`software-injections` for more injection examples
- :doc:`tutorial-noise-psd` for details on noise generation
