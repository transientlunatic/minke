.. This is the Minke documentation file which describes the use of the
   pre-made docker containers which are maintained for Minke.

Using Minke with Docker
=======================

One of the easiest methods for producing MDC sets with Minke is using
the pre-built Minke containers which can be downloaded and run on any
system which can run Docker; this includes Linux, Mac, and Windows.

A Docker container contains all of the software prerequisites for the
tool it contains; the Minke container comes with all of the required
LALSuite code pre-compiled, which increases the speed at which Minke
can be installed on most systems. The container is also entirely
self-contained, so updates to your underlying operating system's
software won't break the installation of Minke.

To start a terminal which has access to the latest version of Minke first run ::
  docker pull lpmn/minke:latest

which will download the container to your machine, and then  ::
  docker run -it lpmn/minke:latest /bin/bash -l

You can also directly run a python script which uses Minke ::
  docker run lpmn/minke:latest python script.py

but it's important to remember that containers are entirely
self-contained, so any files which are written out by a script won't
appear in the file system on the "host" machine. To get around this we
need to link the host file system to the file system in the container.

Connecting your filesystem
--------------------------

The host file system can be connected to the Minke container by adding
an additional flag on the docker run command, which links a directory
on the host machine to a mount point in the container. For example ::
  docker run lpmn/minke:latest -v /home/minke:/minke

would link a directory on the host machine (``home/minke``) to another
in the container (``/minke``). Files saved to ``/minke`` in the
container would be saved to ``/home/minke`` on the computer you ran
the ``docker run`` command on.

Using containers on HTCondor
----------------------------

If you're a member of the LVC there's a good chance you'll want to run
your Minke job on Condor (a job scheduling system for high-performance
computing). Docker is supported through a special *universe* which
runs the job inside a container; an example ``.sub`` file might look
something like this (assuming you've added a shebang pointing to
Python in your script, and made it executable). ::

  universe = docker
  docker_image = lpmn/minke:latest
  executable = /home/minke/script.py
  arguments = sg 10 10
  output = docker_job.out
  error = docker_job.err
  queue

Like before the output of any script will be stored in the container
(which may be desirable for allowing files to be moved between the
compute and head nodes in an efficient manner) but if you want to be
able to read or write to a scratch directory we need to connect the
container to the host machine's file system. We need to add some extra
instructions to the sub file to move data back out of the
container. ::

    should_transfer_files = YES
    when_to_transfer_output = ON_EXIT

