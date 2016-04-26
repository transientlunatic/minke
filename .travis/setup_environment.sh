#! /bin/bash
STARTLOC=$(pwd)
LIBFRAME_VER=v8r26
LALBRANCH=minke
mkdir .local
# CHECK IF WE'RE INSTALLING INTO A VIRTUALENV
if [ -n "$VIRTUAL_ENV" ]; then 
    LOCATION=$VIRTUAL_ENV
    #pip install numpy scipy
else 
    LOCATION=.local; 
fi



# DEFINE THE DIRECTORIES 
LALSUITE_PREFIX=${LOCATION}   # install directory: change as appropriate
LIBFRAME_PREFIX=${LOCATION}
METAIO_PREFIX=${LOCATION}


if [[ $1 = 'GRID' ]]; then
    echo $GRID
    echo "Installing on the LIGO Data Grid".
  LALSUITE_SRCDIR=/usr1/${USER}/repositories; # uncomment for LDG
elif [[ $1 = 'UWM' ]]; then
  LALSUITE_SRCDIR=/people/${USER}/repositories; # uncomment for Nemo, UWM
else
  LALSUITE_SRCDIR=${HOME}/repositories;
fi

# MAKE ALL OF THE REQUIRED DIRECTORIES
mkdir -p ${LALSUITE_SRCDIR}
mkdir -p ${LALSUITE_PREFIX}
mkdir -p ${LIBFRAME_PREFIX}

#
# LIBFRAME
#

# FETCH LIBFRAME
cd ${LALSUITE_SRCDIR}
wget http://lappweb.in2p3.fr/virgo/FrameL/$LIBFRAME_VER.tar.gz  # For the 2015-05-17 release
tar -xf $LIBFRAME_VER.tar.gz
cd $LIBFRAME_VER

# CONFIG BUILD INSTALL LIBFRAME
autoreconf --install
./configure --prefix=${LIBFRAME_PREFIX} --enable-swig-python
make                                   # (build results are stored on the src directory)
make install 
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${LIBFRAME_PREFIX} # to update the libraries paths
export PKG_CONFIG_PATH=${PKG_CONFIG_PATH}:${LIBFRAME_PREFIX}/lib/pkgconfig

#
# MetaIO
#

# Clone the git repo
cd ${LALSUITE_SRCDIR}
git clone git://versions.ligo.org/metaio.git

# Configure make and install
cd metaio

yes "\n" | ./00boot
yes "\n" | ./configure --enable-silent-rules --prefix=${METAIO_PREFIX} --enable-swig-python
make -j 2 --silent || make --silent
make install --silent
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${METAIO_PREFIX}/lib

#
# LALSUITE
#

# Clone the repo
cd ${LALSUITE_SRCDIR}
git clone git://versions.ligo.org/lalsuite.git # replace as appropriate

# Configure make and install
cd ${LALSUITE_SRCDIR}/lalsuite

# Checkout a specific branch of lalsuite if required
if [ -n "$LALBRANCH" ]; then 
    git checkout $LALBRANCH;
fi


./00boot
./configure --prefix=${LALSUITE_PREFIX} --enable-swig-python
make -j
make install

# Now ready the environment (this line needs to be run every time a shell session is started, so add it to the ~/.bashrc file
source ${LALSUITE_PREFIX}/etc/lalsuiterc

#
# PyLAL
#
cd ${LALSUITE_SRCDIR}/lalsuite/pylal
rm -rf build
python setup.py install --prefix=${LALSUITE_PREFIX}

#
# Glue
#
cd ${LALSUITE_SRCDIR}/lalsuite/glue
rm -rf build
python setup.py install --prefix=${LALSUITE_PREFIX}
source ${LALSUITE_PREFIX}/etc/glue-user-env.sh

cd ${STARTLOC}
