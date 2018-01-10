#!/bin/bash

#this script requires wget, python 2.7, virtualenv and pip
#OS-X has trouble with the gpu package of lammps, install XCode8 (not 9) and CUDA

#commented out lines allow for the building of a gpu accelerated version. I will clean this process up at some point...

WDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

while test $# -gt 0; do
        case "$1" in
                -h|--help)
						echo "EMERGENT OPTIMISER INSTALLER"
						echo " "
                        echo "A script for installing the emergent optimiser package. This script will install lammps as a shared local library, and create a python virtual environment with the required pip modules."
                        echo " "
                        echo "You will need bash (with wget), pip, python, and virtualenv to install the program."
                        echo " "
                        echo "Once installed, use activate.sh to enter the virtual environment, and refer to the help in the python runscript for further direction."
                        echo " "
                        echo "options:"
                        echo "-h, --help                show brief help"
                        exit 0
                        ;;
                *)
                        break
                        ;;
        esac
done



CUDA_HOME=/usr/local/cuda
if [ ! -d ${WDIR}/lammps ]; then
	wget -qO- http://lammps.sandia.gov/tars/lammps-stable.tar.gz | tar xvz 
	mv lammps* lammps
fi
#copy src files to lammps folder
cd ${WDIR}/lammps
LAMMPSDIR=$(pwd)
if [ ! -f ${LAMMPSDIR}/src/liblammps.so ]; then
	cd ${LAMMPSDIR}/src
	make clean-all
	make no-all
	make yes-dipole yes-rigid yes-molecule yes-python yes-opt #yes-gpu
	#cd ${LAMMPSDIR}/lib/gpu && make -f Makefile.mpi CUDA_LIB="-L${CUDA_HOME}/lib"
	cp -rf ${WDIR}/src/*.h ${LAMMPSDIR}/src
	cp -rf ${WDIR}/src/*.cpp ${LAMMPSDIR}/src
	make -j4 mpi #LMP_INC="-DLAMMPS_PNG -DLAMMPS_JPEG -DLAMMPS_FFMPEG -DLAMMPS_EXCEPTIONS" JPG_LIB="-lpng -ljpeg" gpu_SYSPATH="-L${CUDA_HOME}/lib"
	make -j4 mpi mode=shlib #LMP_INC="-DLAMMPS_PNG -DLAMMPS_JPEG -DLAMMPS_FFMPEG -DLAMMPS_EXCEPTIONS" JPG_LIB="-lpng -ljpeg" gpu_SYSPATH="-L${CUDA_HOME}/lib"
	make install-python
fi
cd ${WDIR}
if [ ! -d ${WDIR}/venv ]; then
	pip install virtualenv
	virtualenv -p `which python2.7` venv
fi
source ./activate.sh
source ${WDIR}/pipinstall.sh
cd ${LAMMPSDIR}/python
python install.py
cd ${WDIR}

echo "done!"
