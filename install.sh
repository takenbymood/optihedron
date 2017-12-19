#!/bin/bash

#this script requires wget, boost, nvcc, python 2.7 and pip
#OS-X has trouble with the gpu package of lammps, install XCode8 (not 9) and CUDA

#commented out lines allow for the building of a gpu accelerated version. I will clean this process up at some point...


wdir=$(pwd)
CUDA_HOME=/usr/local/cuda
if [ ! -d ${wdir}/lammps ]; then
	wget -qO- http://lammps.sandia.gov/tars/lammps-stable.tar.gz | tar xvz 
	mv lammps* lammps
fi
cd lammps
LAMMPSDIR=$(pwd)
if [ ! -f ${LAMMPSDIR}/src/liblammps.so ]; then
	cd src
	make clean-all
	make no-all
	make yes-rigid yes-molecule yes-python yes-opt #yes-gpu
	#cd ${LAMMPSDIR}/lib/gpu && make -f Makefile.mpi CUDA_LIB="-L${CUDA_HOME}/lib"
	cd ${LAMMPSDIR}/src
	make mpi #LMP_INC="-DLAMMPS_PNG -DLAMMPS_JPEG -DLAMMPS_FFMPEG -DLAMMPS_EXCEPTIONS" JPG_LIB="-lpng -ljpeg" gpu_SYSPATH="-L${CUDA_HOME}/lib"
	make mpi mode=shlib #LMP_INC="-DLAMMPS_PNG -DLAMMPS_JPEG -DLAMMPS_FFMPEG -DLAMMPS_EXCEPTIONS" JPG_LIB="-lpng -ljpeg" gpu_SYSPATH="-L${CUDA_HOME}/lib"
fi
if [ ! -d ${wdir}/venv ]; then
	pip install virtualenv
	cd $wdir
	virtualenv venv
fi
cd ${wdir}
source ${wdir}/activate.sh
source ${wdir}/pipinstall.sh
cd ${LAMMPSDIR}/python
python install.py
cd ${wdir}

echo "done!"
