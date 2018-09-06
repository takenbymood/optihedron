#!/bin/bash

#PBS -N subhedra
#PBS -q batch
#PBS -l nodes=1:ppn=_CORES_
#PBS -l pmem=800M
#PBS -l walltime=240:00:00


#PBS -V
_MODULES_
cd _DIR_/..
source _DIR_/activate.sh
_PRE_ _DIR_/venv/bin/python _DIR_/plammps.py -s _SCRIPT_
deactivate
cd ..