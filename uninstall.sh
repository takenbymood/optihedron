#!/bin/bash

WDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

while test $# -gt 0; do
        case "$1" in
                -h|--help)
						echo "EMERGENT OPTIMISER UNINSTALLER"
						echo " "
                        echo "A script for uninstalling and cleaning the emergent optimiser package."
                        echo " "
                        echo "options:"
                        echo "-h, --help                show brief help"
                        echo "-r, --remove              remove ALL files and folders related to installation (rather than just removing binaries and libraries)"
                        exit 0
                        ;;
                -r|--remove)
						rm -rf ${WDIR}/lammps
                        break
                        ;;
                *)
                        break
                        ;;
        esac
done

if [ -d ${WDIR}/lammps ]; then
	cd ${WDIR}/lammps
	LAMMPSDIR=$(pwd)
	if [ -f ${LAMMPSDIR}/src/liblammps.so ]; then
		rm ${LAMMPSDIR}/src/liblammps.so
		rm ${LAMMPSDIR}/src/lmp_mpi
		rm ${LAMMPSDIR}/src/lmp_serial
		rm ${LAMMPSDIR}/src/lammps
		cd src
		make clean-all
		make no-all
	fi
fi
cd ${WDIR}
if [ -d ${WDIR}/venv ]; then
	rm -rf ${WDIR}/venv
fi