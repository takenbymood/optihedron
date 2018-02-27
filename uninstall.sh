#!/bin/bash

STARTDIR=$(pwd)
cd "$( dirname "${BASH_SOURCE[0]}" )"
WDIR=$(pwd)
LAMMPSDIR="./lammps"

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
						rm -rf $LAMMPSDIR
                        break
                        ;;
                *)
                        break
                        ;;
        esac
done

if [ -d $LAMMPSDIR ]; then
	if [ -f $LAMMPSDIR/src/liblammps.so ]; then
		cd $LAMMPSDIR/src
		rm liblammps.so
		rm lmp_mpi
		rm lmp_serial
		rm lammps
		make clean-all
		cd "$WDIR"
	fi
fi
VENVDIR="./venv"
if [ -d $VENVDIR ]; then
	rm -rf $VENVDIR
fi

cd "${STARTDIR}"
echo "done!"
