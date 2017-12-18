#!/bin/bash

export DYLD_FALLBACK_LIBRARY_PATH=$(pwd)/lammps/src:$DYLD_FALLBACK_LIBRARY_PATH
export TMPDIR=~/tmp


source venv/bin/activate