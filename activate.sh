#!/bin/bash


TK_LIBRARY=/usr/lib/python2.7/lib-tk:/usr/lib/python2.7/site-packages/PIL:/usr/lib
TKPATH=/usr/lib/python2.7/lib-tk:/usr/lib/python2.7/site-packages/PIL:/usr/lib 
TCL_LIBRARY=/usr/lib 
export TCL_LIBRARY TK_LIBRARY TKPATH
export DYLD_FALLBACK_LIBRARY_PATH=$(pwd)/lammps/src:$DYLD_FALLBACK_LIBRARY_PATH
export TMPDIR=~/tmp


source venv/bin/activate