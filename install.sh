#!/bin/bash

# Run SWORD once to compile its dependency: DSSP
./bin/SWORD/bin/SWORD/SWORD &>/dev/null
if [ -f ./bin/SWORD/bin/SWORD/bin/Dssp/dsspcmbi ]
then
    echo "Successfully installed SWORD dependency"
else
    echo "Error: unable to compile necessary dependancies for SWORD"
    exit 1
fi

# Compile MyPMFs
make -C bin/mypmfs-master scoring >/dev/null
if [ -f bin/mypmfs-master/scoring ]
then
    echo "Successfully compiled MyPMFs"
else
    echo "Error: unable to compile necessary dependancies for SWORD"
    exit 1
fi