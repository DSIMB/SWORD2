#!/bin/bash

# Compile DSSP dependency
if [ "$(uname)" == "Darwin" ]; then
    ./bin/Dssp/DsspCompileGCCmacos &>/dev/null
else
    ./bin/Dssp/DsspCompileGCC &>/dev/null
fi
mv dsspcmbi bin/Dssp/dsspcmbi &>/dev/null

if [ -f ./bin/Dssp/dsspcmbi ]
then
    echo "Successfully compiled DSSP dependency"
else
    echo "Error: unable to compile DSSP dependency"
    exit 1
fi

# Compile MyPMFs
make -C bin/mypmfs-master scoring_omp >/dev/null
if [ -f bin/mypmfs-master/scoring_omp ]
then
    echo "Successfully compiled MyPMFs"
else
    echo "Error: unable to compile MyPMFs"
    exit 1
fi

# Compile Peeling
make -C bin/Peeling >/dev/null
if [ -f bin/Peeling/Peeling_omp ]
then
    echo "Successfully compiled Peeling"
else
    echo "Error: unable to compile Peeling"
    exit 1
fi
