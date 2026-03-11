#!/bin/bash

# DSSP is now implemented in pure Rust — no C compilation needed.

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
