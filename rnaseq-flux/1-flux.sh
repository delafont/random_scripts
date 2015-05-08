#!/bin/bash

# Runs FluxSimulator with the given parameters file
# (which includes the path to the fasta)

if [ $# -lt 1 ]
then
    params='simulation.txt'
else
    params=$1
fi

# Simulation
job='flux'
bsub -q long -J $job -e logs/${job}.err -o logs/${job}.out " \
    flux-simulator -t simulator -x -l -s -p $params ;
"

