#!/bin/sh
#$ -NJobName
#$ -M cs390@sussex.ac.uk
#$ -mbea
#$ -cwd
#$ -pe openmpi 64 #eg12-36
#$ -q mps_amd.q
#$ -S /bin/bash
# source modules environment:
module add sge

mpirun -np $NSLOTS -machinefile $TMPDIR/machines /path/to/exec
