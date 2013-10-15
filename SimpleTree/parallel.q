#!/bin/sh
#$ -N boyd
#$ -M cs390@sussex.ac.uk
#$ -mbea
#$ -cwd
#$ -pe openmpi 32 #eg12-36
#$ -q mps_amd.q
#$ -S /bin/bash
# source modules environment:
module add sge

mpirun -np $NSLOTS -machinefile $TMPDIR/machines $1
