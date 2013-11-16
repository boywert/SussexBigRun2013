#!/bin/sh
#$ -N boyd
#$ -M cs390@sussex.ac.uk
#$ -m bea
#$ -j y
#$ -cwd
#$ -pe openmpi 343 #eg12-36
#$ -q mps_amd.q
#$ -S /bin/bash
# source modules environment:
module add sge

mpirun -np $NSLOTS $1
