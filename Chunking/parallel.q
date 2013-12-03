#!/bin/sh
#$ -N boyd
#$ -M cs390@sussex.ac.uk
#$ -m bea
#$ -j y
#$ -cwd
#$ -pe openmpi 216 #eg12-36
#$ -q mps_amd.q
#$ -S /bin/bash
# source modules environment:
module add sge

mpirun -np $NSLOTS ./chunk chunking_parameters.dat
