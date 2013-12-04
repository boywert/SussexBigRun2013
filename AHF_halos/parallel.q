#!/bin/bash
#$ -N boyd
#$ -M cs390@sussex.ac.uk
#$ -m bea
#$ -j y
#$ -cwd
#$ -pe openmpi 8 #eg12-36
#$ -q mps_amd.q
#$ -S /bin/bash
# source modules environment:
module add sge

while read line
do
    redshift=$(printf '%3.3f' $line)
    echo $redshift
done < halofinds

##mpirun -np 8 ../bin/AHF-v1.0-056 AHF.input-template2

