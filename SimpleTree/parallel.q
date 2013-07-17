#!/bin/sh
#$ -NJobName
#$ -M user@sussex.ac.uk
#$ -mbea
#$ -cwd
#$ -pe openmpi NUMBER_OF_CPUS #eg12-36
#$ -q parallel.q
#$ -S /bin/bash
# source modules environment:
./etc/profile.d/modules.sh
module add openmpi/intel/64/1.6.4
mpirun -np $NSLOTS -machinefile $TMPDIR/machines /path/to/exec
