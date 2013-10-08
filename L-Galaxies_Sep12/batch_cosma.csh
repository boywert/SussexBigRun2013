#!/bin/tcsh
#
#BSUB -J L-Galaxies[1-16]
#BSUB -q cosma
#BSUB -P durham
##BSUB -P dp004
#BSUB -n 12
#BSUB -oo logs/batch_500_%I.out
limit stacksize unlimited
limit descriptors unlimited
newgrp millgas
umask 002

@ i = $LSB_JOBINDEX - 1
@ ff = $i
@ lf = $i

# Load modules
source /etc/profile.d/modules.csh
module load intel_comp/2012.0.032
module load platform_mpi/intel_2012.0.032/8.2.0
module load fftw/intel_2012.0.032_platform_mpi_8.2.0/2.1.5
module load gsl/gnu_4.4.6/1.15

# Set working directory
#cd $BSUB_O_WORKDIR

# Create personalised input parameter files
echo FirstFile $ff | cat >  input_batch/input.500_$i
echo LastFile $lf | cat >>  input_batch/input.500_$i
echo MaxMemSize 1500 >> input_batch/input.500_$i
cat input/input.500 >> input_batch/input.500_$i

# Run jobs
./L-Galaxies input_batch/input.500_$i
