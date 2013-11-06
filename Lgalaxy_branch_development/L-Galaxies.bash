################################
#$ -S /bin/bash
#$ -j y
#$ -cwd
#$ -m n
### must be a multiple of 12 up to 48
#$ -pe impi4 60
### up to 24 hours
#$ -l h_rt=24:00:00
#$ -M bhenriques@mpa-garching.mpg.de
module load impi
ulimit -aH
#export I_MPI_DEBUG=2
time mpiexec -n $NSLOTS ./L-Galaxies ./input/input_MRII_W1_PLANCK.par > prog2.out
###time mpiexec -n $NSLOTS ./L-Galaxies ./input/input_mcmc_MR_plus_MRII_W1_planck.par > prog11.out
date
####$TMPDIR
################################

