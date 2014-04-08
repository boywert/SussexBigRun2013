################################
#$ -S /bin/bash
#$ -j y 
#$ -cwd
#$ -m n
### must be a multiple of 12 up to 48
#$ -pe impi4 36
### up to 24 hours
#$ -l h_rt=24:00:00
#$ -M bhenriques@mpa-garching.mpg.de
module load impi
ulimit -aH
#export I_MPI_DEBUG=2

#mpiexec -np $NSLOTS ./L-Galaxies ./input/input_MR_W1_W7.par > prog01.out

#time mpiexec -n $NSLOTS ./L-Galaxies ./input/input_guo11_MR_W1_W1.par > prog01.out

time mpiexec -n $NSLOTS ./L-Galaxies ./input/input_MRII_W1_PLANCK.par > prog02.out 
#time mpiexec -n $NSLOTS ./L-Galaxies ./input/input_guo11_MRII_W1_W1.par > prog02.out 
###time mpiexec -n $NSLOTS ./L-Galaxies ./input/input_mcmc_MR_plus_MRII_W1_planck.par > prog11.out
date
####$TMPDIR
################################


