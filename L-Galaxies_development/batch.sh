#!/bin/bash
# Tell SGE that we are using the bash shell
#$ -S /bin/bash

# Example file to create a task-farm of identical jobs on Apollo

# Do not lmit the stacksize (the maximum memory a job can use)
#ulimit -s unlimited
# Do not limit the number of open files a job can have
#ulimit -n unlimited

# Run the job from the following directory
#cd /mnt/lustre/scratch/petert/L-Galaxies/L-Galaxies

# Created files will have the fs-virgo group
# This feature seems to be disabled on Apollo, so this does not work
#newgrp fs-virgo
# Created files will have rw permission for the group and r for the world
umask 002
mkdir -p logs
rm -f logs/*
# Set pathnames below relative to the current working directory
#$ -cwd
# Say which queue you want to submit to
#$ -q mps.q
#$ -l vf=32G
# Define a task farm of jobs
#$ -t 1-128
# Limit to 20 concurrent jobs
#$ -tc 128
# Join standard error to standard out
#$ -j y
# Give the job a name
#$ -N Lgalaxy_nore
# Name and location of the output file
# SGE will only substitute certain variables here
#$ -o logs/$JOB_NAME_$TASK_ID.log

# The parentheses here allow one to do algebra with shell variables
i=$(($SGE_TASK_ID -1))
echo Running on file $i


# Run jobs

inputfolder="inputs_47Mpc/"
mkdir -p $inputfolder
template="input/input_47mpc_template.par"
exec=./L-Galaxies
OutputDir="/mnt/lustre/scratch/cs390/47Mpc/outputs/no_reionization/"
SimulationDir="/mnt/lustre/scratch/cs390/47Mpc/"
ReionizationOn=0

mkdir -p $OutputDir
mkdir -p $inputfolder
filename="${inputfolder}/input_nr_${i}"


echo "FirstFile" $i > $filename
echo "LastFile" $i >> $filename
echo "OutputDir" $OutputDir >> $filename
echo "SimulationDir" $SimulationDir >> $filename
echo "ReionizationOn" $ReionizationOn >> $filename
cat $template >> $filename

if [ $i -eq 0 ]
then
  make metadata
  mkdir -p "${OutputDir}/inputs/"
  cp $filename "${OutputDir}/inputs/"
  cp python/LGalaxyStruct.py "${OutputDir}/inputs/"
fi


$exec $filename

