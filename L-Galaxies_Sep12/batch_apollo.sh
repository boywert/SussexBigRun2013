#!/bin/bash
# Tell SGE that we are using the bash shell
#$ -S /bin/bash

# Example file to create a task-farm of identical jobs on Apollo

# Do not lmit the stacksize (the maximum memory a job can use)
ulimit -s unlimited
# Do not limit the number of open files a job can have
ulimit -n unlimited

# Run the job from the following directory
#cd /mnt/lustre/scratch/petert/L-Galaxies/L-Galaxies

# Created files will have the fs-virgo group
# This feature seems to be disabled on Apollo, so this does not work
newgrp fs-virgo
# Created files will have rw permission for the group and r for the world
umask 002

# Set pathnames below relative to the current working directory
#$ -cwd
# Say which queue you want to submit to
#$ -q mps_amd.q
# Define a task farm of jobs
#$ -t 216
# Limit to 50 concurrent jobs
#$ -tc 216
# Join standard error to standard out
#$ -j y
# Give the job a name
#$ -N Lgalaxy
# Name and location of the output file
# SGE will only substitute certain variables here
#$ -o logs/$JOB_NAME_$TASK_ID.log

# The parentheses here allow one to do algebra with shell variables
i=$(($SGE_TASK_ID -1))
echo Running on file $i


# Run jobs
./L-Galaxies input_batch/input.500_$i


inputfolder="input_cubep3m/"
template="input_cubep3m/input_template"
exec=./L-Galaxies


filename="input_cubep3m/input_${i}"
echo "FirstFile" $i > $filename
echo "LastFile" $i >> $filename
cat $template >> $filename

$exec $filename