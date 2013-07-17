#!/bin/sh
#$ -N Boyd
#$ -S /bin/sh
# Reference file IO (.e and .o files) to the
# current working directory
#$ -cwd
#Merge the standard out and standard error to one file
#$ -j y
#Use the serial.q
#$ -q serial.q
# If modules are needed, source modules environment:
. /etc/profile.d/modules.sh
# Add any modules you might require:
module add sge
module add gcc
module add shared
#
$1

