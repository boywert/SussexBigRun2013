#!/bin/tcsh

#for ((i=0; i<512; i+=1 ))
for i in 47 48 49
do

# Wrapper to ensure that we create new files only if not already processed
if [ ! -e output/SA_z0.00_$i ]
then

echo File $i does not exist

else

echo File $i already processed

fi
done
