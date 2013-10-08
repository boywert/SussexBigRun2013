#!/bin/bash
#
limit stacksize unlimited
limit descriptors unlimited
umask 002

for ((i=40;i<80;i++)) do

# Create personalised input parameter files
echo FirstFile $i | cat >  input_batch/input.$i
echo LastFile $i | cat >>  input_batch/input.$i
echo MaxMemSize 1000 >> input_batch/input.$i
#cat '/export/virgo/L-Galaxies_Sep12/input/input_MRII_Hen12.par' >> input_batch/input.$i
cat '/export/virgo/L-Galaxies_Sep12/input/input_MRII_Guo10.par' >> input_batch/input.$i

# Run job
./L-Galaxies input_batch/input.$i

done
