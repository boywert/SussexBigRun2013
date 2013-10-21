#!/bin/bash

#export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/pactusers/cs390/local/lib
for i in {30..30}
do
qsub jobscript $( printf '%03d' $i)
done
