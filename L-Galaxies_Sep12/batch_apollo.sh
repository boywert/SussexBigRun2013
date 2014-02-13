#!/bin/bash

firstfile=0
lastfile=215
inputfolder="input_cubep3m/"
template="input_cubep3m/input_template"
exec="./L-Galaxies"
for i in $(seq $firstfile $lastfile)
do
    filename="input_cubep3m/input_${i}"
    echo "FirstFile" $i > $filename
    echo "LastFile" $i >> $filename
    cat $template >> $filename
    qsub -N input_${i} -q mps_amd.q -cwd -j y $exec $filename
done