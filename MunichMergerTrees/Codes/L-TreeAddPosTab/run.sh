#!/bin/bash

export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/pactusers/cs390/local/lib
for i in {0..61}
do
./L-TreeAddPosTab L-TreeAddPosTab.param $( printf '%03d' $i)
done
