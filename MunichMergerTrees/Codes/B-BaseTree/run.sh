#!/bin/bash

export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/users/home/boyd/local/lib
for i in {0..47}
do
./B-BaseTree B-BaseTree.param $( printf '%03d' $i)
done
