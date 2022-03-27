#!/usr/bin/env bash

run=${1%/}

runfoldername=`echo $run | rev | cut -d "/" -f 1 | rev`

mkdir -p $runfoldername
(
    cd $runfoldername
    cp $run/*.py ./
    cp $run/*.sh ./
    cp $run/*.submit ./
    ../combine_multiple_hdf5.py $run
)
