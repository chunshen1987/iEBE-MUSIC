#!/usr/bin/env bash

nev=10
event_folder=$1

cd $event_folder
for (( i=0; i < $nev; i++ ))
do
    mkdir event-$i
    mv epsilon-u-Hydro-t0.4-$i.dat event-$i
done
cd -
