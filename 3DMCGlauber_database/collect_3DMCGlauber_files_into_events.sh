#!/usr/bin/env bash

nev=100
event_folder=$1

cd $event_folder
for (( i=0; i < $nev; i++ ))
do
    mkdir event-$i
    mv strings_event_$i.dat event-$i
done
cd -
