#!/usr/bin/env bash

event_folder=$1
nev=0

(
    cd $event_folder
    for ifile in `ls | grep "strings_event"`
    do
        event_id=`echo $ifile | sed 's/.dat//' | cut -f 3 -d "_"`
        mkdir event-$event_id
        mv strings_event_$event_id.dat event-$event_id
        nev=$((nev + 1))
    done
    echo "There are $nev events in total."
)

