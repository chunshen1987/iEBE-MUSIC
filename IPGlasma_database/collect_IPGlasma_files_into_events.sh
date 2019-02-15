#!/usr/bin/env bash

event_folder=$1

nev=0
tflag="t0.4"

(
    cd $event_folder
    
    for ifile in `ls | grep "epsilon-u-Hydro"`
    do
        event_id=`echo $ifile | sed 's/.dat//' | cut -f 5 -d "-"`
        if [[ -f epsilon-u-Hydro-$tflag-$event_id.dat ]]
        then
            mkdir event-$event_id
            mv epsilon-u-Hydro-$tflag-$event_id.dat event-$event_id
            nev=$((nev + 1))
        fi
    done
    echo "There are $nev events in total."
)
