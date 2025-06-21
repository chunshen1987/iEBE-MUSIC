#!/usr/bin/env bash

evfolder=${1%"/"}
database=${2%"/"}

for iev in `ls --color=none ${evfolder}`
do
    eventid=`echo $iev | sed 's/hydro_results_//' | sed 's$/$$'`
    python3 fetch_Qnvectors_from_hdf5_database.py $database $eventid
    mv Qn_vectors_${eventid}.dat ${evfolder}/$iev/
    mv particle_yield_and_meanpT_${eventid}.dat ${evfolder}/$iev/
    mv Ncoll* ${evfolder}/$iev/ >> /dev/null 2>&1
done
