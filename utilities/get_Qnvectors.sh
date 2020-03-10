#!/usr/bin/env bash

evfolder=${1%"/"}
database=${2%"/"}

for iev in `ls --color=none ${evfolder}`
do
    eventid=`echo $iev | sed 's/hydro_results_//' | sed 's$/$$'`
    ./fetch_Qnvectors_from_hdf5_database.py $database $eventid
    mv Qn_vectors_${eventid}.dat ${evfolder}/$iev/
done
