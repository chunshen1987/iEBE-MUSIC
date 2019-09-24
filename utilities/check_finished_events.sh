#!/usr/bin/env bash

usage="./check_finished_events.sh fromFolder"

fromFolder=$1

if [ -z "$fromFolder" ]
then
    echo $usage
    exit 1
fi

echo "checking events from " $fromFolder

folderName=$fromFolder

event_folder_name="EVENT_RESULTS_"
hydro_folder_name="hydro_results_"
spvn_folder_name="spvn_results_"

total_eventNum=0
collected_eventNum=0
for ijob in `ls --color=none $fromFolder | grep "event" `;
do 
    eventsPath=$fromFolder/$ijob
    for iev in `ls --color=none $eventsPath | grep $event_folder_name`
    do 
        event_id=`echo $iev | rev | cut -f 1 -d "_" | rev`
        hydrostatus="Finished."
        hydrotime="-1"
        if [ -f "$eventsPath/$iev/$hydro_folder_name*$event_id/run.log" ]; then
            hydrotime=`tail -n 3 $eventsPath/$iev/$hydro_folder_name*$event_id/run.log 2>/dev/null | head -n 1 | cut -f 8 -d " "`
            hydrostatus=`tail -n 1 $eventsPath/$iev/$hydro_folder_name*$event_id/run.log 2>/dev/null | cut -f 4 -d " "`
        fi
        if [ "$hydrostatus" == "Finished." ]; then
            if [ -a $eventsPath/$iev/$spvn_folder_name*$event_id/particle_9999_vndata_eta_-2_2.dat ]; then
                echo $iev $hydrostatus $hydrotime
                ((collected_eventNum++))
            fi
        fi
        ((total_eventNum++))
    done
done

echo "Finished events number: " $collected_eventNum " out of " $total_eventNum
