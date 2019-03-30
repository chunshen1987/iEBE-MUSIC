#!/usr/bin/env bash

usage="./collect_events.sh fromFolder toFolder"

fromFolder=$1
toFolder=$2

if [ -z "$fromFolder" ]
then
    echo $usage
    exit 1
fi
if [ -z "$toFolder" ]
then
    echo $usage
    exit 1
fi

echo "collecting events from " $fromFolder " to " $toFolder

folderName=$fromFolder
target_folder=$toFolder/$folderName
mkdir -p $target_folder
target_hydro_folder=$target_folder/HYDRO_RESULTS
mkdir -p $target_hydro_folder
target_urqmd_folder=$target_folder/URQMD_RESULTS
mkdir -p $target_urqmd_folder
target_spvn_folder=$target_folder/SPVN_RESULTS
mkdir -p $target_spvn_folder

event_folder_name="EVENT_RESULTS_"
hydro_folder_name="hydro_results_"
UrQMD_file_name="particle_list_"
spvn_folder_name="spvn_results_"

total_eventNum=0
collected_eventNum=0
for ijob in `ls --color=none $fromFolder | grep "event" `;
do 
    eventsPath=$fromFolder/$ijob
    for iev in `ls --color=none $eventsPath | grep $event_folder_name`
    do 
        echo $iev
        event_id=`echo $iev | cut -f 3 -d "_"`
        hydrostatus=`tail -n 1 $eventsPath/$iev/$hydro_folder_name$event_id/run.log | cut -f 4 -d " "`
        echo $hydrostatus
        if [ "$hydrostatus" == "Finished." ]; then
            if [ -a $eventsPath/$iev/$spvn_folder_name$event_id/particle_9999_vndata_eta_-0.5_0.5.dat ]; then
                mv $eventsPath/$iev/$hydro_folder_name$event_id $target_hydro_folder
                mv $eventsPath/$iev/$UrQMD_file_name$event_id.gz $target_urqmd_folder
                #mv $eventsPath/$iev/$spvn_folder_name$event_id $target_spvn_folder
                #mv $target_hydro_folder/$hydro_folder_name$event_id/eccentricities_evo_eta_-0.5_0.5.dat $target_spvn_folder/$spvn_folder_name$event_id
                #mv $target_hydro_folder/$hydro_folder_name$event_id/momentum_anisotropy_eta_-0.5_0.5.dat $target_spvn_folder/$spvn_folder_name$event_id
                mv $eventsPath/$iev/$spvn_folder_name$event_id.h5 $target_spvn_folder
                ((collected_eventNum++))
            fi
        fi
        ((total_eventNum++))
    done
done

echo "Collected events number: " $collected_eventNum " out of " $total_eventNum

./combine_results_into_hdf5.py $target_spvn_folder
mv SPVN_RESULTS.h5 $target_folder/$folderName.h5
rm -fr $target_spvn_folder
