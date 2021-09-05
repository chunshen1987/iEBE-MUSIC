#!/usr/bin/env bash

usage="./collect_events_onlyh5.sh fromFolder toFolder"

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

fromFolder=${fromFolder%"/"}
toFolder=${toFolder%"/"}

echo "collecting events from " $fromFolder " to " $toFolder

folderName=`echo $fromFolder | rev | cut -d "/" -f 1 | rev`
target_folder=${toFolder}/${folderName}
mkdir -p ${target_folder}
cp ${fromFolder}/parameters_dict_*.py ${target_folder}/
target_hydro_folder=${target_folder}/HYDRO_RESULTS
mkdir -p ${target_hydro_folder}
target_urqmd_folder=${target_folder}/URQMD_RESULTS
mkdir -p ${target_urqmd_folder}
target_spvn_folder=${target_folder}/SPVN_RESULTS
mkdir -p ${target_spvn_folder}

hydro_folder_name="hydro_results_"
UrQMD_file_name="particle_list_"
spvn_folder_name="spvn_results_"

for ijob in `ls --color=none $fromFolder | grep "event" `;
do
    eventsPath=${fromFolder}/${ijob}
    for iev in `ls --color=none ${eventsPath} | grep "RESULTS"`
    do
        mv ${eventsPath}/${iev} $target_spvn_folder
        mv ${eventsPath}/temp/playground/HYDRO_RESULTS/${hydro_folder_name}* $target_hydro_folder
        mv ${eventsPath}/temp/playground/URQMD_RESULTS/${UrQMD_file_name}* $target_urqmd_folder
    done
done

if [ -f ${target_folder}/${folderName}.h5 ]; then
    mv ${target_folder}/${folderName}.h5 ${target_spvn_folder}
fi
./combine_multiple_hdf5.py ${target_spvn_folder}
mv SPVN_RESULTS.h5 ${target_folder}/${folderName}.h5
rm -fr $target_spvn_folder
