#!/usr/bin/env bash

usage="$0 fromFolder toFolder"

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
cp ${fromFolder}/parameter* ${target_folder}/
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
    evId=`echo $ijob | cut -f 2 -d _`
    if [ ! -f ${eventsPath}/RESULTS_${evId}.h5 ]; then
        (cd ${eventsPath}
        mkdir -p temp
        simFolder=`ls | grep "playground"`
        ./collect_events.sh ${simFolder} temp
        mv temp/${simFolder}/${simFolder}.h5 RESULTS_${evId}.h5
        )
    fi
    mv ${eventsPath}/RESULTS_${evId}.h5 $target_spvn_folder
    mv ${eventsPath}/temp/playground/HYDRO_RESULTS/${hydro_folder_name}* $target_hydro_folder 2>/dev/null
    mv ${eventsPath}/temp/playground/URQMD_RESULTS/${UrQMD_file_name}* $target_urqmd_folder 2>/dev/null
done

if [ -f ${target_folder}/${folderName}.h5 ]; then
    mv ${target_folder}/${folderName}.h5 ${target_spvn_folder}
fi
./combine_multiple_hdf5.py ${target_spvn_folder}
mv SPVN_RESULTS.h5 ${target_folder}/${folderName}.h5
rm -fr $target_spvn_folder
