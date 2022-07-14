#!/usr/bin/env bash

usage="./check_and_restart_running_jobs.sh FolderName"

jobFolder=$1

if [ -z "$jobFolder" ]
then
    echo $usage
    exit 1
fi

jobFolder=${jobFolder%"/"}

echo "checking jobs in " ${jobFolder}

for ijob in `ls --color=none ${jobFolder} | grep "event" `;
do
    eventsPath=${jobFolder}/${ijob}
    job_idstr=`cat ${eventsPath}/job_id`
    if [[ ${job_idstr} == "" ]]
    then
        job_status="F"
    else
        output=`squeue --job ${job_idstr} | tail -n 1`
        if [[ ${output} == "" ]]
        then
            if [ -f "${eventsPath}/EVENT_RESULTS_*/*.h5" ]
            then
                job_status="S"
            else
                job_status="F"
            fi
        fi
    fi
    if [[ ${job_status} == "F" ]]
    then
        echo "Job finished, restarting ..."
        (
            cd ${eventsPath}
            sbatch -q requeue submit_job.script | awk {'print $4'} > job_id
        )
    else
        echo ${output}
    fi
done

