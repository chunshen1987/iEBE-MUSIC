#!/usr/bin/env bash

usage="./check_running_jobs.sh FolderName"

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
        output=`qstat -x -n1 ${job_idstr} | tail -n 1`
        job_status=`echo ${output} | awk {'print $10'}`
    fi
    if [[ ${job_status} == "F" ]]
    then
        echo "Job finished, restarting ..."
        (
            cd ${eventsPath}
            qsub -q wsuq submit_job.pbs > job_id
        )
    else
        echo ${output}
    fi
done

