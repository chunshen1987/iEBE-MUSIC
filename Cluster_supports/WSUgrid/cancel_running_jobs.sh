#!/usr/bin/env bash

usage="./cancel_running_jobs.sh fromFolder"

jobFolder=$1

if [ -z "$jobFolder" ]
then
    echo $usage
    exit 1
fi

jobFolder=${jobFolder%"/"}

echo "cancelling jobs in " ${jobFolder}

for ijob in `ls --color=none ${jobFolder} | grep "event" `;
do
    eventsPath=${jobFolder}/${ijob}
    scancel `cat ${eventsPath}/job_id`
done

