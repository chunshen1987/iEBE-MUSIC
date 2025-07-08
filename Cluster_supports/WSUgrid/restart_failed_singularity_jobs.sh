#!/usr/bin/env bash

usage="./restart_failed_singularity_jobs.sh workFolder [queue_name]"

runFolder=${1%/}
queue=$2

if [ -z "$runFolder" ]
then
    echo $usage
    exit 1
fi

if [ -z "$queue" ]
then
    queue="primary"
fi

(
cd $runFolder
for iev in `ls | grep "event"`
do
    if grep -q "FATAL" ${iev}/job.err; then
        jobId=`cat ${iev}/job_id`
        jobStatus=`squeue --job ${jobId} 2>&1`
        if [[ "$jobStatus" == *"error"* ]]; then
            echo "resubmit job in " ${runFolder}/${iev}
            cd $iev
            rm -fr playground*
            rm -fr temp*
            sbatch -q ${queue} submit_job.script | awk {'print $4'} > job_id
            cd ..
        fi
    fi
done
)
