#!/usr/bin/env bash
# available queues are: primary, requeue, express, debug, secondary, gpu

usage="./submite_all_jobs.sh workFolder [queue_name]"

workFolder=$1
queue=$2

if [ -z "$workFolder" ]
then
    echo $usage
    exit 1
fi

if [ -z "$queue" ]
then
    queue="primary"
fi

echo "submit jobs in " ${workFolder}
echo ${workFolder} >> current_sims_list.txt

Numjobs=0
cd ${workFolder}
for ijob in `ls --color=none | grep "event"`;
do
    echo "submit job in " ${workFolder}/${ijob}
    cd ${ijob}
    sbatch -q ${queue} submit_job.script | awk {'print $4'} > job_id
    cd ..
    ((Numjobs++))
done

echo "Submitted " ${Numjobs} " jobs in total. Have a nice day!"
