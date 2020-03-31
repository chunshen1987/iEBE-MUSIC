#!/usr/bin/env bash

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
    queue="wsuq"
fi

echo "submit jobs in " ${workFolder}
echo ${workFolder} >> current_sims_list.txt

Numjobs=0
cd ${workFolder}
for ijob in `ls --color=none | grep "event"`;
do
    echo "submit job in " ${workFolder}/${ijob}
    cd ${ijob}
    qsub -q ${queue} submit_job.pbs > job_id
    cd ..
    ((Numjobs++))
done

echo "Submitted " ${Numjobs} " jobs in total. Have a nice day!"
