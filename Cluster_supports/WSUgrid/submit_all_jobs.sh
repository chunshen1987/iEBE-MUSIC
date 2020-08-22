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
    queue="primary"
fi

partition=$SPRIMARY
if [ "$queue" == "requeue" ]; then
    partition=$SREQUEUE
fi
if [ "$queue" == "debug" ]; then
    partition=$SDEBUG
fi
if [ "$queue" == "secondary" ]; then
    partition=$SSECONDARY
fi
if [ "$queue" == "gpu" ]; then
    partition=$SGPU
fi
if [ "$queue" == "express" ]; then
    partition=$SEXPRESS
fi
echo $partition

echo "submit jobs in " ${workFolder}
echo ${workFolder} >> current_sims_list.txt

Numjobs=0
cd ${workFolder}
for ijob in `ls --color=none | grep "event"`;
do
    echo "submit job in " ${workFolder}/${ijob}
    cd ${ijob}
    sbatch -q ${queue} -p ${partition} submit_job.pbs | awk {'print $4'} > job_id
    cd ..
    ((Numjobs++))
done

echo "Submitted " ${Numjobs} " jobs in total. Have a nice day!"
