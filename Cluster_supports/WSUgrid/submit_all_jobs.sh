#!/usr/bin/env bash

usage="./submite_all_jobs.sh workFolder"

workFolder=$1

if [ -z "$workFolder" ]
then
    echo $usage
    exit 1
fi

echo "submit jobs in " $workFolder
echo $workFolder >> current_sims_list.txt

Numjobs=0
cd $workFolder
for ijob in `ls --color=none | grep "event"`;
do 
    echo "submit job in " $workFolder/$ijob
    cd $ijob
    qsub submit_job.pbs > job_id
    cd ..
    ((Numjobs++))
done

echo "Submitted " $Numjobs " jobs in total. Have a nice day!"
