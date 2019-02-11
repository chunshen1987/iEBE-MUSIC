#!/usr/bin/env bash

usage="./submite_all_jobs.sh workFolder"

workFolder=$1

if [ -z "$workFolder" ]
then
    echo $usage
    exit 1
fi

echo "submit jobs in " $workFolder

Numjobs=0
cd $workFolder
for ijob in `ls --color=none`;
do 
    echo "submit job in " $workFolder/$ijob
    cd $ijob
    qsub submit_job.pbs
    cd ..
    ((Numjobs++))
done

echo "Submitted " $Numjobs " jobs in total. Have a nice day!"
