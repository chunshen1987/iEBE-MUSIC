#!/usr/bin/env bash

parafile=$1
processId=$2
nHydroEvents=$3
nUrQMD=$4
nthreads=$5
seed=$6
bayesFile=$7


# Run the singularity container
export PYTHONIOENCODING=utf-8
export PATH="${PATH}:/usr/lib64/openmpi/bin:/usr/local/gsl/2.5/x86_64/bin"
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:/usr/local/lib:/usr/local/gsl/2.5/x86_64/lib64"

printf "Start time: `/bin/date`\n"
printf "Job is running on node: `/bin/hostname`\n"
printf "system kernel: `uname -r`\n"
printf "Job running as user: `/usr/bin/id`\n"

if [ -z ${bayesFile} ]
then
    /home/iEBE-MUSIC/generate_jobs.py -w playground -c wsugrid -par ${parafile} -id ${processId} -n_hydro ${nHydroEvents} -n_th ${nthreads} -n_urqmd ${nUrQMD} -seed ${seed} --nocopy --continueFlag
else
    /home/iEBE-MUSIC/generate_jobs.py -w playground -c wsugrid -par ${parafile} -id ${processId} -n_hydro ${nHydroEvents} -n_th ${nthreads} -n_urqmd ${nUrQMD} -seed ${seed} -b ${bayesFile} --nocopy --continueFlag
fi


cd playground/event_0
bash submit_job.script
status=$?
if [ $status -ne 0 ]; then
    exit $status
fi
