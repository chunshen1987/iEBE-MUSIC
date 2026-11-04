#!/usr/bin/env bash
# Runs inside the iEBE-MUSIC container (apptainer exec): generates the
# events of one job in workFolder with the code of the container and runs
# them.

workFolder=$1
parafile=$2
processId=$3
nHydroEvents=$4
nUrQMD=$5
nthreads=$6
seed=$7
bayesFile=$8

export PYTHONIOENCODING=utf-8

printf "Start time: `/bin/date`\n"
printf "Job is running on node: `/bin/hostname`\n"
printf "system kernel: `uname -r`\n"
printf "Job running as user: `/usr/bin/id`\n"

if [ -z ${bayesFile} ]
then
    /opt/iEBE-MUSIC/generate_jobs.py -w ${workFolder} -c local -par ${parafile} -id ${processId} -n_hydro ${nHydroEvents} -n_th ${nthreads} -n_urqmd ${nUrQMD} -seed ${seed} --nocopy --continueFlag
else
    /opt/iEBE-MUSIC/generate_jobs.py -w ${workFolder} -c local -par ${parafile} -id ${processId} -n_hydro ${nHydroEvents} -n_th ${nthreads} -n_urqmd ${nUrQMD} -seed ${seed} -b ${bayesFile} --nocopy --continueFlag
fi


cd ${workFolder}/event_0
bash submit_job.script
status=$?
if [ $status -ne 0 ]; then
    exit $status
fi
