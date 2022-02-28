#!/usr/bin/env python3
"""This script generates the job submission script on OSG"""


import sys
from os import path
import random

FILENAME = "singularity.submit"

def print_usage():
    """This function prints out help messages"""
    print("Usage: {} ".format(sys.argv[0].split("/")[-1])
          + "Njobs Nevents_per_job N_threads SingularityImage ParameterFile "
          + "jobId [bayesFile]")


def write_submission_script(para_dict_):
    jobName = "iEBEMUSIC_{}".format(para_dict_["job_id"])
    random_seed = random.SystemRandom().randint(0, 10000000)
    imagePathHeader = "stash:///osgconnect"
    script = open(FILENAME, "w")
    if para_dict_["bayesFlag"]:
        script.write("""universe = vanilla
executable = run_singularity.sh
arguments = {0} $(Process) {1} {2} {3} {4}
""".format(para_dict_["paraFile"], para_dict_["n_events_per_job"],
           para_dict_["n_threads"], random_seed, para_dict_["bayesFile"]))
    else:
        script.write("""universe = vanilla
executable = run_singularity.sh
arguments = {0} $(Process) {1} {2} {3}
""".format(para_dict_["paraFile"], para_dict_["n_events_per_job"],
           para_dict_["n_threads"], random_seed))
    script.write("""
JobBatchName = {0}

should_transfer_files = YES
WhenToTransferOutput = ON_EXIT

+SingularityImage = "{1}"
Requirements = SINGULARITY_CAN_USE_SIF && StringListIMember("stash", HasFileTransferPluginMethods)
""".format(jobName, imagePathHeader + para_dict_["image_with_path"]))

    if para_dict_['bayesFlag']:
        script.write("""
transfer_input_files = {0}, {1}
""".format(para_dict_['paraFile'], para_dict_['bayesFile']))
    else:
        script.write("""
transfer_input_files = {0}
""".format(para_dict_['paraFile']))

    script.write(
            "transfer_checkpoint_files = EVENT_RESULTS_$(Process).tar.gz\n")

    script.write("""
transfer_output_files = playground/event_0/EVENT_RESULTS_$(Process)/spvn_results_$(Process).h5

error = ../log/job.$(Cluster).$(Process).error
output = ../log/job.$(Cluster).$(Process).output
log = ../log/job.$(Cluster).$(Process).log

+JobDurationCategory = "Long"

# remove the failed jobs
periodic_remove = (ExitCode == 73)

checkpoint_exit_code = 85

# Send the job to Held state on failure.
on_exit_hold = (ExitBySignal == True) || (ExitCode != 0 && ExitCode != 73)

# The below are good base requirements for first testing jobs on OSG,
# if you don't have a good idea of memory and disk usage.
request_cpus = {0:d}
request_memory = 2 GB
request_disk = 1 GB

# Queue one job with the above specifications.
queue {1:d}""".format(para_dict_["n_threads"], para_dict_["n_jobs"]))
    script.close()


def write_job_running_script(para_dict_):
    script = open("run_singularity.sh", "w")
    script.write("""#!/usr/bin/env bash

parafile=$1
processId=$2
nHydroEvents=$3
nthreads=$4
seed=$5

# Run the singularity container
export PYTHONIOENCODING=utf-8
export PATH="${PATH}:/usr/lib64/openmpi/bin:/usr/local/gsl/2.5/x86_64/bin"
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:/usr/local/lib:/usr/local/gsl/2.5/x86_64/lib64"

printf "Start time: `/bin/date`\\n"
printf "Job is running on node: `/bin/hostname`\\n"
printf "system kernel: `uname -r`\\n"
printf "Job running as user: `/usr/bin/id`\\n"

""")
    if para_dict_["bayesFlag"]:
        script.write("""bayesFile=$6

/home/iEBE-MUSIC/generate_jobs.py -w playground -c OSG -par ${parafile} -id ${processId} -n_th ${nthreads} -n_urqmd ${nthreads} -n_hydro ${nHydroEvents} -seed ${seed} -b ${bayesFile} --nocopy --continueFlag
""")
    else:
        script.write("""
/home/iEBE-MUSIC/generate_jobs.py -w playground -c OSG -par ${parafile} -id ${processId} -n_th ${nthreads} -n_urqmd ${nthreads} -n_hydro ${nHydroEvents} -seed ${seed} --nocopy --continueFlag
""")

    script.write("""
cd playground/event_0
bash submit_job.pbs
status=$?
if [ $status -ne 0 ]; then
    exit $status
fi
""")
    script.close()


def main(para_dict_):
    write_submission_script(para_dict_)
    write_job_running_script(para_dict_)


if __name__ == "__main__":
    bayesFlag = False
    bayesFile = ""
    try:
        N_JOBS = int(sys.argv[1])
        N_EVENTS_PER_JOBS = int(sys.argv[2])
        N_THREADS = int(sys.argv[3])
        SINGULARITY_IMAGE_PATH = sys.argv[4]
        SINGULARITY_IMAGE = SINGULARITY_IMAGE_PATH.split("/")[-1]
        PARAMFILE = sys.argv[5]
        JOBID = sys.argv[6]
        if len(sys.argv) == 8:
            bayesFile = sys.argv[7]
            bayesFlag = True
    except (IndexError, ValueError) as e:
        print_usage()
        exit(0)

    para_dict = {
        'n_jobs': N_JOBS,
        'n_events_per_job': N_EVENTS_PER_JOBS,
        'n_threads': N_THREADS,
        'image_name': SINGULARITY_IMAGE,
        'image_with_path': SINGULARITY_IMAGE_PATH,
        'paraFile': PARAMFILE,
        'job_id': JOBID,
        'bayesFlag': bayesFlag,
        'bayesFile': bayesFile,
    }

    main(para_dict)

