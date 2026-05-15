#!/usr/bin/env python3
"""This script generates the job submission script on OSG"""


import sys
from os import path, makedirs
import argparse
import random

FILENAME = "singularity.submit"

def print_usage():
    """This function prints out help messages"""
    print("Usage: {} ".format(sys.argv[0].split("/")[-1])
          + "Njobs Nevents_per_job N_threads SingularityImage ParameterFile "
          + "jobId [bayes_file]")


def write_submission_script(para_dict_):
    jobName = "iEBEMUSIC_{}".format(para_dict_["job_name"])
    random_seed = random.SystemRandom().randint(0, 10000000)
    imagePathHeader = "osdf://"
    script = open(FILENAME, "w")
    if para_dict_["bayesFlag"]:
        script.write("""universe = vanilla
executable = run_singularity.sh
arguments = {0} $(Process) {1} {2} {3} {4}
""".format(para_dict_["param_file"], para_dict_["n_events_per_job"],
           para_dict_["n_threads"], random_seed, para_dict_["bayes_file"]))
    else:
        script.write("""universe = vanilla
executable = run_singularity.sh
arguments = {0} $(Process) {1} {2} {3}
""".format(para_dict_["param_file"], para_dict_["n_events_per_job"],
           para_dict_["n_threads"], random_seed))
    script.write("""
JobBatchName = {0}

should_transfer_files = YES
WhenToTransferOutput = ON_EXIT

+SingularityImage = "{1}"
Requirements = SINGULARITY_CAN_USE_SIF && StringListIMember("stash", HasFileTransferPluginMethods)
""".format(jobName, imagePathHeader + para_dict_["singularity_image_path"]))

    if para_dict_['bayesFlag']:
        script.write("""
transfer_input_files = {0}, {1}
""".format(para_dict_['param_file'], para_dict_['bayes_file']))
    else:
        script.write("""
transfer_input_files = {0}
""".format(para_dict_['param_file']))

    script.write(
            "transfer_checkpoint_files = playground/event_0/EVENT_RESULTS_$(Process).tar.gz\n")

    script.write("""
transfer_output_files = playground/event_0/EVENT_RESULTS_$(Process)/spvn_results_$(Process).h5

error = log/job.$(Cluster).$(Process).error
output = log/job.$(Cluster).$(Process).output
log = log/job.$(Cluster).$(Process).log

#+JobDurationCategory = "Long"
max_idle = 1000

# remove the failed jobs
periodic_remove = (ExitCode == 73)

# auto release hold jobs if they are caused by data transfer issues on OSG
periodic_release = ((HoldReasonCode == 13 || HoldReasonCode == 26) && (time() - EnteredCurrentStatus) > 1200 )

checkpoint_exit_code = 85

# Send the job to Held state on failure.
on_exit_hold = (ExitBySignal == True) || (ExitCode != 0 && ExitCode != 73)

# The below are good base requirements for first testing jobs on OSG,
# if you don't have a good idea of memory and disk usage.
request_cpus = {0:d}
request_memory = {1:d} GB
request_disk = 2 GB

# Queue one job with the above specifications.
queue {2:d}""".format(para_dict_["n_threads"], para_dict_["memory_per_job"],
                      para_dict_["n_jobs"]))
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

/opt/iEBE-MUSIC/generate_jobs.py -w playground -c OSG -par ${parafile} -id ${processId} -n_th ${nthreads} -n_urqmd ${nthreads} -n_hydro ${nHydroEvents} -seed ${seed} -b ${bayesFile} --nocopy --continueFlag
""")
    else:
        script.write("""
/opt/iEBE-MUSIC/generate_jobs.py -w playground -c OSG -par ${parafile} -id ${processId} -n_th ${nthreads} -n_urqmd ${nthreads} -n_hydro ${nHydroEvents} -seed ${seed} --nocopy --continueFlag
""")

    script.write("""
cd playground/event_0
mv EVENT_RESULTS_${processId}.tar.gz playground/event_0
bash submit_job.script
status=$?
if [ $status -ne 0 ]; then
    exit $status
fi
""")
    script.close()


def main(para_dict_):
    write_submission_script(para_dict_)
    write_job_running_script(para_dict_)
    logFolderName = "log"
    if not path.exists(logFolderName):
        makedirs(logFolderName)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Welcome to OSG script for the iEBE-MUSIC framework',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-n',
                        '--n_jobs',
                        metavar='',
                        type=int,
                        default=1,
                        help='number of jobs')
    parser.add_argument('-nev',
                        '--n_events_per_job',
                        metavar='',
                        type=int,
                        default=1,
                        help='number of events per job')
    parser.add_argument('-nth',
                        '--n_threads',
                        metavar='',
                        type=int,
                        default=1,
                        help='number of threads used for each job')
    parser.add_argument('-singularity',
                        '--singularity_image_path',
                        metavar='',
                        type=str,
                        default="",
                        help='singularity image path')
    parser.add_argument('-param',
                        '--param_file',
                        metavar='',
                        type=str,
                        default="",
                        help='parameter file')
    parser.add_argument('-jobid',
                        '--job_name',
                        metavar='',
                        type=str,
                        default="test",
                        help='job name')
    parser.add_argument('-bayes',
                        '--bayes_file',
                        metavar='',
                        type=str,
                        default="",
                        help='bayes file')
    parser.add_argument('-mem',
                        '--memory_per_job',
                        metavar='',
                        type=int,
                        default="2",
                        help='memory per job (GB)')

    if len(sys.argv) < 2:
        parser.print_help()
        exit(0)

    para_dict = vars(parser.parse_args())

    para_dict["bayesFlag"] = False
    if para_dict["bayes_file"] != "":
        para_dict["bayesFlag"] = True

    main(para_dict)
