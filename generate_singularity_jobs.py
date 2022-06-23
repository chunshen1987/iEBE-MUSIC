#!/usr/bin/env python3
"""
    This script generate all the running jobs to run with the pre-built
    singularity container
"""

import sys
from os import path, mkdir
import shutil
import subprocess
import argparse
from math import ceil
from glob import glob

support_cluster_list = ["wsugrid", "OSG", "local", "stampede2"]


def write_script_header(cluster, script, n_threads, event_id, walltime,
                        working_folder):
    """This function write the header of the job submission script"""
    mem = 4*n_threads
    if cluster == "wsugrid":
        script.write("""#!/usr/bin/env bash
#SBATCH --job-name event_{0}
#SBATCH -q primary
#SBATCH -N 1
#SBATCH -n {1}
#SBATCH --mem={2:.0f}G
#SBATCH --constraint=intel
#SBATCH -t {3:s}
#SBATCH -e job.err
#SBATCH -o job.log

cd {4:s}
""".format(event_id, n_threads, mem, walltime, working_folder))
    elif cluster in ("local", "OSG"):
        script.write("#!/bin/bash")
    elif cluster == "stampede2":
        script.write("""#!/usr/bin/env bash
module load tacc-singularity
""")
    else:
        print("\U0001F6AB  unrecoginzed cluster name :", cluster)
        print("Available options: ", support_cluster_list)
        exit(1)


def generate_Stampede2_mpi_job_script(folder_name, nodeType, n_nodes, n_jobs,
                                      n_threads, walltime):
    """This function generates job script for Stampede2"""
    working_folder = folder_name

    queueName = '{}-normal'.format(nodeType)
    if nodeType == "knl":
        queueName = 'normal'

    script = open(path.join(working_folder, "submit_MPI_jobs.script"), "w")
    script.write("""#!/bin/bash -l
#SBATCH -J iEBEMUSIC
#SBATCH -o job.o%j
#SBATCH -o job.e%j
#SBATCH -p {0:s}
#SBATCH -N {1:d}
#SBATCH -n {2:d}
#SBATCH -t {3:s}
#SBATCH -A TG-PHY200093

module load python3

export OMP_PROC_BIND=true
export OMP_PLACES=threads
export OMP_NUM_THREADS={4:d}

ibrun python3 job_MPI_wrapper.py

""".format(queueName, n_nodes, n_jobs, walltime, n_threads))
    script.close()


def generate_event_folders(workingFolder, clusterName, eventId,
                           singularityRepoPath, executeScript, parameterFile,
                           eventId0, nHydroEvents, nThreads, seed, wallTime):
    """This function creates the event folder structure"""
    eventFolder = path.join(workingFolder, 'event_{}'.format(eventId))
    mkdir(eventFolder)

    # generate job running script
    executeScriptName = executeScript.split('/')[-1]
    parameterFileName = parameterFile.split('/')[-1]
    script = open(path.join(eventFolder, "submit_job.script"), "w")
    write_script_header(clusterName, script, nThreads, eventId, wallTime,
                        eventFolder)
    script.write("""
singularity exec {0} ./{1} {2} {3} {4} {5} {6}

mkdir -p temp
./collect_events.sh playground temp
mv temp/playground/playground.h5 RESULTS_{7}.h5
""".format(singularityRepoPath, executeScriptName, parameterFileName,
           eventId0, nHydroEvents, nThreads, seed, eventId))
    script.close()

    # copy files
    shutil.copy(executeScript, eventFolder)
    shutil.copy(parameterFile, eventFolder)


def create_a_working_folder(workfolder_path):
    try:
        mkdir(workfolder_path)
    except FileExistsError:
        print("The folder {} exists, do you want to delete it?".format(
            workfolder_path))
        user_answer = input()
        if 'y' in user_answer:
            shutil.rmtree(workfolder_path)
            mkdir(workfolder_path)
        else:
            print("bye~\n")
            exit(0)


def main():
    """This is the main funciton"""
    parser = argparse.ArgumentParser(
        description='\U0000269B Welcome to iEBE-MUSIC package',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-w',
                        '--working_folder_name',
                        metavar='',
                        type=str,
                        default='playground',
                        help='working folder path')
    parser.add_argument('-c',
                        '--cluster_name',
                        metavar='',
                        type=str,
                        default='local',
                        help='name of the cluster')
    parser.add_argument('--node_type',
                        metavar='',
                        type=str,
                        default='SKX',
                        help='name of the queue (work on stampede2)')
    parser.add_argument('-n',
                        '--n_jobs',
                        metavar='',
                        type=int,
                        default=1,
                        help='number of jobs')
    parser.add_argument('-n_hydro',
                        '--n_hydro_per_job',
                        metavar='',
                        type=int,
                        default=1,
                        help='number of hydro events per job to run')
    parser.add_argument('-n_th',
                        '--n_threads',
                        metavar='',
                        type=int,
                        default=1,
                        help='number of threads used for each job')
    parser.add_argument('-par',
                        '--par_dict',
                        metavar='',
                        type=str,
                        default='parameters_dict_user.py',
                        help='user-defined parameter dictionary file')
    parser.add_argument('-singularity',
                        '--singularity',
                        metavar='',
                        type=str,
                        default='iebe-music_latest.sif',
                        help='path of the singularity image')
    parser.add_argument('-exe',
                        '--executeScript',
                        metavar='',
                        type=str,
                        default='Cluster_supports/WSUgrid/run_singularity.sh',
                        help='job running script')
    parser.add_argument('-b',
                        '--bayes_file',
                        metavar='',
                        type=str,
                        default='',
                        help='parameters from bayesian analysis')
    parser.add_argument('-seed',
                        '--random_seed',
                        metavar='',
                        type=int,
                        default='-1',
                        help='Random Seed (-1: according to system time)')
    args = parser.parse_args()

    if len(sys.argv) < 2:
        parser.print_help()
        exit(0)

    # print out all the arguments
    print("="*40)
    print("\U0000269B   Input parameters")
    print("="*40)
    for iarg in vars(args):
        print("\U0000269B   {}: {}".format(iarg, getattr(args, iarg)))
    print("="*40)

    try:
        working_folder_name = args.working_folder_name
        cluster_name = args.cluster_name.lower()
        n_jobs = args.n_jobs
        n_hydro_per_job = args.n_hydro_per_job
        n_threads = args.n_threads
        seed = args.random_seed
        singularityRepoPath = path.abspath(args.singularity)
        executeScript = args.executeScript
        parameterFile = args.par_dict
    except:
        parser.print_help()
        exit(0)

    code_package_path = path.abspath(path.dirname(__file__))
    par_diretory = path.dirname(path.abspath(args.par_dict))
    sys.path.insert(0, par_diretory)
    parameter_dict = __import__(args.par_dict.split('.py')[0].split("/")[-1])
    wallTime = parameter_dict.control_dict['walltime']

    working_folder_name = path.abspath(working_folder_name)
    create_a_working_folder(working_folder_name)

    shutil.copy(args.par_dict, working_folder_name)

    toolbar_width = 40
    sys.stdout.write("\U0001F375  Generating {} jobs [{}]".format(
        n_jobs, " "*toolbar_width))
    sys.stdout.flush()
    sys.stdout.write("\b"*(toolbar_width + 1))
    for i_job in range(n_jobs):
        progress_i = (int(float(i_job + 1)/n_jobs*toolbar_width)
                      - int(float(i_job)/n_jobs*toolbar_width))
        for ii in range(progress_i):
            sys.stdout.write("#")
            sys.stdout.flush()
        generate_event_folders(working_folder_name, cluster_name, i_job,
                               singularityRepoPath, executeScript,
                               parameterFile, i_job*n_hydro_per_job,
                               n_hydro_per_job, n_threads, seed, wallTime)
    sys.stdout.write("\n")
    sys.stdout.flush()

    # copy script to collect final results
    pwd = path.abspath(".")
    script_path = path.join(code_package_path, "utilities")
    shutil.copy(path.join(script_path, 'collect_events_singularity.sh'), pwd)
    shutil.copy(path.join(script_path, 'combine_multiple_hdf5.py'), pwd)

    if cluster_name == "wsugrid":
        shutil.copy(
            path.join(code_package_path,
                      'Cluster_supports/WSUgrid/submit_all_jobs.sh'), pwd)

    if cluster_name == "stampede2":
        nThreadsPerNode = 1
        if args.node_type.lower() == "skx":
            nThreadsPerNode = 96
        elif args.node_type.lower() == "knl":
            nThreadsPerNode = 272
        elif args.node_type.lower() == "icx":
            nThreadsPerNode = 160
        shutil.copy(
            path.join(code_package_path,
                      'Cluster_supports/Stampede2/job_MPI_wrapper.py'),
            working_folder_name)
        n_nodes = max(1, int(n_jobs*n_threads/nThreadsPerNode))
        if n_nodes*nThreadsPerNode < n_jobs*n_threads:
            n_nodes += 1

        generate_Stampede2_mpi_job_script(working_folder_name,
                                          args.node_type.lower(),
                                          n_nodes, n_jobs, n_threads, wallTime)


if __name__ == "__main__":
    main()
