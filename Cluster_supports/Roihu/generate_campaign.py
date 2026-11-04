#!/usr/bin/env python3
"""
    Generates a simulation campaign on Roihu (CSC): one working folder per
    design point (parameter file given with -b), and HyperQueue jobs that run
    the event folders of a group of design points on full nodes.

    Each HyperQueue job runs every event folder of its design points without
    results as one task, then collects each finished design point into
    RESULTS/<design point>/<design point>.h5. Submitting the job again
    continues where it stopped. Calling this script again with more parameter
    files adds design points; existing ones are skipped.
"""

import sys
from os import path, listdir, makedirs
from glob import glob
import subprocess
import argparse
from math import ceil

code_package_path = path.abspath(
    path.join(path.dirname(path.abspath(__file__)), "..", ".."))

# full Roihu-CPU node of the medium partition
cores_per_node = 384


def write_hq_job_script(campaign_folder, job_index, design_points, args):
    """This function writes the HyperQueue job for a group of design points"""
    account = ""
    if args.account != "":
        account = "#SBATCH --account={}\n".format(args.account)
    time_request = ""
    if args.task_time != "":
        time_request = "--time-request={} ".format(args.task_time)
    # full nodes: all their memory; fewer cores (test partition): 2 GiB per
    # core, the most that costs no extra billing units
    memory = "--mem=0"
    if args.cores_per_node < cores_per_node:
        memory = "--mem-per-cpu=2G"
    cleanup = ""
    if args.cleanup:
        cleanup = """
        if [ -f RESULTS/$dp/$dp.h5 ]; then
            rm -fr $dp
        fi"""
    script_name = "hq_job_{:03d}.script".format(job_index)
    script = open(path.join(campaign_folder, script_name), "w")
    script.write("""#!/bin/bash
#SBATCH --job-name=iEBE_hq_{0:03d}
{1}#SBATCH --partition={12:s}
#SBATCH --nodes={2:d}
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task={3:d}
#SBATCH {13:s}
#SBATCH --time={4:s}
#SBATCH --output=hq_job_{0:03d}_%j.out
#SBATCH --error=hq_job_{0:03d}_%j.err

# Runs the event folders of the design points below with HyperQueue and
# collects each finished design point into RESULTS/<design point>/. Submit
# this script again to continue after a time-out, but not while a job of it
# is still running.

designPoints="{5:s}"
image={6:s}

cd {7:s}
module load hyperqueue

export HQ_SERVER_DIR="$PWD/hq-server/$SLURM_JOB_ID"
mkdir -p "$HQ_SERVER_DIR"
hq server start &
until hq job list &> /dev/null ; do sleep 1 ; done

srun --overlap --cpu-bind=none --mpi=none hq worker start --manager slurm \\
    --on-server-lost finish-running --cpus="$SLURM_CPUS_PER_TASK" &
hq worker wait "$SLURM_NTASKS"

# event folders without results of the design points not collected yet
taskList="$HQ_SERVER_DIR/event_folders.txt"
for dp in $designPoints; do
    if [ -f RESULTS/$dp/$dp.h5 ]; then
        continue
    fi
    for eventFolder in $dp/event_*; do
        ls $eventFolder/RESULTS_*.h5 &> /dev/null || echo "$PWD/$eventFolder"
    done
done > "$taskList"
echo "`wc -l < $taskList` event folders to run"

if [ -s "$taskList" ]; then
    hq submit --each-line="$taskList" --cpus={8:d} {9}\\
        --stdout=none --stderr=none -- bash -c \\
        'export OMP_NUM_THREADS={10:d} OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
         cd "$HQ_ENTRY" && bash submit_job.script > job.log 2>&1'
    hq job wait all
fi

hq worker stop all
hq server stop

# collect the design points whose event folders all have results
for dp in $designPoints; do
    if [ -f RESULTS/$dp/$dp.h5 ]; then
        continue
    fi
    nFolders=`ls -d $dp/event_* | wc -l`
    nResults=`ls $dp/event_*/RESULTS_*.h5 2> /dev/null | wc -l`
    if [ "$nFolders" -eq "$nResults" ]; then
        ./collect_results.sh "$image" $dp{11}
    else
        echo "$dp: $nResults of $nFolders event folders finished"
    fi
done
""".format(job_index, account, args.nodes, args.cores_per_node, args.time,
           " ".join(design_points), path.abspath(args.singularity),
           path.abspath(campaign_folder), args.cpus_per_task,
           time_request, args.n_threads, cleanup, args.partition, memory))
    script.close()
    return script_name


def main():
    """This is the main function"""
    parser = argparse.ArgumentParser(
        description='Generate an iEBE-MUSIC campaign on Roihu (CSC)',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-w', '--campaign_folder', type=str, required=True,
                        help='campaign folder (one subfolder per design '
                        'point)')
    parser.add_argument('-par', '--par_dict', type=str, required=True,
                        help='user parameter dictionary file')
    parser.add_argument('-b', '--design_points', type=str, nargs='+',
                        required=True,
                        help='parameter files of the design points; the '
                        'file name without extension names the design point')
    parser.add_argument('-singularity', '--singularity', type=str,
                        required=True, help='path of the Apptainer image')
    parser.add_argument('--account', type=str, default='',
                        help='Slurm account (CSC project)')
    parser.add_argument('-nev', '--n_events', type=int, default=1000,
                        help='number of events per design point')
    parser.add_argument('-n_hydro', '--n_events_per_folder', type=int,
                        default=10,
                        help='events per event folder (one HyperQueue task)')
    parser.add_argument('-n_th', '--n_threads', type=int, default=1,
                        help='threads per event')
    parser.add_argument('--cpus_per_task', type=int, default=0,
                        help='cores HyperQueue reserves per task (default: '
                        'n_threads; more leaves memory for each event)')
    parser.add_argument('--designs_per_job', type=int, default=10,
                        help='design points per HyperQueue job')
    parser.add_argument('--partition', type=str, default='medium',
                        help='partition of the HyperQueue jobs (medium; test '
                        'for trying a setup, at most 15 minutes)')
    parser.add_argument('--cores_per_node', type=int, default=cores_per_node,
                        help='cores per node (fewer than {} only on the test '
                        'partition)'.format(cores_per_node))
    parser.add_argument('--nodes', type=int, default=1,
                        help='full nodes per HyperQueue job (1-6 on medium)')
    parser.add_argument('--time', type=str, default='36:00:00',
                        help='time limit of each HyperQueue job (at most '
                        '36:00:00 on medium)')
    parser.add_argument('--task_time', type=str, default='',
                        help='expected time of one task (e.g. 10h); a task '
                        'is only started if the job has this much time left')
    parser.add_argument('--cleanup', action='store_true',
                        help='delete the working folder of a design point '
                        'after it is collected')
    parser.add_argument('-seed', '--random_seed', type=int, default=-1,
                        help='random seed (-1: random)')
    args = parser.parse_args()
    if args.cpus_per_task <= 0:
        args.cpus_per_task = args.n_threads

    campaign_folder = path.abspath(args.campaign_folder)
    makedirs(campaign_folder, exist_ok=True)
    n_folders = ceil(args.n_events/args.n_events_per_folder)

    new_design_points = []
    for design_file in args.design_points:
        design_point = path.splitext(path.basename(design_file))[0]
        if path.exists(path.join(campaign_folder, design_point)):
            print("{} exists, skipped".format(design_point))
            continue
        status = subprocess.call(
            [sys.executable,
             path.join(code_package_path, "generate_singularity_jobs.py"),
             "-c", "roihu", "-w", design_point,
             "-par", path.abspath(args.par_dict),
             "-b", path.abspath(design_file),
             "-singularity", path.abspath(args.singularity),
             "-n", str(n_folders), "-n_hydro", str(args.n_events_per_folder),
             "-n_th", str(args.n_threads), "-seed", str(args.random_seed),
             "--account", args.account],
            cwd=campaign_folder, stdout=subprocess.DEVNULL)
        if status != 0:
            sys.exit("\U0001F6AB  generating {} failed".format(design_point))
        new_design_points.append(design_point)

    first_index = len(glob(path.join(campaign_folder, "hq_job_*.script")))
    script_list = []
    for i in range(0, len(new_design_points), args.designs_per_job):
        script_list.append(write_hq_job_script(
            campaign_folder, first_index + len(script_list),
            new_design_points[i:i + args.designs_per_job], args))

    print("\U0001F375  {} design points with {} event folders of {} events "
          "each".format(len(new_design_points), n_folders,
                        args.n_events_per_folder))
    if script_list:
        print("\U0001F680  Submit in {} with: for job in {}; do sbatch $job; "
              "done".format(campaign_folder, " ".join(script_list)))


if __name__ == "__main__":
    main()
