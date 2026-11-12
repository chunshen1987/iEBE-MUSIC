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
import re
import shlex
from os import path, makedirs
from glob import glob
from math import ceil
import subprocess
import argparse

code_package_path = path.abspath(
    path.join(path.dirname(path.abspath(__file__)), "..", ".."))

# full Roihu-CPU node of the medium partition
cores_per_node = 384
# memory per core of a full node, in GiB (744 GiB for 384 cores)
memory_per_core = 744/384


def is_hq_duration(duration):
    """Whether HyperQueue reads the duration, e.g. 3h, 90min or 1h 30m"""
    return re.fullmatch(
        r"(\s*\d+\s*(days?|d|hours?|hrs?|h|minutes?|mins?|m|seconds?|secs?|s))+"
        r"\s*", duration) is not None


def print_resources(args, n_folders, n_design_points):
    """Prints how the tasks of a job fill its cores"""
    slots = args.nodes*(args.cores_per_node//args.cpus_per_task)
    designs = min(args.designs_per_job, n_design_points)
    tasks = designs*n_folders
    rounds = ceil(tasks/slots)
    print("\U0001F4CB  Each job: {} design points x {} event folders = {} "
          "tasks on {} slots ({} nodes x {} cores, {} cores per task, "
          "{:.1f} GiB memory per task)".format(
              designs, n_folders, tasks, slots, args.nodes,
              args.cores_per_node, args.cpus_per_task,
              args.cpus_per_task*memory_per_core))
    print("\U0001F4CB  The tasks run in {} rounds, so a job takes about {} "
          "times the median task time (--time {}), and {:.0%} of the "
          "reserved cores are used".format(rounds, rounds, args.time,
                                           tasks/(rounds*slots)))


def write_hq_job_script(campaign_folder, job_index, design_points, args):
    """This function writes the HyperQueue job for a group of design points"""
    account = ""
    if args.account != "":
        account = "#SBATCH --account={}\n".format(args.account)
    time_request = ""
    if args.task_time != "":
        # quoted, so that e.g. 1h 30m stays one argument
        time_request = "--time-request={} ".format(
            shlex.quote(args.task_time))
    # full nodes: all their memory; fewer cores (test partition): 2 GiB per
    # core, the most that costs no extra billing units
    memory = "--mem=0"
    if args.cores_per_node < cores_per_node:
        memory = "--mem-per-cpu=2G"
    cleanup = ""
    if args.cleanup:
        # only after a successful collection
        cleanup = " && rm -fr $dp"
    slots = args.nodes*(args.cores_per_node//args.cpus_per_task)
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
#
# Resources: {2:d} nodes of {3:d} cores. Each event folder is one task on
# {8:d} cores, running its events one after another with {10:d} threads, so
# {14:d} tasks run at the same time; a task that ends frees its cores for
# the next one.{15}

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

# event folders without results of the design points not collected yet; a
# design point with results kept from a failed collection was complete
taskList="$HQ_SERVER_DIR/event_folders.txt"
for dp in $designPoints; do
    if [ -f RESULTS/$dp/$dp.h5 ] || [ -d RESULTS/$dp/SPVN_RESULTS ]; then
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
         cd "$HQ_ENTRY" || exit 1
         start=`date +%s`
         bash submit_job.script > job.log 2>&1
         status=$?
         echo $((`date +%s` - start)) > task_seconds.txt
         exit $status'
    hq job wait all
fi

hq worker stop all
hq server stop

# how long the finished tasks took, for choosing the settings of the next
# jobs (--task_time, --designs_per_job)
for dp in $designPoints; do
    for eventFolder in $dp/event_*; do
        if [ -f $eventFolder/task_seconds.txt ] \\
                && ls $eventFolder/RESULTS_*.h5 &> /dev/null; then
            cat $eventFolder/task_seconds.txt
        fi
    done
done | sort -n | LC_ALL=C awk '{{t[NR] = $1}}
    END {{if (NR > 0) {{
        median = (NR % 2) ? t[(NR + 1)/2] : (t[NR/2] + t[NR/2 + 1])/2
        printf "%d finished tasks took %.2f h (median), %.2f h at most\\n", NR, median/3600, t[NR]/3600}}}}'

# collect the design points whose event folders all have results, and retry
# failed collections; for the others, show why their tasks have no results
for dp in $designPoints; do
    if [ -f RESULTS/$dp/$dp.h5 ]; then
        continue
    fi
    nFolders=`ls -d $dp/event_* | wc -l`
    nResults=`ls $dp/event_*/RESULTS_*.h5 2> /dev/null | wc -l`
    if [ "$nFolders" -eq "$nResults" ] || [ -d RESULTS/$dp/SPVN_RESULTS ]; then
        ./collect_results.sh "$image" $dp{11}
    else
        nStarted=0
        firstFolder=""
        for eventFolder in $dp/event_*; do
            if [ -f $eventFolder/job.log ] \\
                    && ! ls $eventFolder/RESULTS_*.h5 &> /dev/null; then
                nStarted=$((nStarted + 1))
                firstFolder=${{firstFolder:-$eventFolder}}
            fi
        done
        echo "$dp: $nResults of $nFolders event folders finished;" \\
             "$nStarted ran without results (failed, or stopped by the" \\
             "time limit)"
        if [ -n "$firstFolder" ]; then
            echo "  the end of $firstFolder/job.log:"
            tail -n 20 $firstFolder/job.log | sed 's/^/    /'
        fi
    fi
done
""".format(job_index, account, args.nodes, args.cores_per_node, args.time,
           " ".join(design_points), path.abspath(args.singularity),
           path.abspath(campaign_folder), args.cpus_per_task,
           time_request, args.n_threads, cleanup, args.partition, memory,
           slots,
           "" if args.task_time == "" else (
               "\n# A task is only started if the job has {} left."
               .format(args.task_time))))
    script.close()
    return script_name


def main():
    """This is the main function"""
    parser = argparse.ArgumentParser(
        description='Generate an iEBE-MUSIC campaign on Roihu (CSC): one '
        'working folder per design point, and HyperQueue jobs on full nodes '
        'that run its event folders as tasks. See '
        'Cluster_supports/Roihu/README.md for examples.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    what = parser.add_argument_group('what to simulate')
    what.add_argument('-w', '--campaign_folder', type=str, required=True,
                      help='campaign folder: holds a working folder per '
                      'design point, the job scripts and RESULTS/')
    what.add_argument('-par', '--par_dict', type=str, required=True,
                      help='user parameter dictionary, the same for all '
                      'design points')
    what.add_argument('-b', '--design_points', type=str, nargs='+',
                      required=True,
                      help='parameter files of the design points, one per '
                      'design point; the file name without extension names '
                      'the design point, its folder and its results')
    what.add_argument('-singularity', '--singularity', type=str,
                      required=True,
                      help='path of the Apptainer image, used for the events '
                      'and the collection')
    what.add_argument('-nev', '--n_events', type=int, default=1000,
                      help='events per design point')
    what.add_argument('-seed', '--random_seed', type=int, default=-1,
                      help='random seed (-1: random)')
    tasks = parser.add_argument_group(
        'tasks: an event folder is one HyperQueue task, which runs its '
        'events one after another on the same cores')
    tasks.add_argument('-n_hydro', '--n_events_per_folder', type=int,
                       default=10,
                       help='events per event folder (task); -nev must be a '
                       'multiple. Fewer events per folder give more, '
                       'shorter tasks, which fill the nodes more evenly, but '
                       'more files on /scratch')
    tasks.add_argument('-n_th', '--n_threads', type=int, default=1,
                       help='threads per event (OpenMP threads of IP-Glasma '
                       'and MUSIC); also the number of UrQMD runs per event')
    tasks.add_argument('--cpus_per_task', type=int, default=0,
                       help='cores reserved per task (default: -n_th). A core '
                       'of a full node comes with {:.1f} GiB of memory; '
                       'reserve more cores for events that need more '
                       'memory'.format(memory_per_core))
    tasks.add_argument('--task_time', type=str, default='',
                       help='the longest time a task takes (e.g. 3h): '
                       'HyperQueue starts a task only if the job has this '
                       'much time left, so that no task is cut off at the end '
                       'of the job')
    jobs = parser.add_argument_group(
        'jobs: a job script runs the tasks of a group of design points on '
        'its nodes')
    jobs.add_argument('--designs_per_job', type=int, default=10,
                      help='design points per job script; their tasks are '
                      'run as one queue, so that the cores stay busy between '
                      'design points')
    jobs.add_argument('--nodes', type=int, default=1,
                      help='full nodes per job (1-6 on medium); a job runs '
                      'nodes x cores_per_node / cpus_per_task tasks at the '
                      'same time')
    jobs.add_argument('--time', type=str, default='36:00:00',
                      help='time limit of each job (at most 36:00:00 on '
                      'medium); a job that reaches it is submitted again to '
                      'continue')
    jobs.add_argument('--account', type=str, default='',
                      help='Slurm account (CSC project)')
    jobs.add_argument('--partition', type=str, default='medium',
                      help='partition of the jobs: medium (full nodes), or '
                      'test for trying a setup (at most 15 minutes)')
    jobs.add_argument('--cores_per_node', type=int, default=cores_per_node,
                      help='cores per node; fewer than {} only on the test '
                      'partition'.format(cores_per_node))
    jobs.add_argument('--cleanup', action='store_true',
                      help='delete the working folder of a design point after '
                      'its results are collected, to free its files on '
                      '/scratch')
    args = parser.parse_args()
    if args.n_events % args.n_events_per_folder != 0:
        parser.error("-nev ({}) must be a multiple of -n_hydro ({})".format(
            args.n_events, args.n_events_per_folder))
    if args.cpus_per_task <= 0:
        args.cpus_per_task = args.n_threads
    if args.task_time != "" and not is_hq_duration(args.task_time):
        parser.error("--task_time {} is not a duration such as 3h or "
                     "90min".format(args.task_time))

    campaign_folder = path.abspath(args.campaign_folder)
    makedirs(campaign_folder, exist_ok=True)
    n_folders = args.n_events//args.n_events_per_folder

    new_design_points = []
    failed_design_point = ""
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
            # the design points generated so far still get their job scripts
            failed_design_point = design_point
            break
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
    if new_design_points:
        print_resources(args, n_folders, len(new_design_points))
    if script_list:
        print("\U0001F680  Submit in {} with: for job in {}; do sbatch $job; "
              "done".format(campaign_folder, " ".join(script_list)))
    if failed_design_point != "":
        sys.exit("\U0001F6AB  Generating {} failed (see above), so it and "
                 "the design points after it were not generated. Fix it and "
                 "run this command again: the design points generated so far "
                 "are skipped.".format(failed_design_point))


if __name__ == "__main__":
    main()
