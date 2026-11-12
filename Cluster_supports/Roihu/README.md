# Running iEBE-MUSIC on Roihu (CSC) with Apptainer

The jobs run the pre-built iEBE-MUSIC container with Apptainer, so the nodes need only `apptainer`: the events and the collection of results run inside the container. Events run on the node-local disk (`$TMPDIR`), and their results are copied to the event folder at the end.

There are two ways to run:
- **A job array** (sections 2–4) for single runs: each job (event folder) is one task of a Slurm job array on the `small` partition.
- **A campaign** (section 5) for many design points with thousands of events each: HyperQueue runs the event folders of several design points on full nodes.

`<account>` is the Slurm account of your CSC project, as listed by `sacctmgr show assoc user=$USER format=account`; it can be the project name (e.g. `lappi`) or its number (e.g. `project_2001234`).

## 1. Get the code and the container

On the Roihu-CPU login node (`roihu-cpu.csc.fi`):

    git clone https://github.com/chunshen1987/iEBE-MUSIC

    # the Apptainer cache would otherwise fill the 15 GiB home quota
    export APPTAINER_CACHEDIR=/scratch/<project folder>/$USER/.apptainer
    apptainer pull iebe-music_dev.sif docker://chunshen1987/iebe-music:dev

`apptainer build iebe-music_dev.sif docker://...` gives the same image. The image is about 1 GB and can be in your home folder or in `/projappl`; `/scratch` is cleaned of files not used for 90–180 days. Pull it on Roihu-CPU (x86_64): an image built on Roihu-GPU (aarch64) does not run on the CPU nodes. Any other tag of [chunshen1987/iebe-music](https://hub.docker.com/r/chunshen1987/iebe-music) works the same way. Use the same image file for all runs that are compared with each other.

## 2. Generate the jobs

Generate the jobs in the scratch area of the project:

    cd /scratch/<project folder>/$USER
    python3 iEBE-MUSIC/generate_singularity_jobs.py \
        -c roihu --account <account> \
        -singularity <image> \
        -par parameters_dict_user.py -w RUN_NAME \
        -n 100 -n_hydro 5 -n_th 4

- `-n`: number of jobs (array tasks); `-n_hydro`: events per job; `-n_th`: cores per job.
- `-b`: the parameter file of a design point, as on the other clusters.
- **Parameter check:** before writing any job, the generator applies the user dictionary and the parameter file as each job will. A key of the parameter file that no dictionary has, e.g. a misspelled one, stops it. The user dictionary only needs the parameters it changes; the master (`config/parameters_dict_master.py`) has the defaults.
- `--node_type`: the partition, `small` (default; up to 72 h), `longrun` (up to 10 days and 192 cores) or `test` (15 min, for trying a setup).
- The time limit of each job is `walltime` of the `control_dict` in the parameter file.
- `--mem_per_cpu`: memory per core in GiB (default 2). Up to 2 GiB per core costs no extra billing units on Roihu; more memory is billed instead of the cores.
- Without `--account`, give the account at submission: `sbatch --account=<account> ...`.
- Each array task counts as one job against the limits of the project (`sacctmgr show assoc user=$USER format=account,partition,maxjobs,maxsubmit`), e.g. 100 running and 200 submitted jobs on `small`.

## 3. Submit and check the jobs

    sbatch RUN_NAME/submit_jobArr.script
    squeue --me

The output of each job is in `RUN_NAME/job_<array id>_<task>.out` and `.err`. A finished job leaves `RESULTS_<task>.h5` in its event folder.

To rerun failed jobs, resubmit only their tasks; a job whose folder already has its results does nothing:

    sbatch --array=3,17 RUN_NAME/submit_jobArr.script

## 4. Collect the results

The login nodes have no Python with h5py, so the results are collected in the container by `collect_results.sh`. Run it in the folder where the jobs were generated; `generate_singularity_jobs.py` copied it there:

    cd /scratch/<project folder>/$USER
    ./collect_results.sh <image> RUN_NAME

This adds the results of the finished jobs to `RESULTS/RUN_NAME/RUN_NAME.h5`, copies the parameter files next to it, and prints the number of events in the file.

Collect when all jobs are done, and rerun failed jobs (section 3) before. The collection moves the results of each job from its event folder into this file: a job resubmitted after its results were collected runs again, and collecting again adds only the jobs that have finished since. A failed collection (e.g. a damaged file or a full disk) deletes nothing: the results stay in `RESULTS/RUN_NAME/SPVN_RESULTS`, and collecting again retries.

## 5. Campaigns: many design points with HyperQueue

`Cluster_supports/Roihu/generate_campaign.py` runs many events of one or many design points on full `medium` nodes. It uses [HyperQueue](https://docs.csc.fi/apps/hyperqueue/), which runs many small tasks inside one Slurm job.

### 5.1 How a campaign runs

- **Design point:** one parameter file (`-b`), applied on top of the user dictionary (`-par`). Its file name without extension names its working folder in the campaign folder and its results.
- **Event folder = task:** each design point has `-nev` / `-n_hydro` event folders of `-n_hydro` events. Each event folder is one HyperQueue task. A task runs its events one after another on `--cpus_per_task` cores, with `-n_th` threads per event.
- **Slots:** a node has 384 cores, so it runs 384 / `--cpus_per_task` tasks at the same time, e.g. 192 tasks of 2 cores. When a task ends, its cores take the next task.
- **Job:** a job script `hq_job_<n>.script` reserves `--nodes` full nodes for up to `--time`, and runs the tasks of `--designs_per_job` design points as one queue. The tasks run in "rounds" (tasks / slots). The job takes about rounds × the time of a task.
- **Cost:** a job is billed for its full nodes while it runs, whether its cores are busy or not. The settings below keep them busy.

When `generate_campaign.py` writes the jobs, it prints how the tasks fill a job: slots, rounds and the share of cores used.

### 5.2 The options

| Option | What it controls | How to choose it |
|---|---|---|
| `-w` | the campaign folder: working folders of the design points, job scripts, `RESULTS/` | one per campaign |
| `-par` | the user dictionary, the same for all design points | it only needs the parameters it changes |
| `-b` | the parameter files, one per design point | a key that no dictionary has stops the generation |
| `-singularity` | the Apptainer image of the events and of the collection | the same image for all runs that are compared |
| `-nev` | events per design point | a multiple of `-n_hydro` |
| `-n_hydro` | events per event folder, i.e. per task | **1**: the shortest tasks fill the slots most evenly, and a failure costs one event. 2 halves the number of folders and files on `/scratch`. |
| `-n_th` | threads per event (IP-Glasma, MUSIC), and UrQMD runs per event | the transport runs on one core, so 1 or 2 |
| `--cpus_per_task` | cores reserved per task (default `-n_th`) | each core of a full node comes with 1.9 GiB of memory and 1.5 GiB of `$TMPDIR`; reserve more cores for events that need more |
| `--task_time` | HyperQueue starts a task only if the job has this much time left, so no task is cut off at the end of a job | about the longest task time, e.g. twice the median |
| `--designs_per_job` | design points per job; their tasks form one queue, so the cores stay busy between design points | as many as fit into `--time` (5.4) |
| `--nodes` | full nodes per job (1–6 on `medium`) | **1**: more nodes finish sooner at the same cost, but wait longer in the queue |
| `--time` | time limit of a job (at most 36 h on `medium`) | about rounds × median task time, plus a margin. A job that reaches it is simply submitted again. |
| `--cleanup` | delete a design point's working folder after its results are collected | for many design points, to stay within the file quota |
| `--account` | the Slurm account (CSC project) | |
| `--partition test --cores_per_node 4 --time 00:15:00` | a test on 4 cores (2 GiB each) instead of a full node | for trying a setup with a few small events |

### 5.3 Example: one design point with 1000 events

With 2 threads per event and events of about 1.5 h (measure them, see 5.5):

    cd /scratch/<project folder>/$USER
    python3 iEBE-MUSIC/Cluster_supports/Roihu/generate_campaign.py \
        -w RUN_NAME -par parameters_dict_user.py -b designs/0000.dat \
        -singularity <image> --account <account> \
        -nev 1000 -n_hydro 1 -n_th 2 \
        --nodes 1 --time 12:00:00 --task_time 3h
    cd RUN_NAME
    sbatch hq_job_000.script

- **Layout:** 1000 tasks of one event on 192 slots run in 6 rounds, about 9 h with 1.5 h events.
- **Use of the cores:** 87%, since the last round fills only 40 of the 192 slots. This holds for any number of nodes: 2 nodes run 3 rounds in half the time at the same cost. If the number of events is free, a multiple of 192 per node (e.g. 1152) fills every round.
- **`--time 12:00:00`:** 6 rounds of 1.5 h plus a margin. A shorter limit usually starts sooner.

### 5.4 Example: many design points

Put several design points into each 1-node job, so that the cores stay busy between them:

    cd /scratch/<project folder>/$USER
    python3 iEBE-MUSIC/Cluster_supports/Roihu/generate_campaign.py \
        -w CAMPAIGN -par parameters_dict_user.py -b designs/*.dat \
        -singularity <image> --account <account> \
        -nev 1000 -n_hydro 1 -n_th 2 \
        --designs_per_job 3 --nodes 1 --time 36:00:00 --task_time 3h --cleanup
    cd CAMPAIGN
    for job in hq_job_*.script; do sbatch $job; done

- **Layout:** with 1.5 h events, a job runs 3000 tasks in 16 rounds, about 24 h, with 98% of its cores in use.
- **Jobs:** the jobs of the design points run side by side, up to the project's limit of running jobs (100 on `medium`; `sacctmgr show assoc user=$USER format=account,partition,maxjobs,maxsubmit`).
- **`--designs_per_job` by event time:** from the median task time, for 1000 one-event tasks per design point on one node, 2 cores per task, and 36 h with a 15% margin:

  | Median event time | `--designs_per_job` | Rounds | Job length | Cores used |
  |---|---|---|---|---|
  | 0.5 h | 11 | 58 | 29 h | 99% |
  | 1 h | 5 | 27 | 27 h | 96% |
  | 1.5 h | 3 | 16 | 24 h | 98% |
  | 2 h | 2 | 11 | 22 h | 95% |
  | 3–5 h | 1 | 6 | 18–30 h | 87% |
  | more than 5 h | 1, with `--nodes 2` | 3 | 3 × event time | 87% |

  In general: `--designs_per_job` = 0.85 × 36 h × 192 / (1000 × median event time), rounded down.
- **Files on `/scratch`:** with `-n_hydro 1`, each design point has 1000 event folders. Generate the design points in batches, e.g. 50 at a time with `--cleanup`, and check the file count with `csc-workspaces`.

### 5.5 Measure, then adjust

- **Task times:** each job prints the median and the longest time of its finished tasks in `hq_job_<n>_<job id>.out`. Use the median as the event time for `--designs_per_job` (5.4) and `--time`, and about the longest for `--task_time`.
- **Memory:** `seff <job id>` shows the memory a finished job used of the node's 744 GiB. An event must stay within `--cpus_per_task` × 1.9 GiB: 3.9 GiB with 2 cores. IP-Glasma needs 2.4 GiB at 800² and 3.9 GiB at 1024² (`validations/` of IP-Glasma), so use `--cpus_per_task 3` from about 1000².
- **Disk:** a job has 600 GiB of `$TMPDIR` per node, i.e. 3.1 GiB per task with 192 tasks. Saving hydro surfaces or UrQMD files (`control_dict`) can exceed this.

The settings apply to the jobs written next. To change them for design points that are already generated, delete their folders and job scripts, and generate them again.

### 5.6 Restarts, failures and collection

Each job runs the event folders of its design points that have no results yet, then collects every finished design point into `RESULTS/<design point>/<design point>.h5`, together with its parameter files, with `collect_results.sh` (section 4). The output of each task is in `job.log` in its event folder, the output of the job in `hq_job_<n>_<job id>.out`. If a job ends before all its events are done, submit the same script again: finished event folders and collected design points are skipped, and a failed collection is retried without running its events again. Do not run two jobs of the same script at the same time. With `--cleanup`, the working folder of a design point is deleted after its results are collected successfully, which frees its files on `/scratch`.

If tasks end without results, the job's `.out` names the design point, says how many of its tasks ran without results, and shows the end of the `job.log` of one of them.

A design point is collected only when all its event folders have results. To collect one by hand, e.g. when an event folder keeps failing, run in the campaign folder:

    ./collect_results.sh <image> <design point>

A design point collected this way counts as finished: its job no longer runs its remaining event folders.

## Notes

- `$TMPDIR` holds 20 GiB per job on `small` and `longrun`, and 600 GiB per node on `medium`. Saving hydro surfaces or UrQMD files (`control_dict`) for many events can exceed it.
- If a job reaches its time limit, the events running in its `$TMPDIR` are lost; resubmit the task (section 3) or the campaign job (section 5).
- `-n_th` also sets the number of UrQMD runs per event, and with it the number of sampled hadrons: each run samples `afterburner_nIter + 1` times `number_of_particles_needed/10` hadrons (`control_dict` and `iss_dict`). With `-n_th 1`, raise `afterburner_nIter` to keep the statistics of a run with more threads.
- `/scratch` is cleaned of files not used for 90–180 days; copy the results elsewhere in time.
- CSC documentation: [Roihu](https://docs.csc.fi/computing/systems-roihu/), [batch job partitions](https://docs.csc.fi/computing/running/batch-job-partitions/), [Apptainer containers](https://docs.csc.fi/computing/containers/overview/), [HyperQueue](https://docs.csc.fi/apps/hyperqueue/), [billing](https://docs.csc.fi/computing/hpc-billing/).
