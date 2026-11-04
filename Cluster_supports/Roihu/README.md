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

`Cluster_supports/Roihu/generate_campaign.py` makes one working folder per design point (parameter file) in a campaign folder, and HyperQueue jobs that each run the event folders of a group of design points on full `medium` nodes (384 cores each), one HyperQueue task per event folder:

    cd /scratch/<project folder>/$USER
    python3 iEBE-MUSIC/Cluster_supports/Roihu/generate_campaign.py \
        -w CAMPAIGN -par parameters_dict_user.py -b designs/*.dat \
        -singularity <image> --account <account> \
        -nev 1000 -n_hydro 10 -n_th 1 \
        --designs_per_job 10 --nodes 2 --time 36:00:00 --task_time 10h
    cd CAMPAIGN
    for job in hq_job_*.script; do sbatch $job; done

- `-nev`: events per design point, in event folders of `-n_hydro` events (one task each; `-nev` must be a multiple of `-n_hydro`); `-n_th`: threads per event.
- `--cpus_per_task`: cores reserved per task (default `-n_th`). A full node has about 1.9 GiB of memory and 1.5 GiB of `$TMPDIR` per core; reserve more cores for events that need more.
- `--designs_per_job`, `--nodes` (1–6) and `--time` (at most 36 h): the size of each HyperQueue job. `--task_time`: the expected time of one task; HyperQueue starts a task only if the job has that much time left.
- `--partition test --cores_per_node 4 --time 00:15:00`: for trying a setup with a few small events on 4 cores (2 GiB each) instead of a full node.
- The name of a parameter file (without extension) names its design point. Calling the script again with more parameter files adds design points and job scripts; existing design points are skipped. Generating the design points in batches keeps the number of files on `/scratch` within the quota (`csc-workspaces`).

Each job runs the event folders of its design points that have no results yet, then collects every finished design point into `RESULTS/<design point>/<design point>.h5`, together with its parameter files, with `collect_results.sh` (section 4). The output of each task is in `job.log` in its event folder, the output of the job in `hq_job_<n>_<job id>.out`. If a job ends before all its events are done, submit the same script again: finished event folders and collected design points are skipped, and a failed collection is retried without running its events again. Do not run two jobs of the same script at the same time. With `--cleanup`, the working folder of a design point is deleted after its results are collected successfully, which frees its files on `/scratch`.

A design point is collected only when all its event folders have results. To collect one by hand, e.g. when an event folder keeps failing, run in the campaign folder:

    ./collect_results.sh <image> <design point>

A design point collected this way counts as finished: its job no longer runs its remaining event folders.

## Notes

- `$TMPDIR` holds 20 GiB per job on `small` and `longrun`, and 600 GiB per node on `medium`. Saving hydro surfaces or UrQMD files (`control_dict`) for many events can exceed it.
- If a job reaches its time limit, the events running in its `$TMPDIR` are lost; resubmit the task (section 3) or the campaign job (section 5).
- `-n_th` also sets the number of UrQMD runs per event, and with it the number of sampled hadrons: each run samples `afterburner_nIter + 1` times `number_of_particles_needed/10` hadrons (`control_dict` and `iss_dict`). With `-n_th 1`, raise `afterburner_nIter` to keep the statistics of a run with more threads.
- `/scratch` is cleaned of files not used for 90–180 days; copy the results elsewhere in time.
- CSC documentation: [Roihu](https://docs.csc.fi/computing/systems-roihu/), [batch job partitions](https://docs.csc.fi/computing/running/batch-job-partitions/), [Apptainer containers](https://docs.csc.fi/computing/containers/overview/), [HyperQueue](https://docs.csc.fi/apps/hyperqueue/), [billing](https://docs.csc.fi/computing/hpc-billing/).
