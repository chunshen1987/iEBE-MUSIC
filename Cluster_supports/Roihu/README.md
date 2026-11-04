# Running iEBE-MUSIC on Roihu (CSC) with Apptainer

The jobs run the pre-built iEBE-MUSIC container with Apptainer. Each job (event folder) is one task of a Slurm job array on the `small` partition. Its events run on the node-local disk (`$TMPDIR`), and its results are copied to the event folder at the end. The nodes need only `apptainer`: the events and the collection of results run inside the container.

Replace `<project>` with your CSC project, e.g. `project_2001234`.

## 1. Get the code and the container

On the Roihu-CPU login node (`roihu-cpu.csc.fi`):

    cd /projappl/<project>/$USER
    git clone https://github.com/chunshen1987/iEBE-MUSIC

    # the Apptainer cache would otherwise fill the 15 GiB home quota
    export APPTAINER_CACHEDIR=/scratch/<project>/$USER/.apptainer
    apptainer build iebe-music_dev.sif docker://chunshen1987/iebe-music:dev

The image is about 1 GB. Build it on Roihu-CPU (x86_64): an image built on Roihu-GPU (aarch64) does not run on the CPU nodes. Any other tag of [chunshen1987/iebe-music](https://hub.docker.com/r/chunshen1987/iebe-music) works the same way.

## 2. Generate the jobs

Generate the jobs in the scratch area of the project:

    cd /scratch/<project>/$USER
    python3 /projappl/<project>/$USER/iEBE-MUSIC/generate_singularity_jobs.py \
        -c roihu --account <project> \
        -singularity /projappl/<project>/$USER/iebe-music_dev.sif \
        -par parameters_dict_user.py -w RUN_NAME \
        -n 100 -n_hydro 5 -n_th 4

- `-n`: number of jobs (array tasks); `-n_hydro`: events per job; `-n_th`: cores per job.
- `--node_type`: the partition, `small` (default; up to 72 h), `longrun` (up to 10 days and 192 cores) or `test` (15 min, for trying a setup).
- The time limit of each job is `walltime` of the `control_dict` in the parameter file.
- `--mem_per_cpu`: memory per core in GiB (default 2). Up to 2 GiB per core costs no extra billing units on Roihu; more memory is billed instead of the cores.
- Without `--account`, give the project at submission: `sbatch --account=<project> ...`.
- If `python3` is not available, run the script in the container: `apptainer exec --bind=/scratch,/projappl <image> python3 .../generate_singularity_jobs.py ...`.

## 3. Submit and check the jobs

    sbatch RUN_NAME/submit_jobArr.script
    squeue --me

The output of each job is in `RUN_NAME/job_<array id>_<task>.out` and `.err`. A finished job leaves `RESULTS_<task>.h5` in its event folder.

To rerun failed jobs, resubmit only their tasks; a job whose folder already has its results does nothing:

    sbatch --array=3,17 RUN_NAME/submit_jobArr.script

## 4. Collect the results

From the folder where the jobs were generated (it holds `collect_events_singularity.sh` and `combine_multiple_hdf5.py`):

    apptainer exec --bind=/scratch,/projappl <image> \
        ./collect_events_singularity.sh RUN_NAME RESULTS

This combines the results of all jobs into `RESULTS/RUN_NAME/RUN_NAME.h5`.

## Notes

- `$TMPDIR` holds 20 GiB per job on `small` and `longrun`. Saving hydro surfaces or UrQMD files (`control_dict`) for many events per job can exceed it.
- If a job reaches its time limit, the events in its `$TMPDIR` are lost; resubmit the task (step 3).
- CSC documentation: [Roihu](https://docs.csc.fi/computing/systems-roihu/), [batch job partitions](https://docs.csc.fi/computing/running/batch-job-partitions/), [Apptainer containers](https://docs.csc.fi/computing/containers/overview/), [billing](https://docs.csc.fi/computing/hpc-billing/).
