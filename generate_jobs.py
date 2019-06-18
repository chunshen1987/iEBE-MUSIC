#!/usr/bin/env python
"""This script generate all the running jobs."""

import sys
from os import path, mkdir
import shutil
import subprocess


def write_script_header(cluster, script, n_threads,
                        event_id, walltime, working_folder):
    """This function write the header of the job submission script"""
    mem = 4*n_threads
    if cluster == "nersc":
        script.write(
            """#!/bin/bash -l
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -J {0:s}
#SBATCH -t {1:s}
#SBATCH -L SCRATCH
#SBATCH -C haswell
""".format(event_id, walltime))
    elif cluster == "nerscKNL":
        script.write(
            """#!/bin/bash -l
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -J {0:s}
#SBATCH -t {1:s}
#SBATCH -L SCRATCH
#SBATCH -C knl,quad,cache
""".format(event_id, walltime))
    elif cluster == "guillimin":
        script.write(
            """#!/usr/bin/env bash
#PBS -N {0:s}
#PBS -l nodes=1:ppn={1:d}
#PBS -l walltime={2:s}
#PBS -S /bin/bash
#PBS -e test.err
#PBS -o test.log
#PBS -A cqn-654-ad
#PBS -q sw
#PBS -d {3:s}
""".format(event_id, n_threads, walltime, working_folder))
    elif cluster == "McGill":
        script.write(
            """#!/usr/bin/env bash
#PBS -N {0:s}
#PBS -l nodes=1:ppn={1:d}:irulan
#PBS -l walltime={2:s}
#PBS -S /bin/bash
#PBS -e test.err
#PBS -o test.log
#PBS -d {3:s}
""".format(event_id, n_threads, walltime, working_folder))
    elif cluster == "wsugrid":
        script.write(
            """#!/usr/bin/env bash
#PBS -N {0:s}
#PBS -l select=1:ncpus={1:d}:mem={2:.0f}GB:cpu_type=Intel
#PBS -l walltime={3:s}
#PBS -S /bin/bash
#PBS -e test.err
#PBS -o test.log
#PBS -q wsuq

cd {4:s}
""".format(event_id, n_threads, mem, walltime, working_folder))
    elif cluster in "local":
        script.write("#!/usr/bin/env bash")
    else:
        print("\U0001F6AB  unrecoginzed cluster name :", cluster)
        print("Available options: nersc, nerscKNL, wsugrid, local, guillimin, "
              + "McGill")
        exit(1)


def generate_nersc_mpi_job_script(folder_name, n_nodes, n_threads,
                                  n_jobs_per_node):
    """This function generates job script for NERSC"""
    working_folder = folder_name
    walltime = '10:00:00'

    script = open(path.join(working_folder, "submit_MPI_jobs.pbs"), "w")
    script.write(
        """#!/bin/bash -l
#SBATCH --qos=regular
#SBATCH -N {0:d}
#SBATCH -A m1820
#SBATCH -J music
#SBATCH -t {1:s}
#SBATCH -L SCRATCH
#SBATCH -C haswell
#SBATCH --mail-type=ALL
#SBATCH --mail-user=chunshen1987@gmail.com

export OMP_PROC_BIND=true
export OMP_PLACES=threads

num_of_nodes={0:d}
# run all the job
for (( nodeid=1; nodeid <= $num_of_nodes; nodeid++ ))
do
    export OMP_NUM_THREADS={2:d}
    srun -N 1 -n {3:d} -c {2:d} python job_MPI_wrapper.py {3:d} $nodeid &
done
wait
""".format(n_nodes, walltime, n_threads, n_jobs_per_node))
    script.close()


def generate_nerscKNL_mpi_job_script(folder_name, n_nodes, n_threads,
                                     n_jobs_per_node):
    """This function generates job script for NERSC KNL"""
    working_folder = folder_name
    walltime = '10:00:00'

    script = open(path.join(working_folder, "submit_MPI_jobs.pbs"), "w")
    script.write(
        """#!/bin/bash -l
#SBATCH --qos=regular
#SBATCH -N {0:d}
#SBATCH -A m1820
#SBATCH -J music
#SBATCH -t {1:s}
#SBATCH -L SCRATCH
#SBATCH -C knl,quad,cache
#SBATCH --mail-type=ALL
#SBATCH --mail-user=chunshen1987@gmail.com

export OMP_PROC_BIND=true
export OMP_PLACES=cores

num_of_nodes={0:d}
# run all the job
for (( nodeid=1; nodeid <= $num_of_nodes; nodeid++ ))
do
    export OMP_NUM_THREADS={2:d}
    srun -N 1 -n {3:d} -c {2:d} python job_MPI_wrapper.py {3:d} $nodeid &
done
wait
""".format(n_nodes, walltime, n_threads, n_jobs_per_node))
    script.close()


def generate_full_job_script(cluster_name, folder_name, database, initial_type,
                             n_hydro, ev0_id, n_urqmd, n_threads):
    """This function generates full job script"""
    working_folder = folder_name
    event_id = working_folder.split('/')[-1]
    walltime = '100:00:00'

    script = open(path.join(working_folder, "submit_job.pbs"), "w")
    write_script_header(cluster_name, script, n_threads, event_id, walltime,
                        working_folder)
    script.write(
        """
./hydro_plus_UrQMD_driver.py {0:s} {1:s} {2:d} {3:d} {4:d} {5:d} > run.log
""".format(initial_type, database, n_hydro, ev0_id, n_urqmd, n_threads))
    script.close()


def generate_script_hydro(folder_name, nthreads):
    """This function generates script for hydro simulation"""
    working_folder = folder_name

    script = open(path.join(working_folder, "run_hydro.sh"), "w")

    hydro_results_folder = 'hydro_results'
    script.write(
        """#!/usr/bin/env bash

results_folder={0:s}

(
cd MUSIC

""".format(hydro_results_folder))

    if nthreads > 0:
        script.write(
            """
export OMP_NUM_THREADS={0:d}
""".format(nthreads))

    script.write(
        """
# hydro evolution
./MUSIChydro music_input_mode_2 1> run.log 2> run.err
./sweeper.sh $results_folder
)
""")
    script.close()


def generate_script_afterburner(folder_name):
    """This function generates script for hadronic afterburner"""
    working_folder = folder_name

    script = open(path.join(working_folder, "run_afterburner.sh"), "w")
    script.write(
        """#!/usr/bin/env bash

unalias ls

SubEventId=$1

(
cd UrQMDev_$SubEventId

mkdir -p UrQMD_results
rm -fr UrQMD_results/*

for iev in `ls hydro_event | grep "surface"`
do
    cd iSS
    mkdir -p results
    rm -fr results/*
    mv ../hydro_event/$iev results/surface.dat
    mv ../hydro_event/music_input results/music_input
    ./iSS.e
    # turn on global momentum conservation
    #./correct_momentum_conservation.py OSCAR.DAT
    #mv OSCAR_w_GMC.DAT OSCAR.DAT
    cd ../osc2u
    ./osc2u.e < ../iSS/OSCAR.DAT
    mv fort.14 ../urqmd/OSCAR.input
    cd ../urqmd
    ./runqmd.sh
    mv particle_list.dat ../UrQMD_results/particle_list.dat
    rm -fr ../iSS/OSCAR.DAT
    rm -fr OSCAR.input
    cd ..
    ../hadronic_afterburner_toolkit/convert_to_binary.e UrQMD_results/particle_list.dat
    rm -fr UrQMD_results/particle_list.dat
done

rm -fr hydro_event
)
""")
    script.close()


def generate_script_analyze_spvn(folder_name):
    """This function generates script for analysis"""
    working_folder = folder_name

    script = open(path.join(working_folder, "run_analysis_spvn.sh"), "w")
    script.write(
        """#!/usr/bin/env bash

pid=$1

(
    cd hadronic_afterburner_toolkit
    if [ "$pid" == "9999" ]; then
        # charged hadrons
        ./hadronic_afterburner_tools.e run_mode=0 read_in_mode=2 particle_monval=$pid resonance_feed_down_flag=0 distinguish_isospin=0 rap_type=0 rap_min=-0.5 rap_max=0.5 compute_correlation=0 flag_charge_dependence=0
        ./hadronic_afterburner_tools.e run_mode=0 read_in_mode=2 particle_monval=$pid resonance_feed_down_flag=0 distinguish_isospin=0 rap_type=0 rap_min=-1.0 rap_max=-0.1 compute_correlation=0 flag_charge_dependence=0
        ./hadronic_afterburner_tools.e run_mode=0 read_in_mode=2 particle_monval=$pid resonance_feed_down_flag=0 distinguish_isospin=0 rap_type=0 rap_min=0.1 rap_max=1.0 compute_correlation=0 flag_charge_dependence=0
        ./hadronic_afterburner_tools.e run_mode=0 read_in_mode=2 particle_monval=$pid resonance_feed_down_flag=0 distinguish_isospin=0 rap_type=0 rap_min=0.5 rap_max=2.0 compute_correlation=0 flag_charge_dependence=0
        ./hadronic_afterburner_tools.e run_mode=0 read_in_mode=2 particle_monval=$pid resonance_feed_down_flag=0 distinguish_isospin=0 rap_type=0 rap_min=-2.0 rap_max=-0.5 compute_correlation=0 flag_charge_dependence=0
        ./hadronic_afterburner_tools.e run_mode=0 read_in_mode=2 particle_monval=$pid resonance_feed_down_flag=0 distinguish_isospin=0 rap_type=0 rap_min=-1.0 rap_max=1.0 compute_correlation=1 flag_charge_dependence=1 pT_min=0.2 pT_max=2.0
        ./hadronic_afterburner_tools.e run_mode=0 read_in_mode=2 particle_monval=$pid resonance_feed_down_flag=0 distinguish_isospin=0 rap_type=0 rap_min=-2.0 rap_max=2.0 compute_correlation=1 flag_charge_dependence=1 pT_min=0.2 pT_max=2.0
        ./hadronic_afterburner_tools.e run_mode=0 read_in_mode=2 particle_monval=$pid resonance_feed_down_flag=0 distinguish_isospin=0 rap_type=0 rap_min=-1.0 rap_max=1.0 compute_correlation=0 flag_charge_dependence=0
        ./hadronic_afterburner_tools.e run_mode=0 read_in_mode=2 particle_monval=$pid resonance_feed_down_flag=0 distinguish_isospin=0 rap_type=0 rap_min=-2.0 rap_max=2.0 compute_correlation=0 flag_charge_dependence=0
    else
        #./hadronic_afterburner_tools.e run_mode=0 read_in_mode=2 particle_monval=$pid resonance_feed_down_flag=0 distinguish_isospin=1 rap_type=0 rap_min=-0.5 rap_max=0.5 compute_correlation=0 flag_charge_dependence=0
        ./hadronic_afterburner_tools.e run_mode=0 read_in_mode=2 particle_monval=$pid resonance_feed_down_flag=0 distinguish_isospin=1 rap_type=1 rap_min=-0.5 rap_max=0.5 compute_correlation=0 flag_charge_dependence=0
    fi
)
""")
    script.close()


def generate_event_folders(initial_condition_database,
                           initial_condition_type, working_folder,
                           cluster_name, event_id,
                           n_hydro_per_job, n_urqmd_per_hydro, n_threads):
    """This function creates the event folder structure"""
    event_folder = path.join(working_folder, 'event_%d' % event_id)
    mkdir(event_folder)
    shutil.copy('codes/hydro_plus_UrQMD_driver.py', event_folder)
    shutil.copy(path.join('IPGlasma_database',
                          'fetch_IPGlasma_event_from_hdf5_database.py'),
                event_folder)
    shutil.copy(path.join('3DMCGlauber_database',
                          'fetch_3DMCGlauber_event_from_hdf5_database.py'),
                event_folder)
    if initial_condition_database == "self":
        shutil.copytree('codes/3dMCGlauber',
                        path.join(event_folder, '3dMCGlauber'), symlinks=True)
        subprocess.call("ln -s {0:s} {1:s}".format(
                        path.abspath('codes/3dMCGlauber_code/3dMCGlb.e'),
                        path.join(event_folder, "3dMCGlauber/3dMCGlb.e")),
                        shell=True)
        subprocess.call("ln -s {0:s} {1:s}".format(
                        path.abspath('codes/3dMCGlauber_code/eps09'),
                        path.join(event_folder, "3dMCGlauber/eps09")),
                        shell=True)

    generate_full_job_script(cluster_name, event_folder,
                             initial_condition_database,
                             initial_condition_type, n_hydro_per_job,
                             event_id*n_hydro_per_job, n_urqmd_per_hydro,
                             n_threads)

    generate_script_hydro(event_folder, n_threads)

    shutil.copytree('codes/MUSIC', path.join(event_folder, 'MUSIC'),
                    symlinks=True)
    subprocess.call("ln -s {0:s} {1:s}".format(
        path.abspath('codes/MUSIC_code/EOS'),
        path.join(event_folder, "MUSIC/EOS")), shell=True)
    subprocess.call("ln -s {0:s} {1:s}".format(
        path.abspath('codes/MUSIC_code/MUSIChydro'),
        path.join(event_folder, "MUSIC/MUSIChydro")), shell=True)
    generate_script_afterburner(event_folder)
    generate_script_analyze_spvn(event_folder)
    for iev in range(n_urqmd_per_hydro):
        sub_event_folder = path.join(working_folder,
                                     'event_{0:d}'.format(event_id),
                                     'UrQMDev_{0:d}'.format(iev))
        mkdir(sub_event_folder)
        shutil.copytree('codes/iSS', path.join(sub_event_folder, 'iSS'))
        subprocess.call("ln -s {0:s} {1:s}".format(
            path.abspath('codes/iSS_code/iSS_tables'),
            path.join(sub_event_folder, "iSS/iSS_tables")), shell=True)
        subprocess.call("ln -s {0:s} {1:s}".format(
            path.abspath('codes/iSS_code/iSS.e'),
            path.join(sub_event_folder, "iSS/iSS.e")), shell=True)
        shutil.copytree('codes/osc2u', path.join(sub_event_folder, 'osc2u'))
        shutil.copytree('codes/urqmd', path.join(sub_event_folder, 'urqmd'))
        subprocess.call("ln -s {0:s} {1:s}".format(
            path.abspath('codes/urqmd_code/urqmd/urqmd.e'),
            path.join(sub_event_folder, "urqmd/urqmd.e")), shell=True)
    shutil.copytree('codes/hadronic_afterburner_toolkit',
                    path.join(event_folder, 'hadronic_afterburner_toolkit'))
    subprocess.call("ln -s {0:s} {1:s}".format(
        path.abspath(path.join('codes', 'hadronic_afterburner_toolkit_code',
                               'hadronic_afterburner_tools.e')),
        path.join(event_folder, "hadronic_afterburner_toolkit",
                  "hadronic_afterburner_tools.e")), shell=True)
    subprocess.call("ln -s {0:s} {1:s}".format(
        path.abspath('codes/hadronic_afterburner_toolkit_code/EOS'),
        path.join(event_folder, "hadronic_afterburner_toolkit/EOS")),
                    shell=True)


def print_usage():
    """This function prints out help message"""
    print("\U0000269B  Usage: {} ".format(sys.argv[0])
          + "initial_condition_type initial_condition_filename working_folder "
          + "cluster_name n_jobs n_hydro_per_job n_urqmd_per_hydro [n_threads]")
    print("\U0000269B  initial_condition_type: IPGlasma, 3DMCGlauber")
    print("\U0000269B  cluster_name: nersc, wsugrid, local, guillimin, McGill")

def main():
    """This is the main funciton"""
    try:
        initial_condition_type = str(sys.argv[1])
        initial_condition_database = str(sys.argv[2])
        working_folder_name = str(sys.argv[3])
        cluster_name = str(sys.argv[4])
        n_jobs = int(sys.argv[5])
        n_hydro_per_job = int(sys.argv[6])
        n_urqmd_per_hydro = int(sys.argv[7])
    except IndexError:
        print_usage()
        exit(0)

    if len(sys.argv) > 8:
        n_threads = int(sys.argv[8])
    else:
        n_threads = n_urqmd_per_hydro

    if n_threads < n_urqmd_per_hydro:
        print("\U000026A0  "
              + "Warning: n_threads = {} < n_urqmd_per_hydro = {}!".format(
                  n_threads, n_urqmd_per_hydro))
        print("reset n_UrQMD to n_threads.")
        n_urqmd_per_hydro = n_threads

    if initial_condition_type not in ("IPGlasma", "3DMCGlauber"):
        print("\U0001F6AB  "
              + "Do not recognize the initial condition type: {}".format(
                  initial_condition_type))
        exit(1)

    subprocess.call("(cd config; python3 parameters_dict_master.py;)",
                    shell=True)

    if initial_condition_database == "self":
        print("\U0001F375  Generate initial condition on the fly ... ")
    else:
        initial_condition_database = path.abspath(initial_condition_database)
        print("\U0001F375  "
              + "Pre-generated initial conditions from {} ...".format(
                  initial_condition_database))
    working_folder_name = path.abspath(working_folder_name)
    mkdir(working_folder_name)

    toolbar_width = 40
    sys.stdout.write("\U0001F375  Generating {} jobs [{}]".format(
        n_jobs, " " * toolbar_width))
    sys.stdout.flush()
    sys.stdout.write("\b" * (toolbar_width+1))
    for iev in range(n_jobs):
        progress_i = (int(float(iev + 1)/n_jobs*toolbar_width)
                      - int(float(iev)/n_jobs*toolbar_width))
        for ii in range(progress_i):
            sys.stdout.write("#")
            sys.stdout.flush()
        generate_event_folders(initial_condition_database,
                               initial_condition_type, working_folder_name,
                               cluster_name, iev,
                               n_hydro_per_job, n_urqmd_per_hydro, n_threads)
    sys.stdout.write("\n")
    sys.stdout.flush()
    # copy script to collect final results
    pwd = path.abspath(".")
    script_path = "utilities"
    shutil.copy(path.join(script_path, 'collect_events.sh'), pwd)
    shutil.copy(path.join(script_path, 'combine_results_into_hdf5.py'), pwd)
    script_path = "codes/hadronic_afterburner_toolkit_code/ebe_scripts"
    shutil.copy(path.join(script_path, 'average_event_spvn_h5.py'), pwd)

    if cluster_name == "nersc":
        shutil.copy('Cluster_supports/NERSC/job_MPI_wrapper.py',
                    working_folder_name)
        n_nodes = max(1, int(n_jobs*n_threads/64))
        generate_nersc_mpi_job_script(working_folder_name,
                                      n_nodes, n_threads, int(n_jobs/n_nodes))

    if cluster_name == "nerscKNL":
        shutil.copy('Cluster_supports/NERSC/job_MPI_wrapper.py',
                    working_folder_name)
        n_nodes = max(1, int(n_jobs*n_threads/136))
        generate_nerscKNL_mpi_job_script(working_folder_name,
                                         n_nodes, n_threads,
                                         int(n_jobs/n_nodes))

    if cluster_name == "wsugrid":
        shutil.copy('Cluster_supports/WSUgrid/submit_all_jobs.sh', pwd)

if __name__ == "__main__":
    main()
