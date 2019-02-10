#!/usr/bin/env python

import sys
from os import path, mkdir
import shutil
from glob import glob
import subprocess
import random


def write_script_header(cluster, script, n_threads,
                        event_id, walltime, working_folder):
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
#PBS -l select=1:ncpus={1:d}:mem=5GB
#PBS -l walltime={2:s}
#PBS -S /bin/bash
#PBS -e test.err
#PBS -o test.log
#PBS -q wsuq

cd {3:s}
""".format(event_id, n_threads, walltime, working_folder))
    elif cluster == "local":
        script.write("#!/usr/bin/env bash")
    else:
        print("Error: unrecoginzed cluster name :", cluster)
        print("Available options: nersc, wsugrid, local, guillimin, McGill")
        exit(1)


def generate_full_job_script(cluster_name, folder_name,
                             database, n_hydro, ev0_id, n_threads):
    working_folder = folder_name
    event_id = working_folder.split('/')[-1]
    walltime = '35:00:00'

    script = open(path.join(working_folder, "submit_job.pbs"), "w")
    write_script_header(cluster_name, script, n_threads, event_id, walltime,
                        working_folder)
    script.write(
"""
./hydro_plus_UrQMD_driver.py {0:s} {1:d} {2:d} {3:d} > run.log
""".format(database, n_hydro, ev0_id, n_threads))
    script.close()


def generate_script_hydro(folder_name, nthreads):
    working_folder = folder_name

    script = open(path.join(working_folder, "run_hydro.sh"), "w")

    hydro_results_folder = 'hydro_results'
    ppn = nthreads
    script.write(
"""#!/usr/bin/env bash

results_folder={0:s}

(
cd MUSIC

export OMP_NUM_THREADS={1:d}

# hydro evolution
./mpihydro music_input_mode_2
./sweeper.sh $results_folder
)
""".format(hydro_results_folder, ppn))
    script.close()


def generate_script_afterburner(folder_name):
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
    working_folder = folder_name

    script = open(path.join(working_folder, "run_analysis_spvn.sh"), "w")
    script.write(
"""#!/usr/bin/env bash

pid=$1

(
    cd hadronic_afterburner_toolkit
    if [ "$pid" == "9999" ]; then
        # charged hadrons
        ./hadronic_afterburner_tools.e run_mode=0 read_in_mode=2 particle_monval=$pid resonance_feed_down_flag=0 distinguish_isospin=0 rap_type=0
        ./hadronic_afterburner_tools.e run_mode=0 read_in_mode=2 particle_monval=$pid resonance_feed_down_flag=0 distinguish_isospin=0 rap_type=0 rap_min=-1.0 rap_max=-0.1
        ./hadronic_afterburner_tools.e run_mode=0 read_in_mode=2 particle_monval=$pid resonance_feed_down_flag=0 distinguish_isospin=0 rap_type=0 rap_min=0.1 rap_max=1.0
        ./hadronic_afterburner_tools.e run_mode=0 read_in_mode=2 particle_monval=$pid resonance_feed_down_flag=0 distinguish_isospin=0 rap_type=0 rap_min=0.5 rap_max=2.0
        ./hadronic_afterburner_tools.e run_mode=0 read_in_mode=2 particle_monval=$pid resonance_feed_down_flag=0 distinguish_isospin=0 rap_type=0 rap_min=-2.0 rap_max=-0.5
        ./hadronic_afterburner_tools.e run_mode=0 read_in_mode=2 particle_monval=$pid resonance_feed_down_flag=0 distinguish_isospin=0 rap_type=0 rap_min=-1.0 rap_max=1.0 compute_correlation=1 flag_charge_dependence=1 pT_min=0.2 pT_max=2.0
        ./hadronic_afterburner_tools.e run_mode=0 read_in_mode=2 particle_monval=$pid resonance_feed_down_flag=0 distinguish_isospin=0 rap_type=0 rap_min=-2.0 rap_max=2.0 compute_correlation=1 flag_charge_dependence=1 pT_min=0.2 pT_max=2.0
        ./hadronic_afterburner_tools.e run_mode=0 read_in_mode=2 particle_monval=$pid resonance_feed_down_flag=0 distinguish_isospin=0 rap_type=0 rap_min=-1.0 rap_max=1.0
        ./hadronic_afterburner_tools.e run_mode=0 read_in_mode=2 particle_monval=$pid resonance_feed_down_flag=0 distinguish_isospin=0 rap_type=0 rap_min=-2.0 rap_max=2.0
    else
        #./hadronic_afterburner_tools.e run_mode=0 read_in_mode=2 particle_monval=$pid resonance_feed_down_flag=0 distinguish_isospin=1 rap_type=0
        ./hadronic_afterburner_tools.e run_mode=0 read_in_mode=2 particle_monval=$pid resonance_feed_down_flag=0 distinguish_isospin=1 rap_type=1
    fi
)
""")
    script.close()


def copy_IPGlasma_initial_condition(database, event_id, folder):
    time_stamp_str = "0.4"
    file_name = fecth_an_IPGlasma_event(database, time_stamp_str, event_id)
    shutil.move(file_name, folder)

        
def generate_event_folders(initial_condition_database, working_folder,
                           cluster_name, event_id,
                           n_hydro_per_job, n_UrQMD_per_hydro):
    event_folder = path.join(working_folder, 'event_%d' % event_id)
    mkdir(event_folder)
    shutil.copy('codes/hydro_plus_UrQMD_driver.py', event_folder)
    shutil.copy('IPGlasma_database/fetch_IPGlasma_event_from_hdf5_database.py',
                event_folder)
    generate_full_job_script(cluster_name, event_folder,
                             initial_condition_database, n_hydro_per_job,
                             event_id*n_hydro_per_job, n_UrQMD_per_hydro)
    generate_script_hydro(event_folder, n_UrQMD_per_hydro)
    shutil.copytree('codes/MUSIC', path.join(event_folder, 'MUSIC'))
    generate_script_afterburner(event_folder)
    generate_script_analyze_spvn(event_folder)
    for iev in range(n_UrQMD_per_hydro):
        sub_event_folder = path.join(working_folder,
                                     'event_{0:d}'.format(event_id),
                                     'UrQMDev_{0:d}'.format(iev))
        mkdir(sub_event_folder)
        shutil.copytree('codes/iSS',   path.join(sub_event_folder, 'iSS'  ))
        shutil.copytree('codes/osc2u', path.join(sub_event_folder, 'osc2u'))
        shutil.copytree('codes/urqmd', path.join(sub_event_folder, 'urqmd'))
    shutil.copytree('codes/hadronic_afterburner_toolkit', 
                    path.join(event_folder, 'hadronic_afterburner_toolkit'))

    if cluster_name == "nersc":
        shutil.copy('NERSC_supports/job_MPI_wrapper.py', working_folder)
        shutil.copy('NERSC_supports/submit_MPI_jobs.pbs', working_folder)


def print_Usage():
    print("Usage: {} initial_condition working_folder ".format(sys.argv[0])
          + "cluster_name n_jobs n_hydro_per_job n_UrQMD_per_hydro")
                    
if __name__ == "__main__":
    try:
        initial_condition_database = str(sys.argv[1])
        working_folder_name        = str(sys.argv[2])
        cluster_name               = str(sys.argv[3])
        n_jobs                     = int(sys.argv[4])
        n_hydro_per_job            = int(sys.argv[5])
        n_UrQMD_per_hydro          = int(sys.argv[6])
    except IndexError:
        print_Usage()
        exit(0)

    initial_condition_database = path.abspath(initial_condition_database)
    working_folder_name        = path.abspath(working_folder_name)
    mkdir(working_folder_name)
    for iev in range(n_jobs):
        generate_event_folders(initial_condition_database, working_folder_name,
                               cluster_name, iev,
                               n_hydro_per_job, n_UrQMD_per_hydro)

