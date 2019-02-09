#!/usr/bin/env python

import sys
from os import path, mkdir
import shutil
from glob import glob
import subprocess
import random
from IPGlasma_database.fetch_IPGlasma_event_from_hdf5_database import fecth_an_IPGlasma_event


def write_script_header(cluster, script, event_id, walltime, working_folder):
    if cluster == "nersc":
        script.write(
"""#!/bin/bash -l
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -J event_{0:s}
#SBATCH -t {1:s}
#SBATCH -L SCRATCH
#SBATCH -C haswell
""".format(event_id, walltime))
    elif cluster == "guillimin":
        script.write(
"""#!/usr/bin/env bash
#PBS -N event_{0:s}
#PBS -l nodes=1:ppn=1
#PBS -l walltime={1:s}
#PBS -S /bin/bash
#PBS -e test.err
#PBS -o test.log
#PBS -A cqn-654-ad
#PBS -q sw
#PBS -d {2:s}
""".format(event_id, walltime, working_folder))
    elif cluster == "McGill":
        script.write(
"""#!/usr/bin/env bash
#PBS -N event_{0:s}
#PBS -l nodes=1:ppn=1:irulan
#PBS -l walltime={1:s}
#PBS -S /bin/bash
#PBS -e test.err
#PBS -o test.log
#PBS -d {2:s}
""".format(event_id, walltime, working_folder))
    elif cluster == "wsugrid":
        script.write(
"""#!/usr/bin/env bash
#PBS -N event_{0:s}
#PBS -l select=1:ncpus=1:mem=5GB
#PBS -l walltime={1:s}
#PBS -S /bin/bash
#PBS -e test.err
#PBS -o test.log
#PBS -q wsuq

cd {2:s}
""".format(event_id, walltime, working_folder))
    elif cluster == "local":
        script.write("#!/usr/bin/env bash")
    else:
        print("Error: unrecoginzed cluster name :", cluster)
        print("Available options: nersc, wsugrid, local, guillimin, McGill")
        exit(1)


def generate_script_hydro(cluster_name, folder_name, nthreads):
    working_folder = folder_name
    event_id = working_folder.split('/')[-1]
    walltime = '35:00:00'

    script = open(path.join(working_folder, "run_hydro.sh"), "w")

    hydro_results_folder = 'hydro_results'
    ppn = nthreads
    script.write(
"""#!/usr/bin/env bash

results_folder={0:s}

export OMP_NUM_THREADS={1:d}
export OMP_PROC_BIND=spread
export OMP_PLACES=threads

(
cd MUSIC
# hydro evolution
./mpihydro music_input_mode_2 1> mode_2.log 2> mode_2.err
./sweeper.sh $results_folder
)
""".format(hydro_results_folder, ppn))
    script.close()


def generate_script_afterburner(cluster_name, folder_name):
    working_folder = folder_name
    event_id = working_folder.split('/')[-1]
    walltime = '35:00:00'

    script = open(path.join(working_folder, "run_afterburner.sh"), "w")
    script.write(
"""#!/usr/bin/env bash

SubEventId=$1

(
cd UrQMDev_$SubEventId

mkdir UrQMD_results
for iev in `ls hydro_events --color=none | grep "surface"`
do
    event_id=`echo $iev | cut -f 3 -d _ | cut -f 1 -d .`
    cd iSS
    if [ -d "results" ]; then
        rm -fr results
    fi
    mkdir results
    mv ../hydro_events/$iev results/surface.dat
    cp ../hydro_events/music_input_event_$event_id results/music_input
    ./iSS.e >> ../output.log
    mv results/surface.dat ../hydro_events/$iev
    #rm -fr results/sample*
    # turn on global momentum conservation
    ./correct_momentum_conservation.py OSCAR.DAT
    mv OSCAR_w_GMC.DAT OSCAR.DAT
    cd ../osc2u
    ./osc2u.e < ../iSS/OSCAR.DAT >> ../output.log
    mv fort.14 ../urqmd/OSCAR.input
    cd ../urqmd
    ./runqmd.sh >> ../output.log
    mv particle_list.dat ../UrQMD_results/particle_list_$event_id.dat
    rm -fr ../iSS/OSCAR.DAT
    rm -fr OSCAR.input
    cd ..
    ../hadronic_afterburner_toolkit/convert_to_binary.e UrQMD_results/particle_list_$event_id.dat
    rm -fr UrQMD_results/particle_list_$event_id.dat
done
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
    initial_folder = path.join(event_folder, 'Initial')
    mkdir(initial_folder)
    for iev in range(event_id*n_hydro_per_job, (event_id + 1)*n_hydro_per_job):
        copy_IPGlasma_initial_condition(initial_condition_database, iev,
                                        initial_folder)
    generate_script_hydro(cluster_name, event_folder, n_UrQMD_per_hydro)
    shutil.copytree('codes/MUSIC', path.join(event_folder, 'MUSIC'))
    generate_script_afterburner(cluster_name, event_folder)
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

    working_folder_name = path.abspath(working_folder_name)
    mkdir(working_folder_name)
    for iev in range(n_jobs):
        generate_event_folders(initial_condition_database, working_folder_name,
                               cluster_name, iev,
                               n_hydro_per_job, n_UrQMD_per_hydro)

