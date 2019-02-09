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

def generate_script(cluster_name, folder_name):
    working_folder = path.join(path.abspath('./'), folder_name)
    event_id = working_folder.split('/')[-1]
    walltime = '35:00:00'

    script = open(path.join(working_folder, "submit_job.pbs"), "w")
    write_script_header(cluster_name, script, event_id, walltime,
                        working_folder)
    script.close()

def copy_IPGlasma_initial_condition(database, event_id, folder):
    time_stamp_str = "0.4"
    file_name = fecth_an_IPGlasma_event(database, time_stamp_str, event_id)
    shutil.move(file_name, folder)

        
def generate_event_folders(initial_condition_database, working_folder,
                           cluster_name, event_id, nsubev):
    event_folder = path.join(working_folder, 'event_%d' % event_id)
    mkdir(event_folder)
    generate_script(cluster_name, event_folder)
    shutil.copytree('codes/MUSIC', 
                    path.join(path.abspath(event_folder), 'MUSIC'))
    copy_IPGlasma_initial_condition(initial_condition_database, event_id, 
                                    path.join(path.abspath(event_folder),
                                              'MUSIC', 'initial'))
    for iev in range(nsubev):
        sub_event_folder = path.join(working_folder,
                                     'event_{0:d}'.format(event_id),
                                     'subev_{0:d}'.format(iev))
        mkdir(sub_event_folder)
        shutil.copytree('codes/iSS', 
                        path.join(path.abspath(sub_event_folder), 'iSS'))
        shutil.copytree('codes/osc2u', 
                        path.join(path.abspath(sub_event_folder), 'osc2u'))
        shutil.copytree('codes/urqmd', 
                        path.join(path.abspath(sub_event_folder), 'urqmd'))
    shutil.copytree('codes/hadronic_afterburner_toolkit', 
                    path.join(path.abspath(event_folder), 
                    'hadronic_afterburner_toolkit'))

def print_Usage():
    print("Usage: {} initial_condition working_folder ".format(sys.argv[0])
          + "cluster_name n_hydro_ev n_UrQMD_per_hydro")
                    
if __name__ == "__main__":
    try:
        initial_condition_database = str(sys.argv[1])
        working_folder_name        = str(sys.argv[2])
        cluster_name               = str(sys.argv[3])
        n_hydro_ev                 = int(sys.argv[4])
        n_UrQMD_per_hydro          = int(sys.argv[5])
    except IndexError:
        print_Usage()
        exit(0)

    mkdir(working_folder_name)
    for iev in range(n_hydro_ev):
        generate_event_folders(initial_condition_database, working_folder_name,
                               cluster_name, iev, n_UrQMD_per_hydro)

