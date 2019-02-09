#!/usr/bin/env python

import sys
from os import path, mkdir
import shutil
from glob import glob
import subprocess
import random


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
        
def generate_event_folders(cluster_name, working_folder, event_id, nsubev):
    event_folder = path.join(working_folder, 'event_%d' % event_id)
    mkdir(event_folder)
    generate_script(cluster_name, event_folder)
    shutil.copytree('codes/MUSIC', 
                    path.join(path.abspath(event_folder), 'MUSIC'))
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
    print("Usage:")
    print("  %s input_folder working_folder cluster_name num_of_cores nsubev"
          % str(sys.argv[0]))
    print("")
                    
if __name__ == "__main__":
    try:
        from_folder = str(sys.argv[1])
        folder_name = str(sys.argv[2])
        cluster_name = str(sys.argv[3])
        ncore = int(sys.argv[4])
        nsubev = int(sys.argv[5])
    except IndexError:
        print_Usage()
        exit(0)

    mkdir(folder_name)
    for icore in range(ncore):
        generate_event_folders(cluster_name, folder_name, icore, nsubev)
