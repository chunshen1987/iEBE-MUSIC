#!/usr/bin/env python3
"""This script generate all the running jobs."""

import sys
import re
from os import path, mkdir
import shutil
import subprocess
import argparse
from math import ceil
from glob import glob

centrality_list = [(0.00, 0.15, '0-5', 0.05), (0.15, 0.30, '5-10', 0.05),
                   (0.30, 0.45, '10-20', 0.10), (0.45, 0.55, '20-30', 0.10),
                   (0.55, 0.65, '30-40', 0.10), (0.65, 0.75, '40-50', 0.10),
                   (0.75, 0.80, '50-60', 0.10), (0.80, 0.85, '60-70', 0.10),
                   (0.85, 0.90, '70-80', 0.10), (0.90, 0.95, '80-90', 0.10),
                   (0.95, 1.00, '90-100', 0.10)]

known_initial_types = [
    "IPGlasma", "IPGlasma+KoMPoST", "3DMCGlauber_dynamical",
    "3DMCGlauber_consttau"
]

support_cluster_list = [
    'nersc', 'nerscKNL', 'wsugrid', "OSG", "local", "guillimin", "McGill"
]


def write_script_header(cluster, script, n_threads, event_id, walltime,
                        working_folder):
    """This function write the header of the job submission script"""
    mem = 4*n_threads
    if cluster == "nersc":
        script.write("""#!/bin/bash -l
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -J {0:s}
#SBATCH -t {1:s}
#SBATCH -L SCRATCH
#SBATCH -C haswell
""".format(event_id, walltime))
    elif cluster == "nerscKNL":
        script.write("""#!/bin/bash -l
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -J {0:s}
#SBATCH -t {1:s}
#SBATCH -L SCRATCH
#SBATCH -C knl,quad,cache
""".format(event_id, walltime))
    elif cluster == "guillimin":
        script.write("""#!/usr/bin/env bash
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
        script.write("""#!/usr/bin/env bash
#PBS -N {0:s}
#PBS -l nodes=1:ppn={1:d}:irulan
#PBS -l walltime={2:s}
#PBS -S /bin/bash
#PBS -e test.err
#PBS -o test.log
#PBS -d {3:s}
""".format(event_id, n_threads, walltime, working_folder))
    elif cluster == "wsugrid":
        script.write("""#!/usr/bin/env bash
#SBATCH --job-name {0:s}
#SBATCH -q primary
#SBATCH -N 1
#SBATCH -n {1:d}
#SBATCH --mem={2:.0f}G
#SBATCH --constraint=intel
#SBATCH -t {3:s}
#SBATCH -e job.err
#SBATCH -o job.log

cd {4:s}
""".format(event_id, n_threads, mem, walltime, working_folder))
    elif cluster in ("local", "OSG"):
        script.write("#!/bin/bash")
    else:
        print("\U0001F6AB  unrecoginzed cluster name :", cluster)
        print("Available options: ", support_cluster_list)
        exit(1)


def generate_nersc_mpi_job_script(folder_name, n_nodes, n_threads,
                                  n_jobs_per_node, walltime):
    """This function generates job script for NERSC"""
    working_folder = folder_name

    script = open(path.join(working_folder, "submit_MPI_jobs.pbs"), "w")
    script.write("""#!/bin/bash -l
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
                                     n_jobs_per_node, walltime):
    """This function generates job script for NERSC KNL"""
    working_folder = folder_name

    script = open(path.join(working_folder, "submit_MPI_jobs.pbs"), "w")
    script.write("""#!/bin/bash -l
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
                             n_hydro, ev0_id, n_urqmd, n_threads, para_dict,
                             time_stamp):
    """This function generates full job script"""
    working_folder = folder_name
    event_id = working_folder.split('/')[-1]
    walltime = '100:00:00'

    script = open(path.join(working_folder, "submit_job.pbs"), "w")
    write_script_header(cluster_name, script, n_threads, event_id, walltime,
                        working_folder)
    script.write("\nseed_add=${1:-0}\n")
    script.write("""
python3 hydro_plus_UrQMD_driver.py {0:s} {1:s} {2:d} {3:d} {4:d} {5:d} {6} {7} {8} {9} $seed_add {10:s} {11} {12}
""".format(initial_type, database, n_hydro, ev0_id, n_urqmd, n_threads,
           para_dict.control_dict["save_ipglasma_results"],
           para_dict.control_dict["save_kompost_results"],
           para_dict.control_dict["save_hydro_surfaces"],
           para_dict.control_dict["save_UrQMD_files"],
           time_stamp,
           para_dict.control_dict["compute_polarization"],
           para_dict.control_dict["compute_photon_emission"]))
    script.write("""

status=$?
if [ $status -ne 0 ]; then
    exit $status
fi""")
    script.close()


def generate_script_ipglasma(folder_name, nthreads, cluster_name, event_id):
    """This function generates script for IPGlasma simulation"""
    working_folder = folder_name

    script = open(path.join(working_folder, "run_ipglasma.sh"), "w")

    results_folder = 'ipglasma_results'
    script.write("""#!/bin/bash

results_folder={0:s}
evid=$1

(
cd ipglasma

mkdir -p $results_folder
rm -fr $results_folder/*

""".format(results_folder))

    if nthreads > 0:
        script.write("""
export OMP_NUM_THREADS={0:d}
""".format(nthreads))

    if cluster_name != "OSG":
        script.write("sleep {}".format(event_id))
        script.write("""
# IPGlasma evolution (run 1 event)
./ipglasma input 1> run.log 2> run.err
""")
    else:
        script.write("""
# IPGlasma evolution (run 1 event)
./ipglasma input
""")
    script.write("""
for ifile in *.dat
do
    filename=$(echo ${ifile} | sed "s/0.dat/${evid}.dat/")
    cat ${ifile} | sed 's$N/A$0.0$g' | sed 's/Q_s/#Q_s/' > $results_folder/${filename}
    rm -fr ${ifile}
done
mv run.log $results_folder/
mv run.err $results_folder/
)
""")
    script.close()


def generate_script_kompost(folder_name, nthreads, cluster_name):
    """This function generates script for KoMPoST simulation"""
    working_folder = folder_name

    script = open(path.join(working_folder, "run_kompost.sh"), "w")

    hydro_results_folder = 'kompost_results'
    script.write("""#!/bin/bash

results_folder={0:s}

(
cd kompost

mkdir -p $results_folder
rm -fr $results_folder/*

""".format(hydro_results_folder))

    if nthreads > 0:
        script.write("""
export OMP_NUM_THREADS={0:d}
""".format(nthreads))

    if cluster_name != "OSG":
        script.write("""
# KoMPoST EKT evolution
./KoMPoST.exe setup.ini 1> run.log 2> run.err
mv *.txt $results_folder
)
""")
    else:
        script.write("""
# KoMPoST EKT evolution
./KoMPoST.exe setup.ini
mv *.txt $results_folder
)
""")
    script.close()


def generate_script_hydro(folder_name, nthreads, cluster_name):
    """This function generates script for hydro simulation"""
    working_folder = folder_name

    script = open(path.join(working_folder, "run_hydro.sh"), "w")

    hydro_results_folder = 'hydro_results'
    script.write("""#!/bin/bash

results_folder={0:s}

(
cd MUSIC

rm -fr *.dat
rm -fr $results_folder

""".format(hydro_results_folder))

    if nthreads > 0:
        script.write("""
export OMP_NUM_THREADS={0:d}
""".format(nthreads))

    if cluster_name != "OSG":
        script.write("""
# hydro evolution
./MUSIChydro music_input_mode_2 1> run.log 2> run.err
./sweeper.sh $results_folder
)
""")
    else:
        script.write("""
# hydro evolution
./MUSIChydro music_input_mode_2 | tee run.log
./sweeper.sh $results_folder
)
""")
    script.close()


def generate_script_photon(folder_name, nthreads, cluster_name):
    """This function generates script for photon radiation"""
    working_folder = folder_name

    script = open(path.join(working_folder, "run_photon.sh"), "w")

    script.write("""#!/bin/bash
(
cd photonEmission_hydroInterface

""")
    if nthreads > 0:
        script.write("""
export OMP_NUM_THREADS={0:d}
""".format(nthreads))

    if cluster_name != "OSG":
        script.write("""
# perform photon radiation
./hydro_photonEmission.e > run.log
)
""")
    else:
        script.write("""
# perform photon radiation
./hydro_photonEmission.e
)
""")
    script.close()


def generate_script_afterburner(folder_name, cluster_name, HBT_flag, GMC_flag):
    """This function generates script for hadronic afterburner"""
    working_folder = folder_name

    logfile = ""
    if cluster_name != "OSG":
        logfile = " >> run.log"

    script = open(path.join(working_folder, "run_afterburner.sh"), "w")
    script.write("""#!/bin/bash

unalias ls 2>/dev/null

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
    mv ../hydro_event/spectators.dat results/spectators.dat
    if [ $SubEventId = "0" ]; then
    """)
    script.write("    ./iSS.e {0}".format(logfile))
    script.write("""
    else
        ./iSS.e > run.log
    fi
    """)

    if GMC_flag == 1:
        script.write("""
    # turn on global momentum conservation
    ./correct_momentum_conservation.py OSCAR.DAT
    mv OSCAR_w_GMC.DAT OSCAR.DAT
    """)
    script.write("""
    cd ../osc2u
    if [ $SubEventId = "0" ]; then
    """)
    script.write("    ./osc2u.e < ../iSS/OSCAR.DAT {0}".format(logfile))
    script.write("""
    else
        ./osc2u.e < ../iSS/OSCAR.DAT >> run.log
    fi
    mv fort.14 ../urqmd/OSCAR.input
    rm -fr ../iSS/OSCAR.DAT
    cd ../urqmd
    ./runqmd.sh >> run.log
    mv particle_list.dat ../UrQMD_results/particle_list.dat
    rm -fr OSCAR.input
    cd ..
    ../hadronic_afterburner_toolkit/convert_to_binary.e UrQMD_results/particle_list.dat
    rm -fr UrQMD_results/particle_list.dat
""")
    if HBT_flag:
        script.write("""
    cd hadronic_afterburner_toolkit
    mkdir -p results
    cd results; rm -fr *
    ln -s ../../UrQMD_results/particle_list.gz particle_list.dat
    cd ..
""")
        script.write('    if [ $SubEventId = "0" ]; then\n')
        script.write(
            "        ./hadronic_afterburner_tools.e analyze_flow=0 analyze_HBT=1 particle_monval=211 distinguish_isospin=1 event_buffer_size=500000 {0}\n"
            .format(logfile))
        script.write("    else\n")
        script.write(
            "        ./hadronic_afterburner_tools.e analyze_flow=0 analyze_HBT=1 particle_monval=211 distinguish_isospin=1 event_buffer_size=500000 >> run.log\n"
        )
        script.write("    fi\n")
        script.write("    mv results/HBT* ../UrQMD_results/ \n")
    script.write("""
    done

    rm -fr hydro_event
)
""")
    script.close()


def generate_script_analyze_spvn(folder_name, cluster_name, HBT_flag):
    """This function generates script for analysis"""
    working_folder = folder_name

    logfile = ""
    if cluster_name != "OSG":
        logfile = " >> run.log"

    script = open(path.join(working_folder, "run_analysis_spvn.sh"), "w")
    script.write("""#!/bin/bash

(
    cd hadronic_afterburner_toolkit
""")
    script.write(
        "   ./hadronic_afterburner_tools.e analyze_HBT=0 {0}\n".format(logfile))
    if HBT_flag:
        script.write(
            "    python ./average_event_HBT_correlation_function.py .. results\n"
        )
    script.write(")\n")
    script.close()


def generate_event_folders(initial_condition_database, initial_condition_type,
                           package_root_path, code_path, working_folder,
                           cluster_name, event_id, event_id_offset,
                           n_hydro_per_job, n_urqmd_per_hydro, n_threads,
                           time_stamp, para_dict):
    """This function creates the event folder structure"""
    event_folder = path.join(working_folder, 'event_%d' % event_id)
    param_folder = path.join(working_folder, 'model_parameters')
    mkdir(event_folder)
    shutil.copy(path.join(code_path, 'hydro_plus_UrQMD_driver.py'),
                event_folder)
    shutil.copy(
        path.join(package_root_path, 'IPGlasma_database',
                  'fetch_IPGlasma_event_from_hdf5_database.py'), event_folder)
    shutil.copy(
        path.join(package_root_path, '3DMCGlauber_database',
                  'fetch_3DMCGlauber_event_from_hdf5_database.py'),
        event_folder)
    if initial_condition_database == "self":
        if initial_condition_type in ("3DMCGlauber_dynamical",
                                      "3DMCGlauber_consttau"):
            mkdir(path.join(event_folder, '3dMCGlauber'))
            shutil.copyfile(path.join(param_folder, '3dMCGlauber/input'),
                            path.join(event_folder, '3dMCGlauber/input'))
            for link_i in ['3dMCGlb.e', 'eps09', 'tables']:
                subprocess.call("ln -s {0:s} {1:s}".format(
                    path.abspath(
                        path.join(code_path,
                                  '3dMCGlauber_code/{}'.format(link_i))),
                    path.join(event_folder, "3dMCGlauber/{}".format(link_i))),
                                shell=True)
        elif initial_condition_type in ("IPGlasma", "IPGlasma+KoMPoST"):
            generate_script_ipglasma(event_folder, n_threads, cluster_name,
                                     event_id)
            mkdir(path.join(event_folder, 'ipglasma'))
            shutil.copyfile(path.join(param_folder, 'IPGlasma/input'),
                            path.join(event_folder, 'ipglasma/input'))
            link_list = [
                'qs2Adj_vs_Tp_vs_Y_200.in', 'utilities', 'ipglasma',
                'carbon_alpha_3.in', 'carbon_plaintext.in', 'oxygen_alpha_3.in',
                'oxygen_plaintext.in', 'he3_plaintext.in'
            ]
            for link_i in link_list:
                subprocess.call("ln -s {0:s} {1:s}".format(
                    path.abspath(
                        path.join(code_path,
                                  'ipglasma_code/{}'.format(link_i))),
                    path.join(event_folder, "ipglasma/{}".format(link_i))),
                                shell=True)

    generate_full_job_script(cluster_name, event_folder,
                             initial_condition_database,
                             initial_condition_type,
                             n_hydro_per_job, event_id_offset,
                             n_urqmd_per_hydro, n_threads, para_dict,
                             time_stamp)

    if initial_condition_type == "IPGlasma+KoMPoST":
        generate_script_kompost(event_folder, n_threads, cluster_name)
        mkdir(path.join(event_folder, 'kompost'))
        shutil.copyfile(path.join(param_folder, 'KoMPoST/setup.ini'),
                        path.join(event_folder, 'kompost/setup.ini'))
        for link_i in ['EKT', 'KoMPoST.exe']:
            subprocess.call("ln -s {0:s} {1:s}".format(
                path.abspath(
                    path.join(code_path, 'kompost_code/{}'.format(link_i))),
                path.join(event_folder, "kompost/{}".format(link_i))),
                            shell=True)

    # MUSIC
    generate_script_hydro(event_folder, n_threads, cluster_name)

    shutil.copytree(path.join(code_path, 'MUSIC'),
                    path.join(event_folder, 'MUSIC'))
    shutil.copyfile(path.join(param_folder, 'MUSIC/music_input_mode_2'),
                    path.join(event_folder, 'MUSIC/music_input_mode_2'))
    for link_i in ['EOS', 'MUSIChydro']:
        subprocess.call("ln -s {0:s} {1:s}".format(
            path.abspath(path.join(code_path, 'MUSIC_code/{}'.format(link_i))),
            path.join(event_folder, "MUSIC/{}".format(link_i))),
                        shell=True)

    if para_dict.control_dict['compute_photon_emission']:
        # photon
        generate_script_photon(event_folder, n_threads, cluster_name)
        mkdir(path.join(event_folder, 'photonEmission_hydroInterface'))
        shutil.copyfile(path.join(param_folder, 'photonEmission_hydroInterface',
                                  'parameters.dat'),
                        path.join(event_folder, 'photonEmission_hydroInterface',
                                  'parameters.dat'))
        for link_i in ['ph_rates', 'hydro_photonEmission.e']:
            orgFilePath = path.abspath(path.join(code_path,
                                       'photonEmission_hydroInterface_code',
                                       '{}'.format(link_i)))
            trgFilePath = path.join(event_folder,
                                    "photonEmission_hydroInterface",
                                    "{}".format(link_i))
            subprocess.call("ln -s {0:s} {1:s}".format(orgFilePath,
                                                       trgFilePath),
                            shell=True)

    # particlization + hadronic afterburner
    GMC_flag = para_dict.iss_dict['global_momentum_conservation']
    HBT_flag = False
    if "analyze_HBT" in para_dict.hadronic_afterburner_toolkit_dict:
        if para_dict.hadronic_afterburner_toolkit_dict['analyze_HBT'] == 1:
            HBT_flag = True

    generate_script_afterburner(event_folder, cluster_name, HBT_flag, GMC_flag)

    generate_script_analyze_spvn(event_folder, cluster_name, HBT_flag)

    nUrQMDFolders = n_urqmd_per_hydro
    if para_dict.control_dict['compute_polarization']:
        nUrQMDFolders += 1
    for iev in range(nUrQMDFolders):
        sub_event_folder = path.join(working_folder,
                                     'event_{}'.format(event_id),
                                     'UrQMDev_{}'.format(iev))
        mkdir(sub_event_folder)
        mkdir(path.join(sub_event_folder, 'iSS'))
        shutil.copyfile(path.join(param_folder, 'iSS/iSS_parameters.dat'),
                        path.join(sub_event_folder, 'iSS/iSS_parameters.dat'))
        if para_dict.control_dict['compute_polarization'] and iev < nUrQMDFolders:
            f1 = open("temp.dat", "w")
            with open(path.join(sub_event_folder, 'iSS/iSS_parameters.dat')) as f:
                for line in f:
                    line2 = re.sub("calculate_polarization = 1",
                                   "calculate_polarization = 0", line)
                    f1.write(line2)
            f1.close()
            shutil.copyfile("temp.dat", path.join(sub_event_folder,
                                                  'iSS/iSS_parameters.dat'))

        for link_i in ['iSS_tables', 'iSS.e']:
            subprocess.call("ln -s {0:s} {1:s}".format(
                path.abspath(path.join(code_path,
                                       'iSS_code/{}'.format(link_i))),
                path.join(sub_event_folder, "iSS/{}".format(link_i))),
                            shell=True)
        shutil.copytree(path.join(code_path, 'osc2u'),
                        path.join(sub_event_folder, 'osc2u'))
        shutil.copytree(path.join(code_path, 'urqmd'),
                        path.join(sub_event_folder, 'urqmd'))
        subprocess.call("ln -s {0:s} {1:s}".format(
            path.abspath(path.join(code_path, 'urqmd_code/urqmd/urqmd.e')),
            path.join(sub_event_folder, "urqmd/urqmd.e")),
                        shell=True)
        if HBT_flag:
            shutil.copytree(path.join(code_path,
                                      'hadronic_afterburner_toolkit'),
                            path.join(sub_event_folder,
                                      'hadronic_afterburner_toolkit'))
            shutil.copyfile(
                path.join(param_folder,
                          'hadronic_afterburner_toolkit/parameters.dat'),
                path.join(sub_event_folder,
                          'hadronic_afterburner_toolkit/parameters.dat'))
            for link_i in ['hadronic_afterburner_tools.e', 'EOS']:
                subprocess.call("ln -s {0:s} {1:s}".format(
                    path.abspath(
                        path.join(
                            code_path,
                            'hadronic_afterburner_toolkit_code/{}'.format(
                                link_i))),
                    path.join(
                        sub_event_folder,
                        "hadronic_afterburner_toolkit/{}".format(link_i))),
                                shell=True)
    shutil.copytree(path.join(code_path, 'hadronic_afterburner_toolkit'),
                    path.join(event_folder, 'hadronic_afterburner_toolkit'))
    shutil.copyfile(
        path.join(param_folder, 'hadronic_afterburner_toolkit/parameters.dat'),
        path.join(event_folder, 'hadronic_afterburner_toolkit/parameters.dat'))
    for link_i in ['hadronic_afterburner_tools.e', 'EOS']:
        subprocess.call("ln -s {0:s} {1:s}".format(
            path.abspath(
                path.join(
                    code_path,
                    'hadronic_afterburner_toolkit_code/{}'.format(link_i))),
            path.join(event_folder,
                      "hadronic_afterburner_toolkit/{}".format(link_i))),
                        shell=True)


def create_a_working_folder(workfolder_path):
    try:
        mkdir(workfolder_path)
    except FileExistsError:
        print("The folder {} exists, do you want to delete it?".format(
            workfolder_path))
        user_answer = input()
        if 'y' in user_answer:
            shutil.rmtree(workfolder_path)
            mkdir(workfolder_path)
        else:
            print("bye~\n")
            exit(0)


def main():
    """This is the main funciton"""
    parser = argparse.ArgumentParser(
        description='\U0000269B Welcome to iEBE-MUSIC package',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-w',
                        '--working_folder_name',
                        metavar='',
                        type=str,
                        default='playground',
                        help='working folder path')
    parser.add_argument('-c',
                        '--cluster_name',
                        metavar='',
                        type=str,
                        choices=support_cluster_list,
                        default='local',
                        help='name of the cluster')
    parser.add_argument('-n',
                        '--n_jobs',
                        metavar='',
                        type=int,
                        default=1,
                        help='number of jobs')
    parser.add_argument('-n_hydro',
                        '--n_hydro_per_job',
                        metavar='',
                        type=int,
                        default=1,
                        help='number of hydro events per job to run')
    parser.add_argument('-n_urqmd',
                        '--n_urqmd_per_hydro',
                        metavar='',
                        type=int,
                        default=1,
                        help=('number of oversampled UrQMD events '
                              + 'per hydro to run'))
    parser.add_argument('-n_th',
                        '--n_threads',
                        metavar='',
                        type=int,
                        default=1,
                        help='number of threads used for each job')
    parser.add_argument('-par',
                        '--par_dict',
                        metavar='',
                        type=str,
                        default='parameters_dict_user.py',
                        help='user-defined parameter dictionary file')
    parser.add_argument('-b',
                        '--bayes_file',
                        metavar='',
                        type=str,
                        default='',
                        help='parameters from bayesian analysis')
    parser.add_argument('-id',
                        '--job_process_id',
                        metavar='',
                        type=int,
                        default='0',
                        help='Job process id number')
    parser.add_argument('-seed',
                        '--random_seed',
                        metavar='',
                        type=int,
                        default='-1',
                        help='Random Seed (-1: according to system time)')
    parser.add_argument('--nocopy', action='store_true')
    parser.add_argument("--continueFlag", action="store_true")
    args = parser.parse_args()

    if len(sys.argv) < 2:
        parser.print_help()
        exit(0)

    # print out all the arguments
    print("="*40)
    print("\U0000269B   Input parameters")
    print("="*40)
    for iarg in vars(args):
        print("\U0000269B   {}: {}".format(iarg, getattr(args, iarg)))
    print("="*40)

    try:
        working_folder_name = args.working_folder_name
        cluster_name = args.cluster_name
        n_jobs = args.n_jobs
        n_hydro_per_job = args.n_hydro_per_job
        n_urqmd_per_hydro = args.n_urqmd_per_hydro
        n_threads = args.n_threads
        job_id = args.job_process_id
        seed = args.random_seed
    except:
        parser.print_help()
        exit(0)

    code_package_path = path.abspath(path.dirname(__file__))

    if n_threads < n_urqmd_per_hydro:
        print("\U000026A0  "
              + "Warning: n_threads = {} < n_urqmd_per_hydro = {}!".format(
                  n_threads, n_urqmd_per_hydro))
        print("reset n_threads to {}".format(n_urqmd_per_hydro))
        n_threads = n_urqmd_per_hydro

    par_diretory = path.dirname(path.abspath(args.par_dict))
    sys.path.insert(0, par_diretory)
    parameter_dict = __import__(args.par_dict.split('.py')[0].split("/")[-1])

    if cluster_name == "OSG":
        if seed == -1:
            seed = 0
        seed += job_id
        print("seed = ", seed)
        args.nocopy = True

    initial_condition_type = (parameter_dict.control_dict['initial_state_type'])
    if initial_condition_type not in known_initial_types:
        print("\U0001F6AB  "
              + "Do not recognize the initial condition type: {}".format(
                  initial_condition_type))
        exit(1)

    initial_condition_database = ""
    IPGlasma_time_stamp = "0.4"
    if initial_condition_type == "IPGlasma":
        if parameter_dict.ipglasma_dict['type'] == "self":
            initial_condition_database = "self"
        else:
            initial_condition_database = (
                parameter_dict.ipglasma_dict['database_name_pattern'])
        IPGlasma_time_stamp = str(
            parameter_dict.music_dict['Initial_time_tau_0'])
    elif initial_condition_type == "IPGlasma+KoMPoST":
        if parameter_dict.ipglasma_dict['type'] == "self":
            initial_condition_database = "self"
        else:
            initial_condition_database = (
                parameter_dict.ipglasma_dict['database_name_pattern'])
        IPGlasma_time_stamp = str(
            parameter_dict.kompost_dict['KoMPoSTInputs']['tIn'])
    elif initial_condition_type == "3DMCGlauber_consttau":
        initial_condition_database = (
            parameter_dict.mcglauber_dict['database_name'])
        filelist = glob(
            path.join(initial_condition_database, 'nuclear_thickness_TA_*.dat'))
        nev = max(1, len(filelist))
        print("there are {} events found under the folder {}".format(
            nev, initial_condition_database))
        n_jobs = min(nev, n_jobs)
        n_hydro_per_job = int(ceil(nev/n_jobs))
        print("n_jobs = {}, n_hydro_per_job = {}".format(
            n_jobs, n_hydro_per_job))
    else:
        initial_condition_database = (
            parameter_dict.mcglauber_dict['database_name'])

    working_folder_name = path.abspath(working_folder_name)

    if path.exists(working_folder_name) and args.continueFlag:
        return

    create_a_working_folder(working_folder_name)

    shutil.copy(args.par_dict, working_folder_name)
    code_path = path.join(code_package_path, "codes")
    if not args.nocopy:
        code_path = path.join(working_folder_name, "codes")
        shutil.copytree("{}/codes".format(code_package_path), code_path)

    if args.bayes_file != "":
        args.bayes_file = path.join(path.abspath("."), args.bayes_file)
        subprocess.call("(cd {}/config; ".format(code_package_path)
                        + "python3 parameters_dict_master.py "
                        + "-path {} -par {} -b {} -seed {};)".format(
                            working_folder_name, path.abspath(args.par_dict),
                            args.bayes_file, seed),
                        shell=True)
    else:
        subprocess.call(
            "(cd {}/config; ".format(code_package_path)
            + "python3 parameters_dict_master.py "
            + "-path {} -par {} -seed {};)".format(
                working_folder_name, path.abspath(args.par_dict), seed),
            shell=True)

    if (initial_condition_type not in ("IPGlasma", "IPGlasma+KoMPoST")
            or initial_condition_database != "self"):
        parameter_dict.control_dict['save_ipglasma_results'] = False
    if initial_condition_type != "IPGlasma+KoMPoST":
        parameter_dict.control_dict['save_kompost_results'] = False
    if 'compute_polarization' not in parameter_dict.control_dict.keys():
        parameter_dict.control_dict['compute_polarization'] = False
    if 'calculate_polarization' in parameter_dict.iss_dict.keys():
        if parameter_dict.iss_dict['calculate_polarization'] == 1:
            parameter_dict.control_dict['compute_polarization'] = True
    if 'compute_photon_emission' not in parameter_dict.control_dict.keys():
        parameter_dict.control_dict['compute_photon_emission'] = False


    cent_label = "XXX"
    cent_label_pre = cent_label
    if initial_condition_database == "self":
        print("\U0001F375  Generate initial condition on the fly ... ")
    else:
        initial_condition_database = path.abspath(initial_condition_database)
        print("\U0001F375  "
              + "Pre-generated initial conditions from {} ...".format(
                  initial_condition_database.format(cent_label)))

    toolbar_width = 40
    sys.stdout.write("\U0001F375  Generating {} jobs [{}]".format(
        n_jobs, " "*toolbar_width))
    sys.stdout.flush()
    sys.stdout.write("\b"*(toolbar_width + 1))
    event_id_offset = job_id
    n_hydro_rescaled = n_hydro_per_job
    for iev in range(n_jobs):
        progress_i = (int(float(iev + 1)/n_jobs*toolbar_width)
                      - int(float(iev)/n_jobs*toolbar_width))
        for ii in range(progress_i):
            sys.stdout.write("#")
            sys.stdout.flush()
        if (initial_condition_type in ('IPGlasma', 'IPGlasma+KoMPoST')
                and parameter_dict.ipglasma_dict['type'] == 'minimumbias'):
            precent_local = float(iev)/float(n_jobs)
            for cen_min, cen_max, cen_label, cen_precent in centrality_list:
                if precent_local >= cen_min and precent_local < cen_max:
                    cent_label = cen_label
                    rescale_factor = cen_precent/(cen_max - cen_min)
                    n_hydro_rescaled = (max(
                        1, int(n_hydro_per_job*rescale_factor + 0.1)))
                    if cent_label != cent_label_pre:
                        cent_label_pre = cent_label
                        event_id_offset = 0
                    break
        generate_event_folders(initial_condition_database.format(cent_label),
                               initial_condition_type, code_package_path,
                               code_path, working_folder_name, cluster_name,
                               iev, event_id_offset, n_hydro_rescaled,
                               n_urqmd_per_hydro, n_threads,
                               IPGlasma_time_stamp, parameter_dict)
        event_id_offset += n_hydro_rescaled
    sys.stdout.write("\n")
    sys.stdout.flush()

    # copy script to collect final results
    pwd = path.abspath(".")
    script_path = path.join(code_package_path, "utilities")
    shutil.copy(path.join(script_path, 'collect_events.sh'), pwd)
    shutil.copy(path.join(script_path, 'combine_multiple_hdf5.py'), pwd)
    script_path = (path.join(
        code_package_path,
        "codes/hadronic_afterburner_toolkit_code/ebe_scripts"))
    shutil.copy(path.join(script_path, 'average_event_spvn_h5.py'), pwd)

    walltime = '10:00:00'
    if "walltime" in parameter_dict.control_dict.keys():
        walltime = parameter_dict.control_dict["walltime"]
    if cluster_name == "nersc":
        shutil.copy(
            path.join(code_package_path,
                      'Cluster_supports/NERSC/job_MPI_wrapper.py'),
            working_folder_name)
        n_nodes = max(1, int(n_jobs*n_threads/64))
        generate_nersc_mpi_job_script(working_folder_name, n_nodes, n_threads,
                                      int(n_jobs/n_nodes), walltime)

    if cluster_name == "nerscKNL":
        shutil.copy(
            path.join(code_package_path,
                      'Cluster_supports/NERSC/job_MPI_wrapper.py'),
            working_folder_name)
        n_nodes = max(1, int(n_jobs*n_threads/272))
        generate_nerscKNL_mpi_job_script(working_folder_name, n_nodes,
                                         n_threads, int(n_jobs/n_nodes),
                                         walltime)

    if cluster_name == "wsugrid":
        shutil.copy(
            path.join(code_package_path,
                      'Cluster_supports/WSUgrid/submit_all_jobs.sh'), pwd)


if __name__ == "__main__":
    main()
