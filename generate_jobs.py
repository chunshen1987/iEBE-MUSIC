#!/usr/bin/env python3
"""This script generate all the running jobs."""

import sys
from os import path, mkdir
import shutil
import subprocess
import argparse
from math import ceil
from glob import glob

centrality_list = [(0.00, 0.15, '0-5', 0.05), (0.15, 0.28, '5-10', 0.05),
                   (0.28, 0.40, '10-20', 0.10), (0.40, 0.50, '20-30', 0.10),
                   (0.50, 0.60, '30-40', 0.10), (0.60, 0.69, '40-50', 0.10),
                   (0.69, 0.77, '50-60', 0.10), (0.77, 0.85, '60-70', 0.10),
                   (0.85, 0.92, '70-80', 0.10), (0.92, 0.97, '80-90', 0.10),
                   (0.97, 1.00, '90-100', 0.10)]

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
#PBS -N {0:s}
#PBS -l select=1:ncpus={1:d}:mem={2:.0f}GB:cpu_type=Intel
##PBS -l walltime={3:s}
#PBS -S /bin/bash
#PBS -e test.err
#PBS -o test.log
#PBS -q wsuq

cd {4:s}
""".format(event_id, n_threads, mem, walltime, working_folder))
    elif cluster in ("local", "OSG"):
        script.write("#!/usr/bin/env bash")
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
                             n_hydro, ev0_id, n_urqmd, n_threads, kompost_flag,
                             hydro_flag, urqmd_flag, time_stamp):
    """This function generates full job script"""
    working_folder = folder_name
    event_id = working_folder.split('/')[-1]
    walltime = '100:00:00'

    script = open(path.join(working_folder, "submit_job.pbs"), "w")
    write_script_header(cluster_name, script, n_threads, event_id, walltime,
                        working_folder)
    script.write("\nseed_add=${1:-0}\n" + """
python3 hydro_plus_UrQMD_driver.py {0:s} {1:s} {2:d} {3:d} {4:d} {5:d} {6} {7} {8} $seed_add {9:s} > run.log
""".format(initial_type, database, n_hydro, ev0_id, n_urqmd, n_threads,
           kompost_flag, hydro_flag, urqmd_flag, time_stamp))
    script.close()


def generate_script_kompost(folder_name, nthreads):
    """This function generates script for KoMPoST simulation"""
    working_folder = folder_name

    script = open(path.join(working_folder, "run_kompost.sh"), "w")

    hydro_results_folder = 'kompost_results'
    script.write("""#!/usr/bin/env bash

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

    script.write("""
# KoMPoST EKT evolution
./KoMPoST.exe setup.ini 1> run.log 2> run.err
mv *.txt $results_folder
)
""")
    script.close()


def generate_script_hydro(folder_name, nthreads):
    """This function generates script for hydro simulation"""
    working_folder = folder_name

    script = open(path.join(working_folder, "run_hydro.sh"), "w")

    hydro_results_folder = 'hydro_results'
    script.write("""#!/usr/bin/env bash

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

    script.write("""
# hydro evolution
./MUSIChydro music_input_mode_2 1> run.log 2> run.err
./sweeper.sh $results_folder
)
""")
    script.close()


def generate_script_afterburner(folder_name, GMC_flag=0):
    """This function generates script for hadronic afterburner"""
    working_folder = folder_name

    script = open(path.join(working_folder, "run_afterburner.sh"), "w")
    script.write("""#!/usr/bin/env bash

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
    ./iSS.e > run.log
    """)
    if GMC_flag == 1:
        script.write("""
    # turn on global momentum conservation
    ./correct_momentum_conservation.py OSCAR.DAT
    mv OSCAR_w_GMC.DAT OSCAR.DAT
    """)
    script.write("""
    cd ../osc2u
    ./osc2u.e < ../iSS/OSCAR.DAT >> run.log
    mv fort.14 ../urqmd/OSCAR.input
    cd ../urqmd
    ./runqmd.sh >> run.log
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
    script.write("""#!/usr/bin/env bash

(
    cd hadronic_afterburner_toolkit
    ./hadronic_afterburner_tools.e >> run.log
)
        """)
    script.close()


def generate_event_folders(initial_condition_database, initial_condition_type,
                           code_package_path, working_folder, cluster_name,
                           event_id, event_id_offset, n_hydro_per_job,
                           n_urqmd_per_hydro, n_threads, time_stamp,
                           kompost_flag, hydro_flag, urqmd_flag, GMC_flag):
    """This function creates the event folder structure"""
    event_folder = path.join(working_folder, 'event_%d' % event_id)
    param_folder = path.join(working_folder, 'model_parameters')
    mkdir(event_folder)
    shutil.copy(path.join(working_folder, 'codes/hydro_plus_UrQMD_driver.py'),
                event_folder)
    shutil.copy(
        path.join(code_package_path, 'IPGlasma_database',
                  'fetch_IPGlasma_event_from_hdf5_database.py'), event_folder)
    shutil.copy(
        path.join(code_package_path, '3DMCGlauber_database',
                  'fetch_3DMCGlauber_event_from_hdf5_database.py'),
        event_folder)
    if initial_condition_database == "self":
        shutil.copytree(path.join(working_folder, 'codes/3dMCGlauber'),
                        path.join(event_folder, '3dMCGlauber'),
                        symlinks=True)
        shutil.copyfile(path.join(param_folder, '3dMCGlauber/input'),
                        path.join(event_folder, '3dMCGlauber/input'))
        subprocess.call("ln -s {0:s} {1:s}".format(
            path.abspath(
                path.join(working_folder, 'codes/3dMCGlauber_code/3dMCGlb.e')),
            path.join(event_folder, "3dMCGlauber/3dMCGlb.e")),
                        shell=True)
        subprocess.call("ln -s {0:s} {1:s}".format(
            path.abspath(
                path.join(working_folder, 'codes/3dMCGlauber_code/eps09')),
            path.join(event_folder, "3dMCGlauber/eps09")),
                        shell=True)
        subprocess.call("ln -s {0:s} {1:s}".format(
            path.abspath(
                path.join(working_folder, 'codes/3dMCGlauber_code/tables')),
            path.join(event_folder, "3dMCGlauber/tables")),
                        shell=True)

    generate_full_job_script(cluster_name, event_folder,
                             initial_condition_database, initial_condition_type,
                             n_hydro_per_job, event_id_offset,
                             n_urqmd_per_hydro, n_threads, kompost_flag,
                             hydro_flag, urqmd_flag, time_stamp)

    if initial_condition_type == "IPGlasma+KoMPoST":
        generate_script_kompost(event_folder, n_threads)
        shutil.copytree(path.join(working_folder, 'codes/kompost'),
                        path.join(event_folder, 'kompost'),
                        symlinks=True)
        #shutil.copyfile(path.join(param_folder, 'kompost/setup.ini'),
        #                path.join(event_folder, 'kompost/setup.ini'))
        subprocess.call("ln -s {0:s} {1:s}".format(
            path.abspath(path.join(working_folder, 'codes/kompost_code/EKT')),
            path.join(event_folder, "kompost/EKT")),
                        shell=True)
        subprocess.call("ln -s {0:s} {1:s}".format(
            path.abspath(
                path.join(working_folder, 'codes/kompost_code/KoMPoST.exe')),
            path.join(event_folder, "kompost/KoMPoST.exe")),
                        shell=True)

    generate_script_hydro(event_folder, n_threads)

    shutil.copytree(path.join(working_folder, 'codes/MUSIC'),
                    path.join(event_folder, 'MUSIC'),
                    symlinks=True)
    shutil.copyfile(path.join(param_folder, 'MUSIC/music_input_mode_2'),
                    path.join(event_folder, 'MUSIC/music_input_mode_2'))
    subprocess.call("ln -s {0:s} {1:s}".format(
        path.abspath(path.join(working_folder, 'codes/MUSIC_code/EOS')),
        path.join(event_folder, "MUSIC/EOS")),
                    shell=True)
    subprocess.call("ln -s {0:s} {1:s}".format(
        path.abspath(path.join(working_folder, 'codes/MUSIC_code/MUSIChydro')),
        path.join(event_folder, "MUSIC/MUSIChydro")),
                    shell=True)

    generate_script_afterburner(event_folder, GMC_flag)

    generate_script_analyze_spvn(event_folder)

    for iev in range(n_urqmd_per_hydro):
        sub_event_folder = path.join(working_folder,
                                     'event_{0:d}'.format(event_id),
                                     'UrQMDev_{0:d}'.format(iev))
        mkdir(sub_event_folder)
        shutil.copytree(path.join(working_folder, 'codes/iSS'),
                        path.join(sub_event_folder, 'iSS'))
        shutil.copyfile(path.join(param_folder, 'iSS/iSS_parameters.dat'),
                        path.join(sub_event_folder, 'iSS/iSS_parameters.dat'))
        subprocess.call("ln -s {0:s} {1:s}".format(
            path.abspath(path.join(working_folder,
                                   'codes/iSS_code/iSS_tables')),
            path.join(sub_event_folder, "iSS/iSS_tables")),
                        shell=True)
        subprocess.call("ln -s {0:s} {1:s}".format(
            path.abspath(path.join(working_folder, 'codes/iSS_code/iSS.e')),
            path.join(sub_event_folder, "iSS/iSS.e")),
                        shell=True)
        shutil.copytree(path.join(working_folder, 'codes/osc2u'),
                        path.join(sub_event_folder, 'osc2u'))
        shutil.copytree(path.join(working_folder, 'codes/urqmd'),
                        path.join(sub_event_folder, 'urqmd'))
        subprocess.call("ln -s {0:s} {1:s}".format(
            path.abspath(
                path.join(working_folder, 'codes/urqmd_code/urqmd/urqmd.e')),
            path.join(sub_event_folder, "urqmd/urqmd.e")),
                        shell=True)
    shutil.copytree(
        path.join(working_folder, 'codes/hadronic_afterburner_toolkit'),
        path.join(event_folder, 'hadronic_afterburner_toolkit'))
    shutil.copyfile(
        path.join(param_folder, 'hadronic_afterburner_toolkit/parameters.dat'),
        path.join(event_folder, 'hadronic_afterburner_toolkit/parameters.dat'))
    subprocess.call("ln -s {0:s} {1:s}".format(
        path.abspath(
            path.join(working_folder, 'codes',
                      'hadronic_afterburner_toolkit_code',
                      'hadronic_afterburner_tools.e')),
        path.join(event_folder, "hadronic_afterburner_toolkit",
                  "hadronic_afterburner_tools.e")),
                    shell=True)
    subprocess.call("ln -s {0:s} {1:s}".format(
        path.abspath(
            path.join(working_folder,
                      'codes/hadronic_afterburner_toolkit_code/EOS')),
        path.join(event_folder, "hadronic_afterburner_toolkit/EOS")),
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
                        help=('number of oversampled UrQMD events ' +
                              'per hydro to run'))
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
                        '--OSG_process_id',
                        metavar='',
                        type=int,
                        default='0',
                        help='OSG job process id number')
    args = parser.parse_args()

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
        osg_job_id = args.OSG_process_id
    except:
        parser.print_help()
        exit(0)

    code_package_path = path.abspath(path.dirname(__file__))

    if n_threads < n_urqmd_per_hydro:
        print("\U000026A0  " +
              "Warning: n_threads = {} < n_urqmd_per_hydro = {}!".format(
                  n_threads, n_urqmd_per_hydro))
        print("reset n_threads to {}".format(n_urqmd_per_hydro))
        n_threads = n_urqmd_per_hydro

    par_diretory = path.dirname(path.abspath(args.par_dict))
    sys.path.insert(0, par_diretory)
    parameter_dict = __import__(args.par_dict.split('.py')[0].split("/")[-1])
    initial_condition_type = (parameter_dict.control_dict['initial_state_type'])
    if initial_condition_type not in known_initial_types:
        print("\U0001F6AB  " +
              "Do not recognize the initial condition type: {}".format(
                  initial_condition_type))
        exit(1)

    initial_condition_database = ""
    initial_condition_database_name_pattern = ""
    IPGlasma_time_stamp = "0.4"
    if initial_condition_type == "IPGlasma":
        initial_condition_database = (
            parameter_dict.ipglasma['database_name_pattern'])
        IPGlasma_time_stamp = str(
            parameter_dict.music_dict['Initial_time_tau_0'])
    elif initial_condition_type == "IPGlasma+KoMPoST":
        initial_condition_database = (
            parameter_dict.ipglasma['database_name_pattern'])
        IPGlasma_time_stamp = "0.1"
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
    create_a_working_folder(working_folder_name)

    shutil.copy(args.par_dict, working_folder_name)
    shutil.copytree("{}/codes".format(code_package_path),
                    path.join(working_folder_name, "codes"),
                    ignore=shutil.ignore_patterns("*.git*"))
    if args.bayes_file != "":
        args.bayes_file = path.join(path.abspath("."), args.bayes_file)
        subprocess.call("(cd {}/config; ".format(code_package_path) +
                        "python3 parameters_dict_master.py " +
                        "-path {} -par {} -b {};)".format(
                            working_folder_name, path.abspath(args.par_dict),
                            args.bayes_file),
                        shell=True)
    else:
        subprocess.call(
            "(cd {}/config; ".format(code_package_path) +
            "python3 parameters_dict_master.py " + "-path {} -par {};)".format(
                working_folder_name, path.abspath(args.par_dict)),
            shell=True)

    cent_label = "XXX"
    cent_label_pre = cent_label
    if initial_condition_database == "self":
        print("\U0001F375  Generate initial condition on the fly ... ")
    else:
        initial_condition_database = path.abspath(initial_condition_database)
        print("\U0001F375  " +
              "Pre-generated initial conditions from {} ...".format(
                  initial_condition_database.format(cent_label)))

    toolbar_width = 40
    sys.stdout.write("\U0001F375  Generating {} jobs [{}]".format(
        n_jobs, " "*toolbar_width))
    sys.stdout.flush()
    sys.stdout.write("\b"*(toolbar_width + 1))
    event_id_offset = 0
    n_hydro_rescaled = n_hydro_per_job
    for iev in range(n_jobs):
        progress_i = (int(float(iev + 1)/n_jobs*toolbar_width) -
                      int(float(iev)/n_jobs*toolbar_width))
        for ii in range(progress_i):
            sys.stdout.write("#")
            sys.stdout.flush()
        if (initial_condition_type in ('IPGlasma', 'IPGlasma+KoMPoST') and
                parameter_dict.ipglasma['type'] == 'minimumbias'):
            precent_local = float(iev)/float(n_jobs)
            for cen_min, cen_max, cen_label, cen_precent in centrality_list:
                if precent_local >= cen_min and precent_local < cen_max:
                    cent_label = cen_label
                    rescale_factor = cen_precent/(cen_max - cen_min)
                    n_hydro_rescaled = (max(1,
                                            int(n_hydro_per_job*
                                                rescale_factor)))
                    if cent_label != cent_label_pre:
                        cent_label_pre = cent_label
                        event_id_offset = 0
                    break
        if cluster_name == "OSG":
            event_id_offset = osg_job_id
        GMC_flag = parameter_dict.iss_dict['global_momentum_conservation']
        kompost_flag = False
        if initial_condition_type == "IPGlasma+KoMPoST":
            kompost_flag = parameter_dict.control_dict['save_kompost_results']
        hydro_flag = parameter_dict.control_dict['save_hydro_surfaces']
        urqmd_flag = parameter_dict.control_dict['save_UrQMD_files']
        generate_event_folders(initial_condition_database.format(cent_label),
                               initial_condition_type, code_package_path,
                               working_folder_name, cluster_name, iev,
                               event_id_offset, n_hydro_rescaled,
                               n_urqmd_per_hydro, n_threads,
                               IPGlasma_time_stamp, kompost_flag, hydro_flag,
                               urqmd_flag, GMC_flag)
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
