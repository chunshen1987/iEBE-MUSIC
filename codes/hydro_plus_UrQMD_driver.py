#!/usr/bin/env python

from multiprocessing import Pool, cpu_count
from subprocess import call
from os import path, mkdir
from glob import glob
import sys
import shutil
from fetch_IPGlasma_event_from_hdf5_database import fecth_an_IPGlasma_event


def print_Usage():
    print("Usage: {} ".format(sys.argv[0]) + "initial_condition_database "
          + "n_hydro_events hydro_event_id n_threads")


def get_initial_condition(database, nev, idx0):
    time_stamp_str = "0.4"
    for iev in range(idx0, idx0 + nev):
        file_name = fecth_an_IPGlasma_event(database, time_stamp_str, iev)
        yield(file_name)
    

def run_hydro_event(final_results_folder, event_id):
    print("Playing MUSIC ... ")
    call("bash ./run_hydro.sh", shell=True)

    # collect hydro results
    hydro_folder_name = "hydro_results_{}".format(event_id)
    shutil.move("MUSIC/hydro_results", path.join(final_results_folder,
                                                 hydro_folder_name))
    return(hydro_folder_name)


def prepare_surface_files_for_UrQMD(final_results_folder, hydro_folder_name,
                                    n_threads):
    surface_file = glob(path.join(final_results_folder, hydro_folder_name,
                                  "surface*.dat"))
    for iev in range(n_threads):
        hydro_surface_folder = "UrQMDev_{0:d}/hydro_event".format(iev)
        mkdir(hydro_surface_folder)
        call("ln -s {0:s} {1:s}".format(
            path.abspath(surface_file[0]),
            path.join(hydro_surface_folder, "surface.dat")), shell=True)
        shutil.copy(path.join(final_results_folder, hydro_folder_name,
                              "music_input"), hydro_surface_folder)

def run_UrQMD_event(event_id):
    call("bash ./run_afterburner.sh {0:d}".format(event_id), shell=True)

def run_UrQMD_shell(n_threads, final_results_folder, event_id):
    print("Running UrQMD ... ")
    with Pool(processes=n_threads) as pool:
        pool.map(run_UrQMD_event, range(n_threads))

    for iev in range(1, n_threads):
        call("./hadronic_afterburner_toolkit/concatenate_binary_files.e "
             + "UrQMDev_0/UrQMD_results/particle_list.gz "
             + "UrQMDev_{}/UrQMD_results/particle_list.gz".format(iev),
             shell=True)
    UrQMD_results_name = "particle_list_{}.gz".format(event_id)
    shutil.move("UrQMDev_0/UrQMD_results/particle_list.gz",
                path.join(final_results_folder, UrQMD_results_name))
    return(path.join(final_results_folder, UrQMD_results_name))


def run_spvn_analysis(pid):
    call("bash ./run_analysis_spvn.sh {0:s}".format(pid), shell=True)

def run_spvn_analysis_shell(UrQMD_file_path, final_results_folder, event_id):
    spvn_folder = "hadronic_afterburner_toolkit/results"
    mkdir(spvn_folder)
    call("ln -s {0:s} {1:s}".format(
         path.abspath(UrQMD_file_path),
         path.join(spvn_folder, "particle_list.dat")), shell=True)
    # finally collect results
    particle_list = [
            '9999', '211', '-211', '321', '-321', '2212', '-2212',
            '3122', '-3122', '3312', '-3312', '3334', '-3334', '333']
    print("Running spvn analysis ... ")
    with Pool(processes=n_threads) as pool:
        pool.map(run_spvn_analysis, particle_list)

    shutil.move(spvn_folder,
                path.join(final_results_folder,
                          "spvn_results_{0:s}".format(event_id)))


def main(initial_condition_database, n_hydro_events, hydro_event_id0,
         n_threads) :
    print("Number of processors: ", cpu_count())
    
    for ifile in get_initial_condition(initial_condition_database,
                                       n_hydro_events, hydro_event_id0):
        print("Run simulations with {} ... ".format(ifile))
        event_id = ifile.split("/")[-1].split("-")[-1].split(".dat")[0]

        final_results_folder="EVENT_RESULTS_{}".format(event_id)
        mkdir(final_results_folder)
        
        shutil.move(ifile, "MUSIC/initial/epsilon-u-Hydro.dat")

        # first run hydro
        hydro_folder_name = run_hydro_event(final_results_folder, event_id)

        prepare_surface_files_for_UrQMD(final_results_folder,
                                        hydro_folder_name, n_threads)

        # then run UrQMD events in parallel
        UrQMD_file_path = run_UrQMD_shell(n_threads, final_results_folder,
                                          event_id)

        # finally collect results
        run_spvn_analysis_shell(UrQMD_file_path, final_results_folder,
                                event_id)


if __name__ == "__main__":
    try:
        initial_condition_database = str(sys.argv[1])
        n_hydro_events             = int(sys.argv[2])
        hydro_event_id0            = int(sys.argv[3])
        n_threads                  = int(sys.argv[4])
    except IndexError:
        print_Usage()
        exit(0)

    main(initial_condition_database, n_hydro_events, hydro_event_id0,
         n_threads)
