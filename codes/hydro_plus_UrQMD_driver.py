#!/usr/bin/env python

from multiprocessing import Pool, cpu_count
from subprocess import call
from os import path, mkdir
from glob import glob
import sys
import shutil
import h5py
import numpy as np
from fetch_IPGlasma_event_from_hdf5_database import fecth_an_IPGlasma_event
from fetch_3DMCGlauber_event_from_hdf5_database import fecth_an_3DMCGlauber_event


def print_Usage():
    print("Usage: {} ".format(sys.argv[0]) + "initial_condition_database "
          + "initial_condition_type n_hydro_events hydro_event_id n_UrQMD "
          + "n_threads")


def get_initial_condition(database, initial_type, nev, idx0):
    if initial_type == "IPGlasma":
        time_stamp_str = "0.4"
        for iev in range(idx0, idx0 + nev):
            file_name = fecth_an_IPGlasma_event(database, time_stamp_str, iev)
            yield(file_name)
    elif initial_type == "3DMCGlauber":
        for iev in range(idx0, idx0 + nev):
            file_name = "strings_event_{}.dat".format(iev)
            if database == "self":
                call("(cd 3dMCGlauber; ./3dMCGlb.e 1;)", shell=True)
                call("mv 3dMCGlauber/strings_event_0.dat {}".format(file_name),
                     shell=True)
            else:
                file_name = fecth_an_3DMCGlauber_event(database, iev)
            yield(file_name)
    else:
        print("Do not recognize the initial condition type: {}".format(
                                                            initial_type))
        exit(1)
    

def run_hydro_event(final_results_folder, event_id):
    print("Playing MUSIC ... ")
    call("bash ./run_hydro.sh", shell=True)

    # check hydro finishes properly
    ftmp = open("MUSIC/hydro_results/run.log", 'r')
    hydro_status = ftmp.readlines()[-1].split()[3]
    hydro_success = False
    if hydro_status == "Finished.":
        hydro_success = True

    hydro_folder_name = ""
    if hydro_success:
        # collect hydro results
        hydro_folder_name = "hydro_results_{}".format(event_id)
        shutil.move("MUSIC/hydro_results", path.join(final_results_folder,
                                                     hydro_folder_name))
    return(hydro_success, hydro_folder_name)


def prepare_surface_files_for_UrQMD(final_results_folder, hydro_folder_name,
                                    n_UrQMD):
    surface_file = glob(path.join(final_results_folder, hydro_folder_name,
                                  "surface*.dat"))
    for iev in range(n_UrQMD):
        hydro_surface_folder = "UrQMDev_{0:d}/hydro_event".format(iev)
        mkdir(hydro_surface_folder)
        call("ln -s {0:s} {1:s}".format(
            path.abspath(surface_file[0]),
            path.join(hydro_surface_folder, "surface.dat")), shell=True)
        shutil.copy(path.join(final_results_folder, hydro_folder_name,
                              "music_input"), hydro_surface_folder)

def run_UrQMD_event(event_id):
    call("bash ./run_afterburner.sh {0:d}".format(event_id), shell=True)

def run_UrQMD_shell(n_UrQMD, final_results_folder, event_id):
    print("Running UrQMD ... ")
    with Pool(processes=n_UrQMD) as pool:
        pool.map(run_UrQMD_event, range(n_UrQMD))

    for iev in range(1, n_UrQMD):
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

def run_spvn_analysis_shell(UrQMD_file_path, n_threads,
                            final_results_folder, event_id):
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
    
    call("rm {}/particle_list.dat".format(spvn_folder), shell=True)
    shutil.move(spvn_folder,
                path.join(final_results_folder,
                          "spvn_results_{0:s}".format(event_id)))

def zip_results_into_hdf5(final_results_folder, event_id):
    results_name = "spvn_results_{}".format(event_id)
    hydro_info_filepattern = ["eccentricities_evo_eta_*.dat",
                              "momentum_anisotropy_eta_*.dat",
                              "inverse_Reynolds_number_eta_*.dat",
                              "averaged_phase_diagram_trajectory_eta_*.dat"]

    hydrofolder = path.join(final_results_folder,
                            "hydro_results_{}".format(event_id))
    spvnfolder  = path.join(final_results_folder, results_name)

    for ipattern in hydro_info_filepattern:
        hydro_info_list = glob(path.join(hydrofolder, ipattern))
        for ihydrofile in hydro_info_list:
            if path.isfile(ihydrofile):
                shutil.move(ihydrofile, spvnfolder)

    hf = h5py.File("{0}.h5".format(results_name), "w")
    gtemp = hf.create_group("{0}".format(results_name))
    file_list = glob(path.join(spvnfolder, "*"))
    for ifile, file_path in enumerate(file_list):
        file_name = file_path.split("/")[-1]
        dtemp = np.loadtxt(file_path)
        gtemp.create_dataset("{0}".format(file_name), data=dtemp,
                             compression="gzip", compression_opts=9)
    hf.close()
    shutil.move("{}.h5".format(results_name), final_results_folder)


def main(initial_condition_database, initial_condition_type,
         n_hydro_events, hydro_event_id0, n_UrQMD, n_threads):
    print("Number of threads: {}".format(n_threads))
    
    for ifile in get_initial_condition(initial_condition_database,
                                       initial_condition_type,
                                       n_hydro_events, hydro_event_id0):
        print("Run simulations with {} ... ".format(ifile))
        if initial_condition_type == "IPGlasma":
            event_id = ifile.split("/")[-1].split("-")[-1].split(".dat")[0]
            shutil.move(ifile, "MUSIC/initial/epsilon-u-Hydro.dat")
        elif initial_condition_type == "3DMCGlauber":
            event_id = ifile.split("/")[-1].split("_")[-1].split(".dat")[0]
            shutil.move(ifile, "MUSIC/initial/strings.dat")

        final_results_folder="EVENT_RESULTS_{}".format(event_id)
        mkdir(final_results_folder)
        
        # first run hydro
        hydro_success, hydro_folder_name = run_hydro_event(
                                        final_results_folder, event_id)
        
        if not hydro_success:
            # if hydro didn't finish properly, just skip this event
            print("{} did not finsh properly, skipped.".format(
                                        hydro_folder_name))
            continue

        # if hydro finishes properly, we continue to do hadronic transport
        prepare_surface_files_for_UrQMD(final_results_folder,
                                        hydro_folder_name, n_UrQMD)

        # then run UrQMD events in parallel
        UrQMD_file_path = run_UrQMD_shell(n_UrQMD, final_results_folder,
                                          event_id)

        # finally collect results
        run_spvn_analysis_shell(UrQMD_file_path, n_threads,
                                final_results_folder, event_id)

        # zip results into a hdf5 database
        zip_results_into_hdf5(final_results_folder, event_id)


if __name__ == "__main__":
    try:
        initial_condition_type     = str(sys.argv[1])
        initial_condition_database = str(sys.argv[2])
        n_hydro_events             = int(sys.argv[3])
        hydro_event_id0            = int(sys.argv[4])
        n_UrQMD                    = int(sys.argv[5])
        n_threads                  = int(sys.argv[6])
    except IndexError:
        print_Usage()
        exit(0)

    if (initial_condition_type != "IPGlasma"
        and initial_condition_type != "3DMCGlauber"):
        print("Do not recognize the initial condition type: {}".format(
                                                    initial_condition_type))
        exit(1)

    main(initial_condition_database, initial_condition_type,
         n_hydro_events, hydro_event_id0, n_UrQMD, n_threads)
