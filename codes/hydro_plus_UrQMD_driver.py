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
    

def run_UrQMD_event(event_id):
    call("bash ./run_afterburner.sh {0:d}".format(event_id), shell=True)


def main(initial_condition_database, n_hydro_events, hydro_event_id0,
         n_threads) :
    print("Number of processors: ", cpu_count())
    
    final_results_folder="EVENT_RESLUTS"
    mkdir(final_results_folder)

    for ifile in get_initial_condition(initial_condition_database,
                                       n_hydro_events, hydro_event_id0):
        print("Run simulations with {} ... ".format(ifile))
        event_id = ifile.split("/")[-1].split("-")[-1].split(".dat")[0]
        shutil.move(ifile, "MUSIC/initial/epsilon-u-Hydro.dat")

        # first run hydro
        print("Playing MUSIC ... ")
        call("bash ./run_hydro.sh", shell=True)

        # collect hydro results
        hydro_folder_name = "hydro_results_{}".format(event_id)
        shutil.move("MUSIC/hydro_results", path.join(final_results_folder,
                                                     hydro_folder_name))
        
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

        # then run UrQMD events in parallel
        print("Running UrQMD ... ")
        with Pool(processes=n_threads) as pool:
            pool.map(run_UrQMD_event, range(n_threads))


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
