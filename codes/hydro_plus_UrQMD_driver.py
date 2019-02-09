#!/usr/bin/env python

from multiprocessing import Pool, cpu_count
from subprocess import call
from os import path, mkdir
from glob import glob
import sys
import shutil

def print_Usage():
    print("Usage: {} n_UrQMD_ev ".format(sys.argv[0]))


def get_initial_condition():
    filelist = glob("Initial/*")
    for ifile in filelist:
        yield(ifile)
    

def run_UrQMD_event(event_id):
    call("bash ./run_afterburner.sh {0:d}".format(event_id), shell=True)


def main(n_UrQMD_events):
    print("Number of processors: ", cpu_count())
    
    mkdir("event_results")
    for ifile in get_initial_condition():
        print("Run simulations with {} ... ".format(ifile))
        event_id = ifile.split("/")[-1].split("-")[-1].split(".dat")[0]
        shutil.move(ifile, "MUSIC/initial/epsilon-u-Hydro.dat")

        # first run hydro
        print("Playing MUSIC ... ")
        call("bash ./run_hydro.sh", shell=True)

        # collect hydro results
        shutil.move("MUSIC/hydro_results",
                    "event_results/hydro_results_{}".format(event_id))
        
        # then run UrQMD events in parallel
        print("Running UrQMD ... ")
        with Pool(processes=n_UrQMD_events) as pool:
            pool.map(run_UrQMD_event, range(n_UrQMD_events))


if __name__ == "__main__":
    try:
        n_UrQMD_events = int(sys.argv[1])
    except IndexError:
        print_Usage()
        exit(0)

    main(n_UrQMD_events)
