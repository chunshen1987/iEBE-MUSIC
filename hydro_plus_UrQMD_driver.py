#!/usr/bin/env python

from multiprocessing import Pool, cpu_count
from subprocess import call
import sys

def run_UrQMD_event(event_id):
    call("bash ./run_afterburner.sh {0:d}".format(event_id), shell=True)

if __name__ == "__main__":
    print("Number of processors: ", cpu_count())
    n_UrQMD_events = int(sys.argv[1])
    
    # first run hydro
    print("Playing MUSIC ... ")
    call("bash ./run_hydro.sh", shell=True)
    
    # then run UrQMD events in parallel
    print("Running UrQMD ... ")
    with Pool(processes=n_UrQMD_events) as pool:
        pool.map(run_UrQMD_event, range(n_UrQMD_events))
