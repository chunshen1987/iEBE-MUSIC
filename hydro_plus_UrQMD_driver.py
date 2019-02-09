#!/usr/bin/env python

from multiprocessing import Process
from subprocess import call
import sys
import os

def run_UrQMD_event(event_id):
    call("bash ./run_afterburner.sh {0:d}".format(event_id), shell=True)

n_UrQMD_events = int(sys.argv[1])

# first run hydro
#call("bash ./run_hydro.sh", shell=True)

# then run UrQMD events in parallel
for i in range(n_UrQMD_events):
    p = Process(target=run_UrQMD_event, args=(i,))
    p.start()
    p.join()
