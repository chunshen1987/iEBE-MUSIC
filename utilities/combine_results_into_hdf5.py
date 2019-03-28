#!/usr/bin/env python

import sys
import h5py
from os import path, system
from glob import glob

import numpy as np

def print_help():
    print("{0} results_folder".format(sys.argv[0]))

if (len(sys.argv) < 2):
    print_help()
    exit(1)

results_folder = str(sys.argv[1])
results_name = results_folder.split("/")[-1]
if results_name == "":
    results_name = results_folder.split("/")[-2]
results_path = path.abspath(path.join(".", results_folder))
event_list = glob(path.join(results_path, "*.h5"))

for ievent, event_path in enumerate(event_list):
    print("processing {0} ... ".format(event_path))
    event_name = event_path.split("/")[-1]
    hftemp = h5py.File(event_path, "r")
    glist  = list(hftemp.keys())
    hftemp.close()
    for igroup, gtemp in enumerate(glist):
        system('h5copy -i {0} -o {1}.h5 -s {2} -d {2}'.format(
                                            event_path, results_name, gtemp))
