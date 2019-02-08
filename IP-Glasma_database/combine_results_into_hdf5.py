#!/usr/bin/env python

import sys
import h5py
from os import path
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
event_list = glob(path.join(results_path, "*"))

hf = h5py.File("{0}.h5".format(results_name), "w")

for ievent, event_path in enumerate(event_list):
    print("processing {0} ... ".format(event_path))
    event_name = event_path.split("/")[-1]
    gtemp = hf.create_group("{0}".format(event_name))
    file_list = glob(path.join(event_path, "*"))
    for ifile, file_path in enumerate(file_list):
        file_name = file_path.split("/")[-1]
        dtemp = np.loadtxt(file_path)
        gtemp.create_dataset("{0}".format(file_name), data=dtemp,
                             compression="gzip", compression_opts=9)

hf.close()
