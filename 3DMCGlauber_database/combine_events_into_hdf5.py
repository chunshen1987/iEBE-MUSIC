#!/usr/bin/env python
"""
     This script combines pre-generated IP-Glasma events into a hdf5 database.
"""

import sys
import h5py
from os import path
from glob import glob

import numpy as np

def print_help():
    print("{0} results_folder".format(sys.argv[0]))

try:
    results_folder = str(sys.argv[1])
except IndexError:
    print_help()
    exit(1)

results_name = results_folder.split("/")[-1]
if results_name == "":
    results_name = results_folder.split("/")[-2]
results_path = path.abspath(path.join(".", results_folder))

hf = h5py.File("{0}.h5".format(results_name), "w")

# save events summary
event_summary = np.loadtxt(path.join(results_path, "events_summary.dat"))
dset          = hf.create_dataset("events_summary.dat", data = event_summary,
                                  compression="gzip", compression_opts=9)
# save input file
inputfile = np.genfromtxt(path.join(results_path, "input"), dtype='str')
for para_name, para_val in inputfile:
    hf.attrs.create(para_name, np.string_(para_val))

event_list = glob(path.join(results_path, "string*"))
nev = len(event_list)
for ievent, event_path in enumerate(event_list):
    print("processing {0:d}/{1:d} {2} ... ".format(ievent+1, nev,
                                                   results_path))
    file_name = event_path.split("/")[-1]
    dtemp     = np.loadtxt(event_path)
    dset      = hf.create_dataset("{0}".format(file_name), data = dtemp,
                                  compression="gzip", compression_opts=9)
    f = open(event_path)
    header = f.readline().strip('\n')
    dset.attrs.create("header", np.string_(header))
hf.close()
