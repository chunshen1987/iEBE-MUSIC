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
event_list = glob(path.join(results_path, "usedParameters*.dat"))
nev = len(event_list)

print("collect {0} to {1}.h5 ... ".format(results_folder, results_name))
hf = h5py.File("{0}.h5".format(results_name), "w")
for ievent, event_path in enumerate(event_list):
    print("processing {0:d}/{1:d} ... ".format(ievent+1, nev))
    event_id = event_path.split("/")[-1].split("Parameters")[-1].split(".")[0]
    gtemp = hf.create_group("event-{0}".format(event_id))

    file_name = "usedParameters{0}.dat".format(event_id)
    parafile  = open(path.join(results_path, file_name))
    for iline, rawline in enumerate(parafile.readlines()):
        paraline = rawline.strip('\n')
        gtemp.attrs.create("{0}".format(iline), np.string_(paraline))

    file_name = "NcollList{0}.dat".format(event_id)
    filepath = path.join(results_path, file_name)
    if path.isfile(filepath):
        dtemp = np.loadtxt(filepath)
        dtemp = np.nan_to_num(dtemp).reshape(-1, 2)
        dset = gtemp.create_dataset("{0}".format(file_name), data=dtemp,
                                    compression="gzip", compression_opts=9)
    file_name = "NpartList{0}.dat".format(event_id)
    filepath = path.join(results_path, file_name)
    if path.isfile(filepath):
        dtemp = np.loadtxt(path.join(results_path, file_name))
        dtemp = np.nan_to_num(dtemp).reshape(-1, 4)
        dset = gtemp.create_dataset("{0}".format(file_name), data=dtemp,
                                    compression="gzip", compression_opts=9)

    file_name_pattern = "NpartdNdy-t"
    filelist = glob(path.join(results_path, "{0}*-{1}.dat".format(
                                                file_name_pattern, event_id)))
    for ifile, filepath in enumerate(filelist):
        filename = filepath.split("/")[-1]
        dtemp    = np.genfromtxt(filepath, dtype='str')
        data     = np.zeros(len(dtemp))
        for idx in range(len(dtemp)):
            if dtemp[idx] != 'N/A':
                data[idx] = float(dtemp[idx])
            else:
                data[idx] = 0.0
        dset = gtemp.create_dataset("{0}".format(filename), data = data,
                                    compression="gzip", compression_opts=9)

    file_name_pattern = "epsilon-u-Hydro-t"
    filelist = glob(path.join(results_path, "{0}*-{1}.dat".format(
                                                file_name_pattern, event_id)))
    for ifile, filepath in enumerate(filelist):
        filename = filepath.split("/")[-1]
        dtemp    = np.loadtxt(filepath)
        dtemp    = np.nan_to_num(dtemp)
        x_size   = abs(dtemp[0, 1])*2.
        dx       = dtemp[1, 1] - dtemp[0, 1]
        y_size   = abs(dtemp[0, 2])*2.
        dy       = dtemp[1, 1] - dtemp[0, 1]
        nx       = int(np.sqrt(len(dtemp[:, 0])))
        ny       = nx
        data_cut = dtemp[:, 3:]
        dset     = gtemp.create_dataset("{0}".format(filename),
                                        data = data_cut,
                                        compression="gzip", compression_opts=9)
        f = open(filepath)
        header = f.readline().strip('\n')
        dset.attrs.create("header", np.string_(header))
        dset.attrs.create("x_size", x_size)
        dset.attrs.create("y_size", y_size)
        dset.attrs.create("dx", dx)
        dset.attrs.create("dy", dy)
        dset.attrs.create("nx", nx)
        dset.attrs.create("ny", ny)
hf.close()
