#!/usr/bin/env python
"""
     This script fetch some pre-generated IP-Glasma events from a big
     hdf5 database and save them to a smaller hdf5 file
"""

import sys
import h5py
from os import path
from glob import glob

import numpy as np
import random


def print_help():
    print("{0} datafile output_filename nev".format(sys.argv[0]))


def fetch_IPGlasma_events(datafile, outputfile, nev):
    print("fetching {0} events randomly from {1} to {2}.h5 ...".format(
        nev, datafile, outputfile))
    in_data = h5py.File(datafile, "r")
    out_data = h5py.File("{}.h5".format(outputfile), "w")
    event_list = list(in_data.keys())

    if nev > len(event_list):
        print("Error: The original database only has {} events.\n".format(
            len(event_list)))
        print("Error: But {} is requested.\n".format(nev))
        exit(1)

    selected_event_list = random.choices(event_list, k=nev)
    for iev, event_i in enumerate(selected_event_list):
        original_event_id = event_i.split('-')[1]
        in_group = in_data.get(event_i)
        infile_name = "epsilon-u-Hydro-t0.4-{}.dat".format(original_event_id)
        temp_data = in_group.get(infile_name)

        out_group = out_data.create_group("event-{0}".format(iev))
        outfile_name = "epsilon-u-Hydro-t0.4-{}.dat".format(iev)
        dset = out_group.create_dataset("{0}".format(outfile_name),
                                        data = temp_data,
                                        compression="gzip", compression_opts=9)
        data_header = temp_data.attrs["header"].decode('UTF-8').replace('#','')
        dset.attrs.create("header", np.string_(data_header))
        dset.attrs.create("x_size", temp_data.attrs['x_size'])
        dset.attrs.create("y_size", temp_data.attrs['y_size'])
        dset.attrs.create("dx", temp_data.attrs['dx'])
        dset.attrs.create("dy", temp_data.attrs['dy'])
        dset.attrs.create("nx", temp_data.attrs['nx'])
        dset.attrs.create("ny", temp_data.attrs['ny'])

    in_data.close()
    out_data.close()



if __name__ == "__main__":
    try:
        datafile = str(sys.argv[1])
        outputfile = str(sys.argv[2])
        nev = int(sys.argv[3])
        fetch_IPGlasma_events(datafile, outputfile, nev)
    except IndexError:
        print_help()
        exit(1)
