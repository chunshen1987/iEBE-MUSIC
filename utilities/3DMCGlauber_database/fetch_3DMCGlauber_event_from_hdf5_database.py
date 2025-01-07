#! /usr/bin/env python3
"""
     This script fetches individual 3D MC-Glauber event from the hdf5 database.
     It outputs the input file for MUSIC fluid dynamic simulation.
"""

from sys import argv, exit
import numpy as np
import h5py


def print_help():
    print("{0} database_filename event_id".format(argv[0]))


def fetch_an_3DMCGlauber_event(database_path, event_idx):
    print(("fetching an 3DMCGlauber event with "
           + "event id: {} from {}".format(event_idx, database_path)))
    hf = h5py.File(database_path, "r")
    file_name = "strings_event_{0:d}.dat".format(event_idx)
    temp_data = hf.get(file_name)
    data_header = temp_data.attrs["header"].decode('UTF-8').replace('#', '')
    temp_data = np.array(temp_data).reshape(-1, 21)
    np.savetxt(file_name, temp_data, fmt='%.6e', header=data_header)
    return (file_name)


if __name__ == "__main__":
    try:
        database_filename = str(argv[1])
        event_id = int(argv[2])
    except IndexError:
        print_help()
        exit(1)

    fetch_an_3DMCGlauber_event(database_filename, event_id)
