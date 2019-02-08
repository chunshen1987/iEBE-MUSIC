#! /usr/bin/env python
"""
     This script fectches individual IP-Glasma event from the hdf5 database.
     It outputs the input file for MUSIC fluid dynamic simulation.
"""

from sys import argv, exit
import numpy as np
import h5py

def print_help():
    print("{0} database_filename event_id".format(argv[0]))

def fecth_an_IPGlasma_event(database_path, time_stamp, event_idx):
    print(("fectching an IP-Glasma event with "
           + "event id: {} at tau = {} fm from {}".format(event_idx,
                                                          time_stamp,
                                                          database_path))
    )
    hf          = h5py.File(database_path, "r")
    event_name  = "event-{0:d}".format(event_idx)
    event_group = hf.get(event_name)
    file_name   = "epsilon-u-Hydro-t{0:s}-{1:d}.dat".format(time_stamp,
                                                            event_idx)
    temp_data   = event_group.get(file_name)
    data_header = temp_data.attrs["header"].decode('UTF-8').replace('#','')
    np.savetxt(file_name, temp_data, fmt='%.6e', header=data_header)

if __name__ == "__main__":
    try:
        database_filename = str(argv[1])
        event_id          = int(argv[2])
    except IndexError:
        print_help()
        exit(1)
    
    time_stamp_str = "0.4"
    fecth_an_IPGlasma_event(database_filename, time_stamp_str, event_id)
