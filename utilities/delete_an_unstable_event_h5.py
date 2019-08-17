#! /usr/bin/env python3
"""
    This script delete an unstable event in the final collected hdf5 database
"""

from sys import argv
from os import path
import h5py

def main():
    """This is the main function"""
    try:
        data_h5 = path.abspath(argv[1])
        event_id = int(argv[2])
    except IndexError:
        print("Usage: {} database_filename delete_event_id".format(
            str(argv[0])))
        exit(1)
    h5_file = h5py.File(data_h5, "a")
    print("deleting event {} from {} ...".format(event_id, data_h5))
    del h5_file['spvn_results_{}'.format(event_id)]
    h5_file.close()


if __name__ == "__main__":
    main()
