#! /usr/bin/env python3
"""
    This script delete an unstable event in the final collected hdf5 database
"""

from sys import argv
from os import path
import h5py


def check_an_event_is_good(h5_event):
    """This function checks the given event contains all required files"""
    required_files_list = [
        'particle_9999_vndata_eta_-0.5_0.5.dat',
        'particle_9999_vndata_diff_eta_0.5_2.5.dat',
        'particle_9999_vndata_eta_-2.5_2.5.dat',
        'particle_211_vndata_diff_y_-0.5_0.5.dat',
        'particle_321_vndata_diff_y_-0.5_0.5.dat',
        'particle_2212_vndata_diff_y_-0.5_0.5.dat',
        'particle_-211_vndata_diff_y_-0.5_0.5.dat',
        'particle_-321_vndata_diff_y_-0.5_0.5.dat',
        'particle_-2212_vndata_diff_y_-0.5_0.5.dat',
        'particle_3122_vndata_diff_y_-0.5_0.5.dat',
        'particle_3312_vndata_diff_y_-0.5_0.5.dat',
        'particle_3334_vndata_diff_y_-0.5_0.5.dat',
        'particle_-3122_vndata_diff_y_-0.5_0.5.dat',
        'particle_-3312_vndata_diff_y_-0.5_0.5.dat',
        'particle_-3334_vndata_diff_y_-0.5_0.5.dat',
        'particle_333_vndata_diff_y_-0.5_0.5.dat',
    ]
    event_file_list = list(h5_event.keys())
    for ifile in required_files_list:
        if ifile not in event_file_list:
            print("event {} is bad, missing {} ...".format(
                                                        h5_event.name, ifile))
            return False
    return True


def check_events_are_good(h5_filename):
    """This function is a shell to check all the events status in a h5 file"""
    h5_file = h5py.File(h5_filename, "a")
    print("checking {} ...".format(h5_filename))
    event_list = list(h5_file.keys())
    for event_name in event_list:
        test_event = h5_file.get(event_name)
        event_status = check_an_event_is_good(test_event)
        if not event_status:
            print("delete event {} ...".format(event_name))
            del h5_file[event_name]
    print("Finished.")
    h5_file.close()


def delete_an_event_from_hdf5_database(database_name, event_name):
    """This function deletes an event from the hdf5 database"""
    h5_file = h5py.File(database_name, "a")
    print("deleting event {} from {} ...".format(event_name, database_name))
    del h5_file[event_name]
    h5_file.close()


def main():
    """This is the main function"""
    try:
        data_h5 = path.abspath(argv[1])
    except IndexError:
        print("Usage: {} database_filename [delete_event_id]".format(
            str(argv[0])))
        exit(1)

    if len(argv) > 2:
        event_id = int(argv[2])
        event_name = 'spvn_results_{}'.format(event_id)
        delete_an_event_from_hdf5_database(data_h5, event_name)
    else:
        check_events_are_good(data_h5)


if __name__ == "__main__":
    main()
