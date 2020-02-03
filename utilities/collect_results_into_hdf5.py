#!/usr/bin/env python3
"""This is a script to collect one event results into a hdf5 file"""

from os import path, mkdir, remove
from glob import glob
import sys
import shutil
import h5py
import numpy as np


def print_usage():
    """This function prints out help messages"""
    print("Usage: {} ".format(sys.argv[0]) + "results_folder_path ")


def check_an_event_is_good(event_folder):
    """This function checks the given event contains all required files"""
    required_files_list = [
        'particle_9999_vndata_eta_-0.5_0.5.dat',
        'particle_9999_vndata_diff_eta_0.5_2.dat',
        'particle_9999_vndata_eta_-2_2.dat',
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
    event_file_list = glob(path.join(event_folder, "*"))
    for ifile in required_files_list:
        filename = path.join(event_folder, ifile)
        if filename not in event_file_list:
            print("event {} is bad, missing {} ...".format(event_folder,
                                                           filename))
            return False
    return True


def zip_results_into_hdf5(results_folder):
    """This function combines all the results into hdf5"""
    final_results_folder = "/".join(
                            path.abspath(results_folder).split("/")[:-1])
    results_name = results_folder.split("/")[-1]

    spvnfolder = results_folder

    status = check_an_event_is_good(spvnfolder)
    if status:
        print("{} is good, converting results to hdf5".format(spvnfolder))

        hf = h5py.File("{0}.h5".format(results_name), "w")
        gtemp = hf.create_group("{0}".format(results_name))
        file_list = glob(path.join(spvnfolder, "*"))
        for file_path in file_list:
            file_name = file_path.split("/")[-1]
            dtemp = np.loadtxt(file_path)
            h5data = gtemp.create_dataset("{0}".format(file_name), data=dtemp,
                                 compression="gzip", compression_opts=9)
            ftemp = open(file_path, "r")
            header_text = str(ftemp.readline())
            ftemp.close()
            h5data.attrs.create("header", np.string_(header_text))
        hf.close()
        shutil.move("{}.h5".format(results_name), final_results_folder)
    else:
        print("{} is broken, skipped".format(spvnfolder))


if __name__ == "__main__":
    try:
        results_folder = str(sys.argv[1])
        zip_results_into_hdf5(results_folder)
    except IndexError:
        print_usage()
        exit(0)
