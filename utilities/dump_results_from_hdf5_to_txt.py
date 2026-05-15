#!/usr/bin/env python3
"""This is a script to dump results from a hdf5 file to txt files"""

from os import path, mkdir
import sys
import h5py
import numpy as np


def print_usage():
    """This function prints out help messages"""
    print("Usage: {} ".format(sys.argv[0]) + "results.h5")


def dump_results_from_hdf5(resultsFile):
    """This function combines all the results into hdf5"""
    dataPath = "/".join(path.abspath(resultsFile).split("/")[:-1])
    dataFileName = resultsFile.split("/")[-1].split(".h5")[0]
    mkdir(dataFileName)

    hf = h5py.File("{0}.h5".format(dataFileName), "r")
    groupList = list(hf.keys())
    for group_i in groupList:
        resFolder = path.join(dataFileName, group_i)
        mkdir(resFolder)
        fileList = list(hf[group_i].keys())
        for file_i in fileList:
            data = np.nan_to_num(hf[group_i][file_i])
            attrsList = list(hf[group_i][file_i].attrs.keys())
            headerText = ""
            if "header" in attrsList:
                headerText = str(hf[group_i][file_i].attrs.get("header"))
            np.savetxt(path.join(resFolder, file_i),
                       data,
                       fmt="%.6e",
                       delimiter="  ",
                       header=headerText)
        hf.close()


if __name__ == "__main__":
    try:
        resultsFile = str(sys.argv[1])
        dump_results_from_hdf5(resultsFile)
    except IndexError:
        print_usage()
        exit(0)
