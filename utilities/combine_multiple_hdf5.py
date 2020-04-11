#!/usr/bin/env python3
"""This script combine multiple hdf5 data files to one"""

import sys
from os import path, system
from glob import glob
import h5py
import string
import random

def randomString(stringLength=1):
    """Generate a random string of fixed length """
    letters = string.ascii_lowercase
    return ''.join(random.choice(letters) for i in range(stringLength))

def print_help():
    """This function outpus help messages"""
    print("{0} results_folder".format(sys.argv[0]))

if len(sys.argv) < 2:
    print_help()
    exit(1)

RESULTS_FOLDER = str(sys.argv[1])
RESULTS_NAME = RESULTS_FOLDER.split("/")[-1]
if RESULTS_NAME == "":
    RESULTS_NAME = RESULTS_FOLDER.split("/")[-2]
RESULTS_PATH = path.abspath(path.join(".", RESULTS_FOLDER))
EVENT_LIST = glob(path.join(RESULTS_PATH, "*.h5"))

exist_group_keys = []

for ievent, event_path in enumerate(EVENT_LIST):
    print("processing {0} ... ".format(event_path))
    event_name = event_path.split("/")[-1]
    hftemp = h5py.File(event_path, "r")
    glist = list(hftemp.keys())
    hftemp.close()
    for igroup, gtemp in enumerate(glist):
        gtemp2 = gtemp
        random_string_len = 1
        tol = 0
        while gtemp2 in exist_group_keys:
            randomlabel = randomString(random_string_len)
            gtemp2 = "{0}{1}".format(gtemp, randomlabel)
            tol += 1
            if tol > 30:
                random_string_len += 1
                tol = 0
        if gtemp2 != gtemp:
            print("Conflict in mergeing {0}, use {1}".format(gtemp, gtemp2))
        exist_group_keys.append(gtemp2)
        system('h5copy -i {0} -o {1}.h5 -s {2} -d {3}'.format(
            event_path, RESULTS_NAME, gtemp, gtemp2))

