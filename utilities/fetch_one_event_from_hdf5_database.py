#!/usr/bin/env python3

import h5py
import sys
import os
from numpy import *


def help_message():
    print("{0} database_file event_id".format(sys.argv[0]))
    exit(0)


try:
    database_file = str(sys.argv[1])
    event_id = str(sys.argv[2])
except:
    help_message()

h5_data = h5py.File(database_file, "r")
group_name = "spvn_results_{}".format(event_id)
h5_group = h5_data.get(group_name)
print("fetching event {0} from the database {1} ...".format(
    event_id, database_file))

os.mkdir(group_name)

for file_i in list(h5_group.keys()):
    vn_data = h5_group.get(file_i)
    vn_data = nan_to_num(vn_data)

    savetxt("{}/{}".format(group_name, file_i), vn_data,
            fmt="%.6e", delimiter="  ")
