#!/usr/bin/env python3

import pickle
import sys
import string
import random


def help_message():
    print("Usage: {0} dataFile1 dataFile2 ...".format(sys.argv[0]))
    print("merge dataFile2 ... into dataFile1")
    exit(0)


if len(sys.argv) < 3:
    help_message()

database_file1 = str(sys.argv[1])
database_fileList = [str(i) for i in sys.argv[2:]]

with open(database_file1, "rb") as pf:
    data1 = pickle.load(pf)
print(f"File {database_file1} has {len(list(data1.keys()))} entries ...")

letters = string.ascii_lowercase
for i, database_file2 in enumerate(database_fileList):
    keyList = set(data1.keys())
    with open(database_file2, "rb") as pf:
        data2 = pickle.load(pf)
    print(f"File {database_file2} has {len(list(data2.keys()))} entries ...")
    key2 = set(data2.keys())
    if len(keyList.intersection(key2)) > 0:
        addLetter = letters[i]
        for key_i in list(data2.keys()):
            newKey = f"{key_i}_{addLetter}"
            data2[newKey] = data2.pop(key_i)
    data1 = data1 | data2

with open("combined.pkl", "wb") as pf:
    pickle.dump(data1, pf)
