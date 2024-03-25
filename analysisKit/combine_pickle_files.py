#!/usr/bin/env python3

import pickle
import sys
import string
import random


def help_message():
    print("Usage: {0} dataFile1 dataFile2".format(sys.argv[0]))
    print("merge dataFile2 into dataFile1")
    exit(0)


def randomString(stringLength=1):
    """Generate a random string of fixed length """
    letters = string.ascii_lowercase
    return ''.join(random.choice(letters) for i in range(stringLength))


try:
    database_file1 = str(sys.argv[1])
    database_file2 = str(sys.argv[2])
except IndexError:
    help_message()


with open(database_file1, "rb") as pf:
    data1 = pickle.load(pf)

with open(database_file2, "rb") as pf:
    data2 = pickle.load(pf)

keyList = list(data1.keys())
for key_i in list(data2.keys()):
    keyName = key_i
    while keyName in keyList:
        randomlabel = randomString()
        keyName = f"{keyName}{randomlabel}"
    data1[keyName] = data2[key_i]
    keyList.append(keyName)

fileName = "combined_{}".format(database_file1.split("/")[-1])
with open(fileName, "wb") as pf:
    pickle.dump(data1, pf)
