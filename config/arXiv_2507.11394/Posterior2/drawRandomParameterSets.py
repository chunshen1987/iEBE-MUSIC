#!/usr/bin/env python3
"""
    This script translates the posterior chain files in the parameters
    can be read in by the iEBE-MUSIC package.
"""

import pickle
import sys
import numpy as np
import random
import subprocess


with open("posteriorChain.pkl", 'rb') as f:
    data = pickle.load(f)

setName = 'chain'
nParamSets = data[setName].shape[0]

try:
    nSets = int(sys.argv[1])
    ecm = float(sys.argv[2])
    paramFile = str(sys.argv[3])
except:
    print("Usage: drawRandomParameterSets.py <nSets> <ecm> <paramFile>")
    exit(1)

randomList = random.sample(range(0, nParamSets), nSets)
print("Using parameter sets: ", randomList)

for iparam in randomList:
    fileName = f"{paramFile}_{iparam}.txt"
    subprocess.call(
        f"python3 parameterGenerator.py {iparam} 0 {fileName} {ecm}",
        shell=True
    )

