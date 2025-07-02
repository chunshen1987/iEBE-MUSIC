import pickle
import sys
from os import path
from glob import glob

import numpy as np

resultFolder = sys.argv[1]

with open(path.join(resultFolder, "HBTresults.pickle"), "rb") as pf:
    HBTdict = pickle.load(pf)

outdata = {}
outdata["HBT"] = HBTdict

fileList = glob(path.join(resultFolder, "resFolder", "*"))
for f in fileList:
    data = np.loadtxt(f)
    outdata[path.basename(f)] = data

with open(path.join(resultFolder, f"{path.basename(resultFolder)}.pickle"),
          "wb") as pf:
    pickle.dump(outdata, pf)
