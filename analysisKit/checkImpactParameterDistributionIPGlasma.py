import h5py
import numpy as np
import matplotlib.pyplot as plt
import sys

h5File = str(sys.argv[1])

data = h5py.File(h5File, 'r')

b = []
for iev in data.keys():
    b.append(float(data[iev].attrs['39'].decode("utf-8").split(" ")[2]))

data.close()

plt.hist(b, 50, density=True)
plt.savefig("impactParameterDistribution.png", dpi=300)
