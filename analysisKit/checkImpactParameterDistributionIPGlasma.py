import h5py
import numpy as np
import matplotlib.pyplot as plt
import sys

h5File = str(sys.argv[1])

data = h5py.File(h5File, 'r')

# the lines of usedParameters<id>.dat are stored as attributes; the impact
# parameter is the line "b = <b> fm" ("# b = <b> fm" since IP-Glasma 2.0)
b = []
for iev in data.keys():
    for line in data[iev].attrs.values():
        words = line.decode("utf-8").lstrip("# ").split()
        if words[:2] == ["b", "="]:
            b.append(float(words[2]))
            break

data.close()

plt.hist(b, 50, density=True)
plt.savefig("impactParameterDistribution.png", dpi=300)
