#!/usr/bin/env python3

import sys
from os import path
import pickle
import numpy as np


def help_message():
    print("Usage: {0} database_file".format(sys.argv[0]))
    exit(0)


centralityRange = 1.
Reg_centrality_cut_list = [
    0., 5., 10., 20., 30., 40., 50., 60., 70., 80., 90., 100.
]
centralityCutList = Reg_centrality_cut_list
#centralityCutList = [0, 1, 2, 3, 4, 6, 8, 10, 15, 20, 30, 40, 50, 60,
#                     70, 80, 90, 100]
dNcutList = []  # pre-defined Nch cut if simulation is not minimum bias


def computeJKMeanandErr(dataArr):
    nev, nEta = dataArr.shape
    dataMean = np.mean(dataArr, axis=0)
    dataErr = np.sqrt((nev - 1)/nev*np.sum((dataArr - dataMean)**2, axis=0))
    return dataMean, dataErr


def calculate_dNdeta(etaArr, dNdetaArr, outputFileName: str) -> None:
    """
        this function calculates the rapidity distribution of dN/deta

        etaArr[nev, neta]
    """
    etaMean = np.sum(etaArr*dNdetaArr, axis=0)/np.sum(dNdetaArr, axis=0)
    dNMean = np.mean(dNdetaArr, axis=0)
    dNerr = np.std(dNdetaArr, axis=0)/np.sqrt(len(dNdetaArr[:, 0]))

    if path.isfile(outputFileName):
        f = open(outputFileName, 'a')
    else:
        f = open(outputFileName, 'w')
        f.write("# eta  dN/deta  dN/deta_err\n")
    for ieta in range(len(etaMean)):
        f.write("{:.3f}  {:.5e}  {:.5e}\n".format(etaMean[ieta], dNMean[ieta],
                                                  dNerr[ieta]))
    f.close()


def calculate_dETdeta(etaArr, dETdetaArr, outputFileName: str) -> None:
    """
        this function calculates the rapidity distribution of dET/deta

        etaArr[nev, neta]
    """
    etaMean = np.sum(etaArr*dETdetaArr, axis=0)/np.sum(dETdetaArr, axis=0)
    dETMean = np.mean(dETdetaArr, axis=0)
    dETerr = np.std(dETdetaArr, axis=0)/np.sqrt(len(dETdetaArr[:, 0]))

    if path.isfile(outputFileName):
        f = open(outputFileName, 'a')
    else:
        f = open(outputFileName, 'w')
        f.write("# eta  dET/deta  dET/deta_err\n")
    for ieta in range(len(etaMean)):
        f.write("{:.3f}  {:.5e}  {:.5e}\n".format(etaMean[ieta], dETMean[ieta],
                                                  dETerr[ieta]))
    f.close()


try:
    database_file = str(sys.argv[1])
except IndexError:
    help_message()

with open(database_file, "rb") as pf:
    data = pickle.load(pf)

dNdyList = []
for event_name in data.keys():
    Nch = data[event_name]['Nch']
    dNdyList.append(Nch)
dNdyList = -np.sort(-np.array(dNdyList))
print(f"Number of good events: {len(dNdyList)}")

for icen in range(len(centralityCutList) - 1):
    if centralityCutList[icen + 1] < centralityCutList[icen]:
        continue
    selected_events_list = []

    dN_dy_cut_high = dNdyList[int(len(dNdyList)*centralityCutList[icen]/100.)]
    dN_dy_cut_low = dNdyList[min(
        len(dNdyList) - 1, int(len(dNdyList)*centralityCutList[icen + 1]/100.))]

    if len(dNcutList) == len(centralityCutList):
        dN_dy_cut_high = dNcutList[icen]
        dN_dy_cut_low = dNcutList[icen + 1]

    for event_name in data.keys():
        Nch = data[event_name]['Nch']
        if Nch > dN_dy_cut_low and Nch <= dN_dy_cut_high:
            selected_events_list.append(event_name)

    nev = len(selected_events_list)
    if nev <= 0:
        continue

    cenLabel = "{:02d}-{:02d}".format(int(centralityCutList[icen]),
                                      int(centralityCutList[icen + 1]))
    print("analysis {}%-{}% nev = {}...".format(
        centralityCutList[icen]*centralityRange,
        centralityCutList[icen + 1]*centralityRange, nev))
    print("dNdy: {0:.2f} - {1:.2f}".format(dN_dy_cut_low, dN_dy_cut_high))

    dNdetaArr = []
    dETdetaArr = []
    etaArr = []
    for event_name in selected_events_list:
        etaArr.append(data[event_name]['etaArr'])
        dNdetaArr.append(data[event_name]['dNch/deta'])
        dETdetaArr.append(data[event_name]['dET/deta'])

    etaArr = np.array(etaArr)
    dNdetaArr = np.array(dNdetaArr)
    dETdetaArr = np.array(dETdetaArr)

    calculate_dNdeta(etaArr, dNdetaArr, f"STAR_dNdeta_C{cenLabel}.txt")
    calculate_dETdeta(etaArr, dETdetaArr, f"STAR_dETdeta_C{cenLabel}.txt")
