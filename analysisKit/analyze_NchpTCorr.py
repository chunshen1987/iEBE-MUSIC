#!/usr/bin/env python3

import sys
from os import path
import pickle
import numpy as np


def help_message():
    print("Usage: {0} database_file".format(sys.argv[0]))
    exit(0)


centralityRange = 1.
centralityCutList = [0, 5]
dNcutList = []  # pre-defined Nch cut if simulation is not minimum bias


def computeJKMeanandErr(dataArr):
    nev = len(dataArr)
    dataMean = np.mean(dataArr)
    dataErr = np.sqrt((nev - 1)/nev*np.sum((dataArr - dataMean)**2))
    return dataMean, dataErr


def calculate_NchpTCorr(NchArr, meanpTArr) -> None:
    """
        this function calculates the histrogram of mean pT vs Nch
    """
    normalizedNchArr = NchArr/np.mean(NchArr)
    normalizedmeanpTArr = meanpTArr/np.mean(meanpTArr)

    NchBins = np.linspace(0.8, 1.3, 11)
    nBins = len(NchBins) - 1
    pTHist = np.zeros(nBins)
    NchHist = np.zeros(nBins)
    pTErr = np.zeros(nBins)
    NchErr = np.zeros(nBins)
    nevArr = []
    for ibin in range(nBins):
        idx = ((normalizedNchArr >= NchBins[ibin]) &
               (normalizedNchArr < NchBins[ibin + 1]))
        nev = len(normalizedNchArr[idx])
        nevArr.append(nev)
        if nev > 1:
            pTHist[ibin] = np.mean(normalizedmeanpTArr[idx])
            NchHist[ibin] = np.mean(normalizedNchArr[idx])
            pTErr[ibin] = np.std(normalizedmeanpTArr[idx])/np.sqrt(nev)
            NchErr[ibin] = np.std(normalizedNchArr[idx])/np.sqrt(nev)

    outputFileName = "Nch_pT.dat"
    f = open(outputFileName, 'w')
    f.write("# Nch  Nch_err  <pT>  <pT>_err  nev\n")
    for ibin in range(nBins):
        f.write(f"{NchHist[ibin]:.4e}  {NchErr[ibin]:.4e}  "
                + f"{pTHist[ibin]:.4e}  {pTErr[ibin]:.4e}  {nevArr[ibin]}\n")
    f.close()


try:
    database_file = str(sys.argv[1])
except IndexError:
    help_message()

with open(database_file, "rb") as pf:
    data = pickle.load(pf)

dNdyList = []
for event_name in data.keys():
    dNdyList.append(data[event_name]['Nch'])
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
        if (data[event_name]['Nch'] > dN_dy_cut_low
                and data[event_name]['Nch'] <= dN_dy_cut_high):
            selected_events_list.append(event_name)

    nev = len(selected_events_list)
    if nev <= 0:
        continue

    cenLabel = (centralityCutList[icen]
                + centralityCutList[icen + 1])/2.*centralityRange
    print("analysis {}%-{}% nev = {}...".format(
        centralityCutList[icen]*centralityRange,
        centralityCutList[icen + 1]*centralityRange, nev))
    print("dNdy: {0:.2f} - {1:.2f}".format(dN_dy_cut_low, dN_dy_cut_high))

    NchArr = []
    meanpTArr = []
    for event_name in selected_events_list:
        NchArr.append(data[event_name]['Nch'])
        meanpTArr.append(data[event_name]['mean_pT_ch'])
    NchArr = np.array(NchArr)
    meanpTArr = np.array(meanpTArr)
    calculate_NchpTCorr(NchArr, meanpTArr)
