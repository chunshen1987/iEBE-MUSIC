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
    nev = len(dataArr)
    dataMean = np.mean(dataArr)
    dataErr = np.sqrt((nev - 1)/nev*np.sum((dataArr - dataMean)**2))
    return dataMean, dataErr


def calculate_pTfluct(dataArr1, dataArr2,
                      outputFileHeader: str, cenLabel: str) -> None:
    """
        paper: https://arxiv.org/pdf/2411.09334v2
        this function calculates the moments of pT fluctuation
        mean:            <pT>,
        normalized-std:  sqrt(<delta pT^2>)/<pT>,

        dataArr = [Nch, <pT>, Vn, totalN]
    """
    nev = len(dataArr1[:, 0])
    dN1 = np.real(dataArr1[:, -1])
    dN2 = np.real(dataArr2[:, -1])

    deltaPT_1 = dN1*np.real(dataArr1[:, 1] - np.mean(dataArr1[:, 1]))
    deltaPT_2 = dN2*np.real(dataArr2[:, 1] - np.mean(dataArr2[:, 1]))

    N2_weight = dN1*dN2

    var_dPT = deltaPT_1*deltaPT_2

    # calcualte observables with Jackknife resampling method
    meanPT_array = np.zeros(nev)
    varPT_array = np.zeros(nev)
    for iev in range(nev):
        array_idx = [True]*nev
        array_idx[iev] = False
        array_idx = np.array(array_idx)

        meanPT_array[iev] = (
            (np.mean(dataArr1[array_idx, 1]) + np.mean(dataArr2[array_idx, 1]))
            /2.)
        varPT = np.mean(var_dPT[array_idx]/N2_weight[array_idx])
        varPT_array[iev] = np.sqrt(varPT)/meanPT_array[iev]

    meanPTMean, meanPTErr = computeJKMeanandErr(meanPT_array)
    varPTMean, varPTErr = computeJKMeanandErr(varPT_array)

    pTfluctResults = [
        meanPTMean, meanPTErr, varPTMean, varPTErr,
    ]

    dN_mean = np.real(np.mean(dataArr1[:, 0] + dataArr2[:, 0]))
    dN_err = np.std(dataArr1[:, 0] + dataArr2[:, 0])/np.sqrt(nev)

    outputFileName = outputFileHeader + "pTFluct.dat"
    if path.isfile(outputFileName):
        f = open(outputFileName, 'a')
    else:
        f = open(outputFileName, 'w')
        f.write("# cen  Nch  <pT>  sqrt(var<dpT>)/<pT>\n")
    f.write("{:.3f}  {:.5e}  {:.5e}".format(cenLabel, dN_mean, dN_err))
    for ires in pTfluctResults:
        f.write("  {:.5e}".format(ires))
    f.write("\n")
    f.close()


try:
    database_file = str(sys.argv[1])
except IndexError:
    help_message()

with open(database_file, "rb") as pf:
    data = pickle.load(pf)

dNdyDict = {}
for event_name in data.keys():
    if 'global' not in event_name:
        Nch = np.real(
              data[event_name]['ALICE_V0A_eta_2p8_5p1_pT_0_4'][0]
            + data[event_name]['ALICE_V0C_eta_-3p7_-1p7_pT_0_4'][0]
        )
        dNdyDict[event_name] = Nch
dNdyList = -np.sort(-np.array(list(dNdyDict.values())))
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

    for event_name in dNdyDict.keys():
        if (dNdyDict[event_name] > dN_dy_cut_low
                and dNdyDict[event_name] <= dN_dy_cut_high):
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

    QnArr1 = []
    QnArr2 = []
    for event_name in selected_events_list:
        QnArr1.append(data[event_name]['ALICE_eta_0_0p8_pT_0p15_2'])
        QnArr2.append(data[event_name]['ALICE_eta_-0p8_0_pT_0p15_2'])
    QnArr1 = np.array(QnArr1)
    QnArr2 = np.array(QnArr2)
    calculate_pTfluct(QnArr1, QnArr2, "ALICE", cenLabel)
