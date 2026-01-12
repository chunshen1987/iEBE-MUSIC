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
centralityCutList = [0., 20., 40., 60]
#centralityCutList = [0, 1, 2, 3, 4, 6, 8, 10, 15, 20, 30, 40, 50, 60,
#                     70, 80, 90, 100]
dNcutList = []  # pre-defined Nch cut if simulation is not minimum bias


def computeJKMeanandErr(dataArr):
    nev, nEta = dataArr.shape
    dataMean = np.mean(dataArr, axis=0)
    dataErr = np.sqrt((nev - 1)/nev*np.sum((dataArr - dataMean)**2, axis=0))
    return dataMean, dataErr


def calculate_rnpT(pTArr, dataTrig, dataAsso,
                   nOrder: int, outputFileName: str) -> None:
    """
        this function calculates the flow pT decorrelation
        r_n(pTa, pTb) = <Qn(pTa)Qn(pTb)^*>/sqrt(<|Qn(pTa)|^2><|Qn(pTb|^2>)

        dataTrig.shape = (nev, npTbins, len(vnVectors))
        vnVectors = [Nch, <pT>, Vn, ET, totalN]
    """
    nev, npTbins, nQn = dataTrig.shape
    QnTrigArr = dataTrig[:, :, nOrder + 1]*dataTrig[:, :, -1]
    QnAssoArr = dataAsso[:, :, nOrder + 1]*dataAsso[:, :, -1]

    # calcualte observables with Jackknife resampling method
    rnpT_array = np.zeros([nev, int(npTbins*(npTbins - 1)/2)])
    for iev in range(nev):
        array_idx = [True]*nev
        array_idx[iev] = False
        array_idx = np.array(array_idx)

        pTidx = 0
        for ipT in range(npTbins - 1):
            for jpT in range(0, ipT + 1):
                if jpT == ipT:
                    rnpT_array[iev, pTidx] = 1.
                else:
                    rnpT_array[iev, pTidx] = (
                        np.real(np.mean(
                            QnTrigArr[array_idx, ipT]
                            *np.conj(QnAssoArr[array_idx, jpT]), axis=0))
                        / np.sqrt(
                            np.mean(np.abs(QnTrigArr[array_idx, ipT])**2,
                                    axis=0)
                            *np.mean(np.abs(QnAssoArr[array_idx, jpT])**2,
                                     axis=0))
                    )
                pTidx += 1

    rnMean, rnErr = computeJKMeanandErr(rnpT_array)

    if path.isfile(outputFileName):
        f = open(outputFileName, 'a')
    else:
        f = open(outputFileName, 'w')
        f.write("# pT^trig (GeV)  pT^asso (GeV)  r_n  r_n_err\n")
    pTidx = 0
    for ipT in range(npTbins - 1):
        pTtrigMid = (pTArr[ipT] + pTArr[ipT + 1]) / 2.
        for jpT in range(0, ipT + 1):
            pTassoMid = (pTArr[jpT] + pTArr[jpT + 1]) / 2.
            f.write("{:.3f}  {:.3f}  {:.5e}  {:.5e}\n".format(
                pTtrigMid, pTassoMid, rnMean[pTidx], rnErr[pTidx]))
            pTidx += 1
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

    cenLabel = "{:02d}-{:02d}".format(int(centralityCutList[icen]),
                                      int(centralityCutList[icen + 1]))
    print("analysis {}%-{}% nev = {}...".format(
        centralityCutList[icen]*centralityRange,
        centralityCutList[icen + 1]*centralityRange, nev))
    print("dNdy: {0:.2f} - {1:.2f}".format(dN_dy_cut_low, dN_dy_cut_high))

    rn_pTArr = data["global"]["rn_pTArr"]
    QnArrTrig = []
    QnArrAsso = []
    for event_name in selected_events_list:
        QnArrTrig.append(data[event_name]['rn_ch_vnpT_trig'])
        QnArrAsso.append(data[event_name]['rn_ch_vnpT_asso'])

    QnArrTrig = np.array(QnArrTrig)
    QnArrAsso = np.array(QnArrAsso)

    calculate_rnpT(rn_pTArr, QnArrTrig, QnArrAsso,
                   2, f"ALICE_r2pT_C{cenLabel}.txt")
    calculate_rnpT(rn_pTArr, QnArrTrig, QnArrAsso,
                   3, f"ALICE_r3pT_C{cenLabel}.txt")
