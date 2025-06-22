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


def calculate_vneta(dataArr, dataRef, etaRef, nOrder,
                    outputFileName: str) -> None:
    """
        this function calculates the rapidity distribution of Vn
        v_n =  <Q_n(eta) (Q_n^Ref)*>/sqrt(Q_n^Ref * (Q_n^Ref)*)

        dataArr = [Nch, <pT>, Vn, totalN]
    """
    nev, nQn, nEta = dataArr.shape
    etaArr = np.linspace(-7, 7, nEta)
    nQn = nQn - 3

    etaRefMin = etaRef[0]
    etaRefMax = etaRef[1]
    etaRef1Interp = np.linspace(etaRefMin, etaRefMax, 16)
    etaRef2Interp = np.linspace(-etaRefMax, -etaRefMin, 16)
    QnRef1 = []
    QnRef2 = []
    dNRef1 = []
    dNRef2 = []
    for iev in range(nev):
        Qn1_interp = np.interp(etaRef1Interp, etaArr,
                               dataRef[iev, -1, :]*dataRef[iev, nOrder + 1, :])
        Qn2_interp = np.interp(etaRef2Interp, etaArr,
                               dataRef[iev, -1, :]*dataRef[iev, nOrder + 1, :])
        Q01_interp = np.interp(etaRef1Interp, etaArr, dataRef[iev, -1, :])
        Q02_interp = np.interp(etaRef2Interp, etaArr, dataRef[iev, -1, :])
        QnRef1.append(np.sum(Qn1_interp))
        QnRef2.append(np.sum(Qn2_interp))
        dNRef1.append(np.sum(Q01_interp))
        dNRef2.append(np.sum(Q02_interp))

    QnRef1 = np.array(QnRef1).reshape((nev, 1))
    QnRef2 = np.array(QnRef2).reshape((nev, 1))
    dNRef1 = np.array(dNRef1).reshape((nev, 1))
    dNRef2 = np.array(dNRef2).reshape((nev, 1))

    Qneta = dataArr[:, nOrder + 1, :]*dataArr[:, -1, :]
    dNeta = dataArr[:, -1, :]

    vnNum = np.zeros([nev, nEta])
    n2Num = np.zeros([nev, nEta])
    posIdx = (etaArr >= 0)
    vnNum[:, posIdx] = np.real(Qneta[:, posIdx]*np.conj(QnRef2))
    n2Num[:, posIdx] = np.real(dNeta[:, posIdx]*dNRef2) + 1e-16
    negIdx = (etaArr < 0)
    vnNum[:, negIdx] = np.real(Qneta[:, negIdx]*np.conj(QnRef1))
    n2Num[:, negIdx] = np.real(dNeta[:, negIdx]*dNRef1) + 1e-16

    vnDen = np.real(QnRef1*np.conj(QnRef2))
    n2Den = np.real(dNRef1*dNRef2) + 1e-16

    # calcualte observables with Jackknife resampling method
    vnEta_array = np.zeros([nev, nEta])
    for iev in range(nev):
        array_idx = [True]*nev
        array_idx[iev] = False
        array_idx = np.array(array_idx)

        # average weighted by number of particle pairs
        #vnEta_array[iev, :] = (
        #    (np.mean(vnNum[array_idx, :], axis=0)
        #     / np.mean(n2Num[array_idx, :], axis=0))
        #    / (np.sqrt(np.mean(vnDen[array_idx])
        #               / np.mean(n2Den[array_idx]) + 1e-16))
        #)

        # STAR definition in https://arxiv.org/pdf/nucl-ex/0409033
        vnEta_array[iev, :] = ((np.mean(vnNum[array_idx, :], axis=0)
                                /np.mean(dNeta[array_idx, :], axis=0))/
                               (np.sqrt(np.mean(vnDen[array_idx]) + 1e-16)))

    vnMean, vnErr = computeJKMeanandErr(vnEta_array)

    if path.isfile(outputFileName):
        f = open(outputFileName, 'a')
    else:
        f = open(outputFileName, 'w')
        f.write("# eta  v_n(eta)  v_n(eta)_err\n")
    for ieta in range(nEta):
        f.write("{:.3f}  {:.5e}  {:.5e}\n".format(etaArr[ieta], vnMean[ieta],
                                                  vnErr[ieta]))
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

    cenLabel = "{:02d}-{:02d}".format(int(centralityCutList[icen]),
                                      int(centralityCutList[icen + 1]))
    print("analysis {}%-{}% nev = {}...".format(
        centralityCutList[icen]*centralityRange,
        centralityCutList[icen + 1]*centralityRange, nev))
    print("dNdy: {0:.2f} - {1:.2f}".format(dN_dy_cut_low, dN_dy_cut_high))

    QnArrEta = []
    QnArrEtapTw = []
    for event_name in selected_events_list:
        QnArrEta.append(data[event_name]['chVneta_pT_0p15_2'])
        QnArrEtapTw.append(data[event_name]['chVneta_pTw_pT_0p15_2'])

    QnArrEta = np.array(QnArrEta)
    QnArrEtapTw = np.array(QnArrEtapTw)

    # use pT weighted vn for reference, https://arxiv.org/pdf/nucl-ex/0409033
    calculate_vneta(QnArrEta, QnArrEtapTw, [2.4, 4.2], 2,
                    f"STAR_v2eta_FTPC_C{cenLabel}.txt")
    calculate_vneta(QnArrEta, QnArrEtapTw, [0.5, 1], 2,
                    f"STAR_v2eta_TPC_C{cenLabel}.txt")
