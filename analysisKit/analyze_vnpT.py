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
# centralityCutList = [0, 1, 2, 3, 4, 6, 8, 10, 15, 20, 30, 40, 50, 60,
#                      70, 80, 90, 100]
dNcutList = []  # pre-defined Nch cut if simulation is not minimum bias


def computeJKMeanandErr(dataArr):
    nev, nEta = dataArr.shape
    dataMean = np.mean(dataArr, axis=0)
    dataErr = np.sqrt((nev - 1)/nev*np.sum((dataArr - dataMean)**2, axis=0))
    return dataMean, dataErr


def calculate_vnpT(pTArr, poiSpVn, dataRef, etaRef, nOrder,
                   outputFileName: str) -> None:
    """
        this function compute the v_n(p_T) according to the scalar product
        method
    """
    nev, nQn, nEta = dataRef.shape
    etaArr = np.linspace(-7, 7, nEta)
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
                               dataRef[iev, -1, :]*dataRef[iev, nOrder+1, :])
        Qn2_interp = np.interp(etaRef2Interp, etaArr,
                               dataRef[iev, -1, :]*dataRef[iev, nOrder+1, :])
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

    nev, nQn, npT = poiSpVn.shape
    nQn = nQn - 1
    dNpT = np.real(poiSpVn[:, 0, :]*pTArr.reshape((1, npT))) + 1e-16
    QnpT = dNpT*poiSpVn[:, nOrder, :]

    vnpTNum = np.real(QnpT*np.conj(QnRef1 + QnRef2))
    n2Num = dNpT*(dNRef1 + dNRef2) + 1e-16
    vnpTDen = np.real(QnRef1*np.conj(QnRef2))
    n2Den = dNRef1*dNRef2 + 1e-16

    vnpT_arr = np.zeros([nev, npT])
    for iev in range(nev):
        array_idx = [True]*nev
        array_idx[iev] = False
        array_idx = np.array(array_idx)

        vnpT_arr[iev, :] = (
            np.mean(vnpTNum[array_idx, :], axis=0)
                /np.mean(n2Num[array_idx, :], axis=0)
            /(np.sqrt(np.mean(vnpTDen[array_idx], axis=0)
                    /np.mean(n2Den[array_idx], axis=0)) + 1e-16)
        )

    vnpT_mean, vnpT_err = computeJKMeanandErr(vnpT_arr)
    dNpT_mean = np.mean(dNpT/(pTArr + 1e-16), axis=0)
    dNpT_err = np.std(dNpT/(pTArr + 1e-16), axis=0)/np.sqrt(nev)

    results = np.array([pTArr, dNpT_mean, dNpT_err, vnpT_mean, vnpT_err])
    np.savetxt(outputFileName,
               results.transpose(),
               fmt="%.4e",
               delimiter="  ",
               header="pT (GeV)  dN/d2pT  dN/d2pT_err  vn(pT)  vn(pT)_err")


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

    cenLabel = "{:d}-{:d}".format(
        int(centralityCutList[icen]*centralityRange),
        int(centralityCutList[icen + 1]*centralityRange))
    print("analysis {}%-{}% nev = {}...".format(
        centralityCutList[icen]*centralityRange,
        centralityCutList[icen + 1]*centralityRange, nev))
    print(f"dNdy: {dN_dy_cut_low:.2f} - {dN_dy_cut_high:.2f}")

    chargedpTDiff = []
    QnArrEtapTw = []
    pTArr = []
    for event_name in selected_events_list:
        chargedpTDiff.append(data[event_name]['ch_pTArr'])
        QnArrEtapTw.append(data[event_name]['chVneta_pTw_pT_0p15_2'])
        pTArr = data[event_name]['pTArr']
    chargedpTDiff = np.array(chargedpTDiff)
    QnArrEtapTw = np.array(QnArrEtapTw)

    calculate_vnpT(pTArr, chargedpTDiff, QnArrEtapTw, [0.5, 1.], 2,
                   f"v2pT_ChargedHadron_C{cenLabel}.dat")
