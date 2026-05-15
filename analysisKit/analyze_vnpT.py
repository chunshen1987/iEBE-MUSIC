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
centralityCutList = [0., 20., 40., 60]
# centralityCutList = [0, 1, 2, 3, 4, 6, 8, 10, 15, 20, 30, 40, 50, 60,
#                      70, 80, 90, 100]
dNcutList = []  # pre-defined Nch cut if simulation is not minimum bias


def computeJKMeanandErr(dataArr):
    nev, nEta = dataArr.shape
    dataMean = np.mean(dataArr, axis=0)
    dataErr = np.sqrt((nev - 1)/nev*np.sum((dataArr - dataMean)**2, axis=0))
    return dataMean, dataErr


def calculate_vnSP(pTArr, poiSpVn1, poiSpVn2, etaArr, dataRef, etaRef,
                   nOrder: int, outputFileName: str) -> None:
    """
        this function compute the v_n{SP}(p_T) according to the scalar-product
        method
    """
    nev, nQn, nEta = dataRef.shape
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
        dNRef1.append(np.real(np.sum(Q01_interp)))
        dNRef2.append(np.real(np.sum(Q02_interp)))

    QnRef1 = np.array(QnRef1).reshape((nev, 1))
    QnRef2 = np.array(QnRef2).reshape((nev, 1))
    dNRef1 = np.array(dNRef1).reshape((nev, 1))
    dNRef2 = np.array(dNRef2).reshape((nev, 1))

    nev, nQn, npT = poiSpVn1.shape
    dNpT1 = np.real(poiSpVn1[:, -1, :])
    QnpT1 = dNpT1*poiSpVn1[:, nOrder + 1, :]
    dNpT2 = np.real(poiSpVn2[:, -1, :])
    QnpT2 = dNpT2*poiSpVn2[:, nOrder + 1, :]

    vnpTNum = np.real(QnpT1*np.conj(QnRef1) + QnpT2*np.conj(QnRef2))
    n2Num = dNpT1*dNRef1 + dNpT2*dNRef2 + 1e-16
    vnpTDen = np.real(QnRef1*np.conj(QnRef2))
    n2Den = dNRef1*dNRef2 + 1e-16

    vnpT_arr = np.zeros([nev, npT])
    for iev in range(nev):
        array_idx = [True]*nev
        array_idx[iev] = False
        array_idx = np.array(array_idx)

        vnpT_arr[iev, :] = (np.mean(vnpTNum[array_idx, :], axis=0)
                            /np.mean(n2Num[array_idx, :], axis=0)/(np.sqrt(
                                np.mean(vnpTDen[array_idx], axis=0)
                                /np.mean(n2Den[array_idx], axis=0)) + 1e-16))

    vnpT_mean, vnpT_err = computeJKMeanandErr(vnpT_arr)

    results = np.array([pTArr, vnpT_mean, vnpT_err])
    np.savetxt(outputFileName,
               results.transpose(),
               fmt="%.4e",
               delimiter="  ",
               header="pT (GeV)  vn{SP}(pT)  vn{SP}(pT)_err")


def calculate_vn2PC(pTArr, poiSpVn1, poiSpVn2, nOrder: int,
                    outputFileName: str) -> None:
    """
        this function compute the v_n[2](p_T) according to the 2PC method
    """
    nev, nQn, npT = poiSpVn1.shape
    dNpT1 = np.real(poiSpVn1[:, -1, :])
    QnpT1 = dNpT1*poiSpVn1[:, nOrder + 1, :]
    dNpT2 = np.real(poiSpVn2[:, -1, :])
    QnpT2 = dNpT2*poiSpVn2[:, nOrder + 1, :]

    vnpTNum = np.real(QnpT1*np.conj(QnpT2))
    n2Num = dNpT1*dNpT2 + 1e-16

    vnpT_arr = np.zeros([nev, npT])
    for iev in range(nev):
        array_idx = [True]*nev
        array_idx[iev] = False
        array_idx = np.array(array_idx)

        vnpT_arr[iev, :] = np.nan_to_num(
            np.sqrt(
                np.mean(vnpTNum[array_idx, :], axis=0)/
                (np.mean(n2Num[array_idx, :], axis=0) + 1e-16)))

    vnpT_mean, vnpT_err = computeJKMeanandErr(vnpT_arr)

    results = np.array([pTArr, vnpT_mean, vnpT_err])
    np.savetxt(outputFileName,
               results.transpose(),
               fmt="%.4e",
               delimiter="  ",
               header="pT (GeV)  vn[2](pT)  vn[2](pT)_err")


try:
    database_file = str(sys.argv[1])
except IndexError:
    help_message()

with open(database_file, "rb") as pf:
    data = pickle.load(pf)

dNdyDict = {}
for event_name in data.keys():
    if 'global' not in event_name:
        Nch = np.real(data[event_name]['ALICE_V0A_eta_2p8_5p1_pT_0_4'][0]
                      + data[event_name]['ALICE_V0C_eta_-3p7_-1p7_pT_0_4'][0])
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

    cenLabel = "{:d}-{:d}".format(
        int(centralityCutList[icen]*centralityRange),
        int(centralityCutList[icen + 1]*centralityRange))
    print("analysis {}%-{}% nev = {}...".format(
        centralityCutList[icen]*centralityRange,
        centralityCutList[icen + 1]*centralityRange, nev))
    print(f"dNdy: {dN_dy_cut_low:.2f} - {dN_dy_cut_high:.2f}")

    chargedpTDiff1 = []
    chargedpTDiff2 = []
    QnArrEta = []
    pTArr = data['global']['pTArr']
    etaArr = data['global']['etaArr']
    for event_name in selected_events_list:
        chargedpTDiff1.append(data[event_name]['chVnpT_eta_-0p8_-0p6'])
        chargedpTDiff2.append(data[event_name]['chVnpT_eta_0p6_0p8'])
        QnArrEta.append(data[event_name]['chVneta_pT_0p2_3'])
    chargedpTDiff1 = np.array(chargedpTDiff1)
    chargedpTDiff2 = np.array(chargedpTDiff2)
    QnArrEta = np.array(QnArrEta)

    calculate_vnSP(pTArr, chargedpTDiff1, chargedpTDiff2, etaArr, QnArrEta,
                   [0.6, 0.8], 2, f"v2pT_SP_ChargedHadron_C{cenLabel}.dat")
    calculate_vn2PC(pTArr, chargedpTDiff1, chargedpTDiff2, 2,
                    f"v2pT_2PC_ChargedHadron_C{cenLabel}.dat")
    calculate_vnSP(pTArr, chargedpTDiff1, chargedpTDiff2, etaArr, QnArrEta,
                   [0.6, 0.8], 3, f"v3pT_SP_ChargedHadron_C{cenLabel}.dat")
    calculate_vn2PC(pTArr, chargedpTDiff1, chargedpTDiff2, 3,
                    f"v3pT_2PC_ChargedHadron_C{cenLabel}.dat")
