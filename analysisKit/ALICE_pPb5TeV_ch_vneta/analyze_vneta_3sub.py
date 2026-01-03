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
centralityCutList = [0, 5, 10, 20, 40]
dNcutList = []  # pre-defined Nch cut if simulation is not minimum bias


def computeJKMeanandErr(dataArr):
    nev, nEta = dataArr.shape
    dataMean = np.mean(dataArr, axis=0)
    dataErr = np.sqrt((nev - 1)/nev*np.sum((dataArr - dataMean)**2, axis=0))
    return dataMean, dataErr


def calculate_vneta_3sub(etaArr, dataArr, dataRef,
                         etaRef1, etaRef2, etaRef3, etaRef4,
                         nOrder: int, outputFileName: str) -> None:
    """
        this function calculates the rapidity distribution of Vn with
        the 3x2PC method
        v_n =  sqrt{<Q_n(eta) (Q_n^Ref_1)*><Q_n(eta) (Q_n^Ref_2)*>
                    /<Q_n^Ref_1 (Q_n^Ref_2)*>}

        dataArr = [Nch, <pT>, Vn, totalN]
    """
    nev, nQn, nEta = dataArr.shape
    nQn = nQn - 3

    etaRef1Interp = np.linspace(etaRef1[0], etaRef1[1], 16)
    etaRef2Interp = np.linspace(etaRef2[0], etaRef2[1], 16)
    etaRef3Interp = np.linspace(etaRef3[0], etaRef3[1], 16)
    etaRef4Interp = np.linspace(etaRef4[0], etaRef4[1], 16)
    QnRef1 = []; dNRef1 = [];
    QnRef2 = []; dNRef2 = [];
    QnRef3 = []; dNRef3 = [];
    QnRef4 = []; dNRef4 = [];
    for iev in range(nev):
        Qn1_interp = np.interp(etaRef1Interp, etaArr,
                               dataRef[iev, -1, :]*dataRef[iev, nOrder + 1, :])
        QnRef1.append(np.sum(Qn1_interp))
        Q01_interp = np.interp(etaRef1Interp, etaArr, dataRef[iev, -1, :])
        dNRef1.append(np.sum(Q01_interp))
        Qn2_interp = np.interp(etaRef2Interp, etaArr,
                               dataRef[iev, -1, :]*dataRef[iev, nOrder + 1, :])
        QnRef2.append(np.sum(Qn2_interp))
        Q02_interp = np.interp(etaRef2Interp, etaArr, dataRef[iev, -1, :])
        dNRef2.append(np.sum(Q02_interp))
        Qn3_interp = np.interp(etaRef3Interp, etaArr,
                               dataRef[iev, -1, :]*dataRef[iev, nOrder + 1, :])
        QnRef3.append(np.sum(Qn3_interp))
        Q03_interp = np.interp(etaRef3Interp, etaArr, dataRef[iev, -1, :])
        dNRef3.append(np.sum(Q03_interp))
        Qn4_interp = np.interp(etaRef4Interp, etaArr,
                               dataRef[iev, -1, :]*dataRef[iev, nOrder + 1, :])
        QnRef4.append(np.sum(Qn4_interp))
        Q04_interp = np.interp(etaRef4Interp, etaArr, dataRef[iev, -1, :])
        dNRef4.append(np.sum(Q04_interp))

    QnRef1 = np.array(QnRef1).reshape((nev, 1))
    QnRef2 = np.array(QnRef2).reshape((nev, 1))
    QnRef3 = np.array(QnRef3).reshape((nev, 1))
    QnRef4 = np.array(QnRef4).reshape((nev, 1))
    dNRef1 = np.array(dNRef1).reshape((nev, 1))
    dNRef2 = np.array(dNRef2).reshape((nev, 1))
    dNRef3 = np.array(dNRef3).reshape((nev, 1))
    dNRef4 = np.array(dNRef4).reshape((nev, 1))

    Qneta = dataArr[:, nOrder + 1, :]*dataArr[:, -1, :]
    dNeta = dataArr[:, -1, :]

    vnNum = np.zeros([nev, nEta])
    n2Num = np.zeros([nev, nEta])
    poiIdx = np.abs(etaArr) < 1.5
    vnNum[:, poiIdx] = np.real(
        Qneta[:, poiIdx]*np.conj(QnRef3)*Qneta[:, poiIdx]*np.conj(QnRef4))
    n2Num[:, poiIdx] = (
        np.real(dNeta[:, poiIdx]*dNRef3 * dNeta[:, poiIdx]*dNRef4) + 1e-16)
    poiIdx = etaArr < -1.5
    vnNum[:, poiIdx] = np.real(
        Qneta[:, poiIdx]*np.conj(QnRef2)*Qneta[:, poiIdx]*np.conj(QnRef4))
    n2Num[:, poiIdx] = (
        np.real(dNeta[:, poiIdx]*dNRef2 * dNeta[:, poiIdx]*dNRef4) + 1e-16)
    poiIdx = etaArr > 1.5
    vnNum[:, poiIdx] = np.real(
        Qneta[:, poiIdx]*np.conj(QnRef1)*Qneta[:, poiIdx]*np.conj(QnRef3))
    n2Num[:, poiIdx] = (
        np.real(dNeta[:, poiIdx]*dNRef1 * dNeta[:, poiIdx]*dNRef3) + 1e-16)

    vnDenMid = np.real(QnRef3*np.conj(QnRef4))
    n2DenMid = np.real(dNRef3*dNRef4) + 1e-16
    vnDenBack = np.real(QnRef2*np.conj(QnRef4))
    n2DenBack = np.real(dNRef2*dNRef4) + 1e-16
    vnDenFoward = np.real(QnRef1*np.conj(QnRef3))
    n2DenFoward = np.real(dNRef1*dNRef3) + 1e-16

    # calcualte observables with Jackknife resampling method
    vnEta_array = np.zeros([nev, nEta])
    for iev in range(nev):
        array_idx = [True]*nev
        array_idx[iev] = False
        array_idx = np.array(array_idx)

        # average weighted by number of particle pairs
        poiIdx = np.abs(etaArr) < 1.5
        vnEta_array[iev, poiIdx] = np.sqrt(
            (np.mean((vnNum[array_idx, :])[:, poiIdx], axis=0)
            / np.mean((n2Num[array_idx, :])[:, poiIdx], axis=0))
            / (np.mean(vnDenMid[array_idx])
               / np.mean(n2DenMid[array_idx]))
        )
        poiIdx = etaArr < -1.5
        vnEta_array[iev, poiIdx] = np.sqrt(
            (np.mean((vnNum[array_idx, :])[:, poiIdx], axis=0)
            / np.mean((n2Num[array_idx, :])[:, poiIdx], axis=0))
            / (np.mean(vnDenBack[array_idx])
               / np.mean(n2DenBack[array_idx]))
        )
        poiIdx = etaArr > 1.5
        vnEta_array[iev, poiIdx] = np.sqrt(
            (np.mean((vnNum[array_idx, :])[:, poiIdx], axis=0)
            / np.mean((n2Num[array_idx, :])[:, poiIdx], axis=0))
            / (np.mean(vnDenFoward[array_idx])
               / np.mean(n2DenFoward[array_idx]))
        )
    vnEta_array = np.nan_to_num(vnEta_array)
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

dNdyDict = {}
for event_name in data.keys():
    if 'global' not in event_name:
        Nch = np.real(
              data[event_name]['ALICE_V0A_eta_-5p1_-2p8_pT_0_4'][0])
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

    etaArr = data['global']['etaArr']
    QnArrEta = []
    for event_name in selected_events_list:
        QnArrEta.append(data[event_name]['chVneta_pT_0_4'])

    QnArrEta = np.array(QnArrEta)

    calculate_vneta_3sub(etaArr, QnArrEta, QnArrEta,
                         [-0.4, 0], [0, 0.4], [-3.1, -2.9], [2.9, 3.1],
                         2, f"ALICE_v2eta_3sub_C{cenLabel}.txt")
