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
#centralityCutList = Reg_centrality_cut_list
centralityCutList = [0., 10., 20., 40., 60., 80]
# centralityCutList = [0, 1, 2, 3, 4, 6, 8, 10, 15, 20, 30, 40, 50, 60,
#                      70, 80, 90, 100]
dNcutList = []  # pre-defined Nch cut if simulation is not minimum bias


def computeJKMeanandErr(dataArr):
    nev, nEta = dataArr.shape
    dataMean = np.mean(dataArr, axis=0)
    dataErr = np.sqrt((nev - 1)/nev*np.sum((dataArr - dataMean)**2, axis=0))
    return dataMean, dataErr


def calculate_vnpT(pTArr, poiSpVn, etaArr, EPD_QnArr, etaRef,
                   QnArr_C, QnArr_D, nOrder,
                   outputFileName: str) -> None:
    """
        this function compute the v_n(p_T) according to the scalar product
        method
    """
    nev, nQn, nEta = EPD_QnArr.shape
    nev, nQn, npT = poiSpVn.shape
    pTArr = pTArr[:npT]
    etaRef1Interp = np.linspace(etaRef[0], etaRef[1], 16)
    QnRef1 = []
    for iev in range(nev):
        Qn1_interp = np.interp(etaRef1Interp, etaArr,
                               EPD_QnArr[iev, -1, :]*EPD_QnArr[iev, nOrder + 1, :])
        QnRef1.append(np.sum(Qn1_interp))

    QnRef1 = np.array(QnRef1).reshape((nev, 1))

    nev, nQn, npT = poiSpVn.shape
    dNpT = np.real(poiSpVn[:, 0, :]*pTArr.reshape((1, npT))) + 1e-16
    QnpT = dNpT*poiSpVn[:, nOrder + 1, :]

    QnC = (QnArr_C[:, nOrder + 1]*QnArr_C[:, -1]).reshape((nev, 1))
    QnD = (QnArr_D[:, nOrder + 1]*QnArr_D[:, -1]).reshape((nev, 1))

    vnpTNum = np.real(QnpT*np.conj(QnRef1))
    vnBC = np.real(QnRef1*np.conj(QnC))
    vnBD = np.real(QnRef1*np.conj(QnD))
    vnCD = np.real(QnC*np.conj(QnD))

    vnpT_arr = np.zeros([nev, npT])
    for iev in range(nev):
        array_idx = [True]*nev
        array_idx[iev] = False
        array_idx = np.array(array_idx)

        vnpT_arr[iev, :] = (np.mean(vnpTNum[array_idx, :], axis=0)
                            /np.mean(dNpT[array_idx, :], axis=0)
                            /(np.sqrt(np.mean(vnBC[array_idx], axis=0)
                                      *np.mean(vnBD[array_idx], axis=0)
                                      /np.mean(vnCD[array_idx], axis=0)))
        )

    vnpT_mean, vnpT_err = computeJKMeanandErr(vnpT_arr)
    dNpT_mean = np.mean(dNpT/(pTArr + 1e-16), axis=0)
    dNpT_err = np.std(dNpT/(pTArr + 1e-16), axis=0)/np.sqrt(nev)

    results = np.array([pTArr, dNpT_mean, dNpT_err, vnpT_mean, vnpT_err])
    np.savetxt(outputFileName,
               results.transpose(),
               fmt="%.4e",
               delimiter="  ",
               header=("pT (GeV)  dN/(pTdpT)  dN/(pTdpT)_err  "
                       + f"vn(pT)  vn(pT)_err (n = {nOrder})"))

def calculate_vnpT_3sub(pTArr, poiSpVn, etaArr, QnRefArr,
                        etaRef1, etaRef2, etaRef3, nOrder: int, method: int,
                        outputFileName: str) -> None:
    """
        this function compute the v_n(p_T) according to
            method == 0: the scalar-product method
            method == 1: the event-plane metho
        where the denominator is computed with the 3-subevent method.
    """
    nev, nQn, nEta = QnRefArr.shape
    nev, nQn, npT = poiSpVn.shape
    pTArr = pTArr[:npT]
    etaRef1Interp = np.linspace(etaRef1[0], etaRef1[1], 16)
    etaRef2Interp = np.linspace(etaRef2[0], etaRef2[1], 16)
    etaRef3Interp = np.linspace(etaRef3[0], etaRef3[1], 16)
    QnRef1 = []
    QnRef2 = []
    QnRef3 = []
    for iev in range(nev):
        Qn1_interp = np.interp(etaRef1Interp, etaArr,
                               QnRefArr[iev, -1, :]*QnRefArr[iev, nOrder + 1, :])
        QnRef1.append(np.sum(Qn1_interp))
        Qn2_interp = np.interp(etaRef2Interp, etaArr,
                               QnRefArr[iev, -1, :]*QnRefArr[iev, nOrder + 1, :])
        QnRef2.append(np.sum(Qn2_interp))
        Qn3_interp = np.interp(etaRef3Interp, etaArr,
                               QnRefArr[iev, -1, :]*QnRefArr[iev, nOrder + 1, :])
        QnRef3.append(np.sum(Qn3_interp))

    QnRef1 = np.array(QnRef1).reshape((nev, 1))
    QnRef2 = np.array(QnRef2).reshape((nev, 1))
    QnRef3 = np.array(QnRef3).reshape((nev, 1))

    nev, nQn, npT = poiSpVn.shape
    dNpT = np.real(poiSpVn[:, 0, :]*pTArr.reshape((1, npT))) + 1e-16
    QnpT = dNpT*poiSpVn[:, nOrder + 1, :]

    if method == 0:
        vnpTNum = np.real(QnpT*np.conj(QnRef1))
        resPsi1 = QnRef1*np.conj(QnRef2)
        resPsi2 = QnRef1*np.conj(QnRef3)
        resPsi3 = QnRef2*np.conj(QnRef3)
    else:
        vnpTNum = np.real(QnpT*np.conj(QnRef1)/np.abs(QnRef1))
        resPsi1 = QnRef1*np.conj(QnRef2)/(np.abs(QnRef1)*np.abs(QnRef2))
        resPsi2 = QnRef1*np.conj(QnRef3)/(np.abs(QnRef1)*np.abs(QnRef3))
        resPsi3 = QnRef2*np.conj(QnRef3)/(np.abs(QnRef2)*np.abs(QnRef3))

    vnpT_arr = np.zeros([nev, npT])
    for iev in range(nev):
        array_idx = [True]*nev
        array_idx[iev] = False
        array_idx = np.array(array_idx)

        vnpT_arr[iev, :] = (np.mean(vnpTNum[array_idx, :], axis=0)
                            /np.mean(dNpT[array_idx, :], axis=0)
                            /(np.sqrt(np.mean(resPsi1[array_idx], axis=0)
                                      *np.mean(resPsi2[array_idx], axis=0)
                                      /np.mean(resPsi3[array_idx], axis=0)))
        )

    vnpT_mean, vnpT_err = computeJKMeanandErr(vnpT_arr)
    dNpT_mean = np.mean(dNpT/(pTArr + 1e-16), axis=0)
    dNpT_err = np.std(dNpT/(pTArr + 1e-16), axis=0)/np.sqrt(nev)

    results = np.array([pTArr, dNpT_mean, dNpT_err, vnpT_mean, vnpT_err])
    np.savetxt(outputFileName,
               results.transpose(),
               fmt="%.4e",
               delimiter="  ",
               header=("pT (GeV)  dN/(pTdpT)  dN/(pTdpT)_err  "
                       + f"vn(pT)  vn(pT)_err (n = {nOrder})"))

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
    EPDQnArr = []
    chQnEtaArr = []
    QnArrC = []
    QnArrD = []
    pTArr = []
    etaArr = []
    for event_name in selected_events_list:
        chargedpTDiff.append(data[event_name]['chVnpT_eta_-0p5_0p5'])
        EPDQnArr.append(data[event_name]['allParticles_Vneta'])
        chQnEtaArr.append(data[event_name]['chVneta_pT_0p2_2'])
        QnArrC.append(data[event_name]['STAR_eta_-1p5_-0p5_pT_0p2_2'])
        QnArrD.append(data[event_name]['STAR_eta_0p5_1p5_pT_0p2_2'])
        etaArr = data[event_name]['etaArr']
        pTArr = data[event_name]['pTArr']
    chargedpTDiff = np.array(chargedpTDiff)
    EPDQnArr = np.array(EPDQnArr)
    chQnEtaArr = np.array(chQnEtaArr)
    QnArrC = np.array(QnArrC)
    QnArrD = np.array(QnArrD)

    # STAR
    calculate_vnpT_3sub(pTArr, chargedpTDiff, etaArr, EPDQnArr,
                        [-5.1, -2.13], [-1.5, -0.5], [0.5, 1.5], 2, 0,
                        f"v2pT_SP_TPCEPDE_ChargedHadron_C{cenLabel}.dat")
    calculate_vnpT_3sub(pTArr, chargedpTDiff, etaArr, EPDQnArr,
                        [-5.1, -2.13], [-1.5, -0.5], [0.5, 1.5], 3, 0,
                        f"v3pT_SP_TPCEPDE_ChargedHadron_C{cenLabel}.dat")

    # PHENIX
    calculate_vnpT_3sub(pTArr, chargedpTDiff, etaArr, chQnEtaArr,
                        [-3.9, -3.1], [-3., -1], [-0.35, 0.35], 2, 1,
                        f"v2pT_EP_PHENIX_ChargedHadron_C{cenLabel}.dat")
    calculate_vnpT_3sub(pTArr, chargedpTDiff, etaArr, chQnEtaArr,
                        [-3.9, -3.1], [-3., -1], [-0.35, 0.35], 3, 1,
                        f"v3pT_EP_PHENIX_ChargedHadron_C{cenLabel}.dat")
