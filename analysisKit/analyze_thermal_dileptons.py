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


def calculate_dilepton_dNdM(MInvArr, data_dN,
                            outputFileName: str) -> None:
    """
        this function calculate dilepton dN/dM, v2(M)
    """
    nev, nM = data_dN.shape
    dNdM_mean = np.mean(data_dN, axis=0)
    dNdM_err = np.std(data_dN, axis=0)/np.sqrt(nev)
    results = np.array([MInvArr, dNdM_mean, dNdM_err])
    np.savetxt(outputFileName,
               results.transpose(),
               fmt="%.4e",
               delimiter="  ",
               header="M  dN/dM  dN/dM_err")


def calculate_dilepton_dv2dM(MInvArr, photon_dN, etaArr, photon_v2,
                             dataRef, etaRef, nOrder,
                             outputFileName: str) -> None:
    """
        this function compute the v_n(M) according to the scalar product
        method
    """
    nev, nQn, nEta = dataRef.shape
    nev, nM = photon_dN.shape
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

    dNdM = photon_dN
    dilepton_Q2dM = photon_v2*photon_dN

    vnpTNum = np.real(dilepton_Q2dM*np.conj(QnRef1 + QnRef2))
    n2Num = dNdM*(dNRef1 + dNRef2) + 1e-16
    vnpTDen = np.real(QnRef1*np.conj(QnRef2))
    n2Den = dNRef1*dNRef2 + 1e-16

    vndM_arr = np.zeros([nev, nM])
    for iev in range(nev):
        array_idx = [True]*nev
        array_idx[iev] = False
        array_idx = np.array(array_idx)

        vndM_arr[iev, :] = (np.mean(vnpTNum[array_idx, :], axis=0)
                            /np.mean(n2Num[array_idx, :], axis=0)/(np.sqrt(
                                np.mean(vnpTDen[array_idx], axis=0)
                                /np.mean(n2Den[array_idx], axis=0)) + 1e-16))

    vndM_mean, vndM_err = computeJKMeanandErr(vndM_arr)
    dNdM_mean = np.mean(dNdM, axis=0)
    dNdM_err = np.std(dNdM, axis=0)/np.sqrt(nev)

    results = np.array([MInvArr, dNdM_mean, dNdM_err, vndM_mean, vndM_err])
    np.savetxt(outputFileName,
               results.transpose(),
               fmt="%.4e",
               delimiter="  ",
               header="pT (GeV)  dN/dM  dN/dM_err  vn(M)  vn(M)_err")


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
    cenBinMid = (centralityCutList[icen]
                 + centralityCutList[icen + 1])/2.*centralityRange
    print("analysis {}%-{}% nev = {}...".format(
        centralityCutList[icen]*centralityRange,
        centralityCutList[icen + 1]*centralityRange, nev))
    print(f"dNdy: {dN_dy_cut_low:.2f} - {dN_dy_cut_high:.2f}")

    dilepton_dNdM = []
    dilepton_dv2dM = []
    QnArrEta = []
    Ncoll = []
    MInvArr = data[selected_events_list[0]]['dilepton_MInv'][:, 0]
    etaArr = data[selected_events_list[0]]['etaArr']
    nM = len(MInvArr); neta = len(etaArr)
    for event_name in selected_events_list:
        Ncoll.append(data[event_name]['Ncoll'])
        dilepton_dNdM.append(data[event_name]['dilepton_MInv'][:, 1])
        dilepton_dv2dM.append(data[event_name]['dilepton_MInv'][:, 4]
                                 + 1j*data[event_name]['dilepton_MInv'][:, 5])
        QnArrEta.append(data[event_name]['chVneta_pT_0p15_2'])
    Ncoll = np.array(Ncoll)
    dilepton_dNdM = np.array(dilepton_dNdM).reshape(-1, nM)
    dilepton_dv2dM = np.array(dilepton_dv2dM).reshape(-1, nM)
    QnArrEta = np.array(QnArrEta)

    if icen == 0:
        f = open("Ncoll.dat", 'w')
        f.write("# centrality  Ncoll  Ncoll_err\n")
    else:
        f = open("Ncoll.dat", 'a')
    f.write(f"{cenBinMid} {np.mean(Ncoll):.3e} {np.std(Ncoll)/np.sqrt(nev):.3e}\n")

    calculate_dilepton_dNdM(MInvArr, dilepton_dNdM,
                            f"dilepton_dNdM_C{cenLabel}.dat")
    calculate_dilepton_dv2dM(MInvArr, dilepton_dNdM, etaArr, dilepton_dv2dM,
                             QnArrEta, [2.1, 5.1], 2,
                             f"dilepton_dv2dM_C{cenLabel}.dat")
