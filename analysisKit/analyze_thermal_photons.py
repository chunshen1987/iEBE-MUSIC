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


def calculate_photon_dNdy(yArr, pTArr, data_ypTdiff,
                          outputFileName: str) -> None:
    """
        this function calculate the photon dN/dy
    """
    dpT = pTArr[1] - pTArr[0]
    nev = data_ypTdiff.shape[0]
    dNdy = np.sum(data_ypTdiff, axis=2)*dpT
    dNdy_mean = np.mean(dNdy, axis=0)
    dNdy_err = np.std(dNdy, axis=0)/np.sqrt(nev)
    results = np.array([yArr, dNdy_mean, dNdy_err])
    np.savetxt(outputFileName,
               results.transpose(),
               fmt="%.4e",
               delimiter="  ",
               header="y  dN/dy  dN/dy_err")


def calculate_photon_pTSpectra(yArr, pTArr, data_ypTdiff,
                               outputFileName: str) -> None:
    """
        this function calculate the photon pT spectra
    """
    nev = data_ypTdiff.shape[0]
    idx = np.abs(yArr) < 0.1
    dy = yArr[1] - yArr[0]
    Yinterval = len(yArr[idx])*dy
    dNd2pT = np.sum(data_ypTdiff[:, idx, :], axis=1)*dy/Yinterval/pTArr/(2*np.pi)
    dNd2pT_mean = np.mean(dNd2pT, axis=0)
    dNd2pT_err = np.std(dNd2pT, axis=0)/np.sqrt(nev)
    results = np.array([pTArr, dNd2pT_mean, dNd2pT_err])
    np.savetxt(outputFileName,
               results.transpose(),
               fmt="%.4e",
               delimiter="  ",
               header="pT(GeV)  dN/(2pi dy pT dpT)  dN/(2pi dy pT dpT)_err")


def calculate_photon_vnpT(yArr, pTArr, photon_dN, photon_v2, rapInterval,
                          dataRef, etaRef, nOrder,
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

    idx = (yArr < rapInterval[1]) & (yArr > rapInterval[0])
    dy = yArr[1] - yArr[0]
    dNd2pT = np.sum(photon_dN[:, idx, :], axis=1)
    QnpT = np.sum(photon_v2[:, idx, :]*photon_dN[:, idx, :], axis=1)
    Yinterval = len(yArr[idx])*dy

    vnpTNum = np.real(QnpT*np.conj(QnRef1 + QnRef2))
    n2Num = dNd2pT*(dNRef1 + dNRef2) + 1e-16
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
    dNd2pT = dNd2pT*dy/Yinterval/(2*np.pi)/pTArr
    dNpT_mean = np.mean(dNd2pT, axis=0)
    dNpT_err = np.std(dNd2pT, axis=0)/np.sqrt(nev)

    results = np.array([pTArr, dNpT_mean, dNpT_err, vnpT_mean, vnpT_err])
    np.savetxt(outputFileName,
               results.transpose(),
               fmt="%.4e",
               delimiter="  ",
               header="pT (GeV)  dN/(2pi dy pT dpT)  dN/(2pi dy pT dpT)_err  vn(pT)  vn(pT)_err")


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

    photon_dN_ypTDiff = []
    photon_v2_ypTDiff = []
    QnArrEta = []
    Ncoll = []
    pTArr = data[selected_events_list[0]]['photon_pTArr']
    yArr = data[selected_events_list[0]]['photon_yArr']
    npT = len(pTArr); ny = len(yArr)
    for event_name in selected_events_list:
        Ncoll.append(data[event_name]['Ncoll'])
        photon_dN_ypTDiff.append(data[event_name]['photon_ypTdiff'][:, 0])
        photon_v2_ypTDiff.append(data[event_name]['photon_ypTdiff'][:, 3]
                                 + 1j*data[event_name]['photon_ypTdiff'][:, 4])
        QnArrEta.append(data[event_name]['chVneta_pT_0p15_2'])
    Ncoll = np.array(Ncoll)
    photon_dN_ypTDiff = np.array(photon_dN_ypTDiff).reshape(-1, ny, npT)
    photon_v2_ypTDiff = np.array(photon_v2_ypTDiff).reshape(-1, ny, npT)
    QnArrEta = np.array(QnArrEta)

    if icen == 0:
        f = open("Ncoll.dat", 'w')
        f.write("# centrality  Ncoll  Ncoll_err\n")
    else:
        f = open("Ncoll.dat", 'a')
    f.write(
        f"{cenBinMid} {np.mean(Ncoll):.3e} {np.std(Ncoll)/np.sqrt(nev):.3e}\n")

    calculate_photon_dNdy(yArr, pTArr, photon_dN_ypTDiff,
                          f"thermalphoton_dNdy_C{cenLabel}.dat")
    calculate_photon_pTSpectra(yArr, pTArr, photon_dN_ypTDiff,
                               f"thermalphoton_pTSpectra_C{cenLabel}.dat")
    calculate_photon_vnpT(yArr, pTArr, photon_dN_ypTDiff, photon_v2_ypTDiff,
                          [-0.5, 0.5], QnArrEta, [0.5, 1.], 2,
                          f"thermalphoton_v2pT_C{cenLabel}.dat")
