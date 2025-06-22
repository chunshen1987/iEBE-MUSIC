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


def calculate_rhon(dataArr1, dataArr2, dataArr3, outputFileHeader: str,
                   cenLabel: str) -> None:
    """
        this function calculates the rho_n correlator between vn and <pT>
        rho_n = (<\hat{delta} vn^2  \hat{delta} pT>)
                /sqrt{<(\hat{delta} vn^2)^2><\hat{delta} pT^2>}
        \hat{delta}O = delta O - delta O*delta N/<(N - <N>)^2>*delta N

        dataArr = [Nch, <pT>, Vn, totalN]
    """
    nev = len(dataArr1[:, 0])
    dN1 = np.real(dataArr1[:, -1])
    dN2 = np.real(dataArr2[:, -1])
    dN3 = np.real(dataArr3[:, -1])

    deltaPT_1 = dN1*np.real(dataArr1[:, 1] - np.mean(dataArr1[:, 1]))

    deltaPT_2 = dN2*np.real(dataArr2[:, 1] - np.mean(dataArr2[:, 1]))
    Q2_2 = dN2*dataArr2[:, 3]
    Q3_2 = dN2*dataArr2[:, 4]
    Q4_2 = dN2*dataArr2[:, 5]
    Q6_2 = dN2*dataArr2[:, 7]
    Q8_2 = dN2*dataArr2[:, 9]

    deltaPT_3 = dN3*np.real(dataArr3[:, 1] - np.mean(dataArr3[:, 1]))
    Q2_3 = dN3*dataArr3[:, 3]
    Q3_3 = dN3*dataArr3[:, 4]
    Q4_3 = dN3*dataArr3[:, 5]
    Q6_3 = dN3*dataArr3[:, 7]
    Q8_3 = dN3*dataArr3[:, 9]

    N3_weight = dN1*dN2*dN3
    cov_Q2dPT = np.real(Q2_2*np.conj(Q2_3))*deltaPT_1
    cov_Q3dPT = np.real(Q3_2*np.conj(Q3_3))*deltaPT_1
    cov_Q4dPT = np.real(Q4_2*np.conj(Q4_3))*deltaPT_1

    N2_weight = dN2*dN3
    var_dPT = deltaPT_2*deltaPT_3
    Q2_N2 = np.real(Q2_2*np.conj(Q2_3))
    Q3_N2 = np.real(Q3_2*np.conj(Q3_3))
    Q4_N2 = np.real(Q4_2*np.conj(Q4_3))

    N4_weight = dN2*(dN2 - 1)*dN3*(dN3 - 1)
    Q2_N4 = (np.real(Q2_2*np.conj(Q2_3))**2 - np.real(Q4_2*np.conj(Q2_3*Q2_3))
             - np.real(Q4_3*np.conj(Q2_2*Q2_2)) + np.real(Q4_2*np.conj(Q4_3)))
    Q3_N4 = (np.real(Q3_2*np.conj(Q3_3))**2 - np.real(Q6_2*np.conj(Q3_3*Q3_3))
             - np.real(Q6_3*np.conj(Q3_2*Q3_2)) + np.real(Q6_2*np.conj(Q6_3)))
    Q4_N4 = (np.real(Q4_2*np.conj(Q4_3))**2 - np.real(Q8_2*np.conj(Q4_3*Q4_3))
             - np.real(Q8_3*np.conj(Q4_2*Q4_2)) + np.real(Q8_2*np.conj(Q8_3)))

    # calcualte observables with Jackknife resampling method
    varPT_array = np.zeros(nev)
    rho2_array = np.zeros(nev)
    covQ2pT_array = np.zeros(nev)
    varQ2_array = np.zeros(nev)
    rho3_array = np.zeros(nev)
    covQ3pT_array = np.zeros(nev)
    varQ3_array = np.zeros(nev)
    rho4_array = np.zeros(nev)
    covQ4pT_array = np.zeros(nev)
    varQ4_array = np.zeros(nev)
    for iev in range(nev):
        array_idx = [True]*nev
        array_idx[iev] = False
        array_idx = np.array(array_idx)

        varPT_array[iev] = (np.mean(var_dPT[array_idx])
                            /np.mean(N2_weight[array_idx]))
        varQ2_array[iev] = (
            np.mean(Q2_N4[array_idx])/np.mean(N4_weight[array_idx]) -
            (np.mean(Q2_N2[array_idx])/np.mean(N2_weight[array_idx]))**2.)
        covQ2pT_array[iev] = np.mean(cov_Q2dPT[array_idx]/N3_weight[array_idx])
        rho2_array[iev] = (covQ2pT_array[iev]
                           /np.sqrt(varQ2_array[iev]*varPT_array[iev]))

        varQ3_array[iev] = (
            np.mean(Q3_N4[array_idx])/np.mean(N4_weight[array_idx]) -
            (np.mean(Q3_N2[array_idx])/np.mean(N2_weight[array_idx]))**2.)
        covQ3pT_array[iev] = np.mean(cov_Q3dPT[array_idx]/N3_weight[array_idx])
        rho3_array[iev] = (covQ3pT_array[iev]
                           /np.sqrt(varQ3_array[iev]*varPT_array[iev]))
        varQ4_array[iev] = (
            np.mean(Q4_N4[array_idx])/np.mean(N4_weight[array_idx]) -
            (np.mean(Q4_N2[array_idx])/np.mean(N2_weight[array_idx]))**2.)
        covQ4pT_array[iev] = np.mean(cov_Q4dPT[array_idx]/N3_weight[array_idx])
        rho4_array[iev] = (covQ4pT_array[iev]
                           /np.sqrt(varQ4_array[iev]*varPT_array[iev]))

    varPTMean, varPTErr = computeJKMeanandErr(varPT_array)
    varQ2Mean, varQ2Err = computeJKMeanandErr(varQ2_array)
    covQ2pTMean, covQ2pTErr = computeJKMeanandErr(covQ2pT_array)
    rho2Mean, rho2Err = computeJKMeanandErr(rho2_array)
    varQ3Mean, varQ3Err = computeJKMeanandErr(varQ3_array)
    covQ3pTMean, covQ3pTErr = computeJKMeanandErr(covQ3pT_array)
    rho3Mean, rho3Err = computeJKMeanandErr(rho3_array)
    varQ4Mean, varQ4Err = computeJKMeanandErr(varQ4_array)
    covQ4pTMean, covQ4pTErr = computeJKMeanandErr(covQ4pT_array)
    rho4Mean, rho4Err = computeJKMeanandErr(rho4_array)
    rho2Results = [
        rho2Mean, rho2Err, covQ2pTMean, covQ2pTErr, varQ2Mean, varQ2Err,
        varPTMean, varPTErr
    ]
    rho3Results = [
        rho3Mean, rho3Err, covQ3pTMean, covQ3pTErr, varQ3Mean, varQ3Err,
        varPTMean, varPTErr
    ]
    rho4Results = [
        rho4Mean, rho4Err, covQ4pTMean, covQ4pTErr, varQ4Mean, varQ4Err,
        varPTMean, varPTErr
    ]

    dN_mean = np.real(np.mean(dataArr2[:, 0] + dataArr3[:, 0]))
    dN_err = np.std(dataArr2[:, 0] + dataArr3[:, 0])/np.sqrt(nev)

    outputFileName = outputFileHeader + "v2pTCorr.dat"
    if path.isfile(outputFileName):
        f = open(outputFileName, 'a')
    else:
        f = open(outputFileName, 'w')
        f.write("# cen  Nch  rho_2  cov(Q2^2 pT)  var(Q2^2)  var(PT)\n")
    f.write("{:.3f}  {:.5e}  {:.5e}".format(cenLabel, dN_mean, dN_err))
    for ires in rho2Results:
        f.write("  {:.5e}".format(ires))
    f.write("\n")
    f.close()

    outputFileName = outputFileHeader + "v3pTCorr.dat"
    if path.isfile(outputFileName):
        f = open(outputFileName, 'a')
    else:
        f = open(outputFileName, 'w')
        f.write("# cen  Nch  rho_3  cov(Q3^2 pT)  var(Q3^2)  var(PT)\n")
    f.write("{:.3f}  {:.5e}  {:.5e}".format(cenLabel, dN_mean, dN_err))
    for ires in rho3Results:
        f.write("  {:.5e}".format(ires))
    f.write("\n")
    f.close()

    outputFileName = outputFileHeader + "v4pTCorr.dat"
    if path.isfile(outputFileName):
        f = open(outputFileName, 'a')
    else:
        f = open(outputFileName, 'w')
        f.write("# cen  Nch  rho_4  cov(Q4^2 pT)  var(Q4^2)  var(PT)\n")
    f.write("{:.3f}  {:.5e}  {:.5e}".format(cenLabel, dN_mean, dN_err))
    for ires in rho4Results:
        f.write("  {:.5e}".format(ires))
    f.write("\n")
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

    QnArr1 = []
    QnArr2 = []
    QnArr3 = []
    for event_name in selected_events_list:
        QnArr1.append(data[event_name]['ALICE_eta_-0p4_0p4'])
        QnArr2.append(data[event_name]['ALICE_eta_-0p8_-0p4'])
        QnArr3.append(data[event_name]['ALICE_eta_0p4_0p8'])
    QnArr1 = np.array(QnArr1)
    QnArr2 = np.array(QnArr2)
    QnArr3 = np.array(QnArr3)
    calculate_rhon(QnArr1, QnArr2, QnArr3, "ALICE", cenLabel)
