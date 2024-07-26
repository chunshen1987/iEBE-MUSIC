#!/usr/bin/env python3

import sys
from os import path
import pickle
import numpy as np


def help_message():
    print("Usage: {0} database_file".format(sys.argv[0]))
    exit(0)


centralityRange = 1.
Reg_centrality_cut_list = [0., 5., 10., 20., 30., 40., 50.,
                           60., 70., 80., 90., 100.]
centralityCutList = Reg_centrality_cut_list
#centralityCutList = [0, 1, 2, 3, 4, 6, 8, 10, 15, 20, 30, 40, 50, 60,
#                     70, 80, 90, 100]
dNcutList = []    # pre-defined Nch cut if simulation is not minimum bias


def computeJKMeanandErr(dataArr):
    nev = len(dataArr)
    dataMean = np.mean(dataArr)
    dataErr = np.sqrt((nev - 1)/nev*np.sum((dataArr - dataMean)**2))
    return dataMean, dataErr


def calculate_pTfluct(dataArr1, dataArr2, dataArr3,
                      outputFileHeader: str, cenLabel: str) -> None:
    """
        this function calculates the moments of pT fluctuation
        mean:      <pT>,
        variance:  <delta pT^2>/<pT>^2,
        sknewness: <delta pT^3>/std(delta pT)^3
        kurtosis:  <delta pT^4>/std(delta pT)^4

        dataArr = [Nch, <pT>, Vn, totalN]
    """
    nev = len(dataArr1[:, 0])
    dN1 = np.real(dataArr1[:, -1])
    dN2 = np.real(dataArr2[:, -1])
    dN3 = np.real(dataArr3[:, -1])

    deltaPT_1 = dN1*np.real(dataArr1[:, 1] - np.mean(dataArr1[:, 1]))
    deltaPT_2 = dN2*np.real(dataArr2[:, 1] - np.mean(dataArr2[:, 1]))
    deltaPT_3 = dN3*np.real(dataArr3[:, 1] - np.mean(dataArr3[:, 1]))

    N2_weight = dN1*dN2
    N3_weight = dN1*dN2*dN3
    N4_weight = dN1*dN1*dN2*dN2

    var_dPT = deltaPT_1*deltaPT_2
    sknew_dPT = deltaPT_1*deltaPT_2*deltaPT_3
    kurtosis_dPT = deltaPT_1*deltaPT_2*deltaPT_1*deltaPT_2

    # calcualte observables with Jackknife resampling method
    meanPT_array = np.zeros(nev)
    varPT_array = np.zeros(nev)
    sknewPT_array = np.zeros(nev)
    kurtosisPT_array = np.zeros(nev)
    for iev in range(nev):
        array_idx      = [True]*nev
        array_idx[iev] = False
        array_idx      = np.array(array_idx)

        meanPT_array[iev] = ((np.mean(dataArr1[array_idx, 1])
                              + np.mean(dataArr2[array_idx, 1]))/2.)
        varPT = np.mean(var_dPT[array_idx])/np.mean(N2_weight[array_idx])
        varPT_array[iev] = varPT/meanPT_array[iev]**2.
        sknewPT_array[iev] = (np.mean(sknew_dPT[array_idx])
                              /np.mean(N3_weight[array_idx])/varPT**1.5)
        kurtosisPT_array[iev] = (np.mean(kurtosis_dPT[array_idx])
                                 /np.mean(N4_weight[array_idx])/varPT**2)

    meanPTMean, meanPTErr = computeJKMeanandErr(meanPT_array)
    varPTMean, varPTErr = computeJKMeanandErr(varPT_array)
    sknewPTMean, sknewPTErr = computeJKMeanandErr(sknewPT_array)
    kurtosisPTMean, kurtosisPTErr = computeJKMeanandErr(kurtosisPT_array)

    pTfluctResults = [meanPTMean, meanPTErr, varPTMean, varPTErr,
                      sknewPTMean, sknewPTErr, kurtosisPTMean, kurtosisPTErr]

    dN_mean = np.real(np.mean(dataArr1[:, 0] + dataArr2[:, 0]))
    dN_err = np.std(dataArr1[:, 0] + dataArr2[:, 0])/np.sqrt(nev)

    outputFileName = outputFileHeader + "pTFluct.dat"
    if path.isfile(outputFileName):
        f = open(outputFileName, 'a')
    else:
        f = open(outputFileName, 'w')
        f.write("# cen  Nch  <pT>  var<dpT>/<pT>^2  <dpT^3>/<dpT^2>^1.5  "
                + "<dpT^4>/<dpT^2>^2\n")
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

dNdyList = []
for event_name in data.keys():
    dNdyList.append(data[event_name]['Nch'])
dNdyList = - np.sort(-np.array(dNdyList))
print(f"Number of good events: {len(dNdyList)}")

for icen in range(len(centralityCutList) - 1):
    if centralityCutList[icen+1] < centralityCutList[icen]:
        continue
    selected_events_list = []

    dN_dy_cut_high = dNdyList[
        int(len(dNdyList)*centralityCutList[icen]/100.)
    ]
    dN_dy_cut_low = dNdyList[
        min(len(dNdyList)-1,
            int(len(dNdyList)*centralityCutList[icen+1]/100.))
    ]

    if len(dNcutList) == len(centralityCutList):
        dN_dy_cut_high = dNcutList[icen]
        dN_dy_cut_low = dNcutList[icen+1]

    for event_name in data.keys():
        if (data[event_name]['Nch'] > dN_dy_cut_low
                and data[event_name]['Nch'] <= dN_dy_cut_high):
            selected_events_list.append(event_name)

    nev = len(selected_events_list)
    if nev <= 0:
        continue

    cenLabel = (centralityCutList[icen] +
                centralityCutList[icen+1])/2.*centralityRange
    print("analysis {}%-{}% nev = {}...".format(
            centralityCutList[icen]*centralityRange,
            centralityCutList[icen+1]*centralityRange, nev))
    print("dNdy: {0:.2f} - {1:.2f}".format(dN_dy_cut_low, dN_dy_cut_high))

    QnArr1 = []
    QnArr2 = []
    QnArr3 = []
    for event_name in selected_events_list:
        QnArr1.append(data[event_name]['ALICE_eta_-0p8_-0p4'])
        QnArr2.append(data[event_name]['ALICE_eta_0p4_0p8'])
        QnArr3.append(data[event_name]['ALICE_eta_-0p4_0p4'])
    QnArr1 = np.array(QnArr1)
    QnArr2 = np.array(QnArr2)
    QnArr3 = np.array(QnArr3)
    calculate_pTfluct(QnArr1, QnArr2, QnArr3, "ALICE", cenLabel)
