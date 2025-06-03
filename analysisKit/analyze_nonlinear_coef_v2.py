#!/usr/bin/env python3

import pickle
import numpy as np
from os import path
import sys


def help_message():
    print("Usage: {0} database_file".format(sys.argv[0]))
    exit(0)


Reg_centrality_cut_list = [
    0., 5., 10., 20., 30., 40., 50., 60., 70., 80., 90., 100.
]
centralityCutList = Reg_centrality_cut_list
dNcutList = []  # pre-defined Nch cut if simulation is not minimum bias


def calculateNonLinearResponseV2_2sub(vn_data_array1, vn_data_array2,
                                      outputFileNameV2, cenLabel):
    """
        This function computes the non-linear response coefficients
        for V4 = V6L + chi_422 V2^2 + chi_413 V1*V3
    """
    nev = len(vn_data_array1[:, 0])

    V1_1 = vn_data_array1[:, 2]
    V2_1 = vn_data_array1[:, 3]

    V1_2 = vn_data_array2[:, 2]
    V2_2 = vn_data_array2[:, 3]

    chi_211_num = V2_1*(np.conj(V1_2)**2) + V2_2*(np.conj(V1_1)**2)
    V12_V12 = 2*np.real((V1_1*np.conj(V1_2))**2)

    chi_211_JK = np.zeros(nev, dtype=complex)
    for iev in range(nev):
        array_idx = [True]*nev
        array_idx[iev] = False
        array_idx = np.array(array_idx)

        num_JK1 = np.mean(chi_211_num[array_idx])
        den_JK11 = np.mean(V12_V12[array_idx])

        chi_211_JK[iev] = num_JK1/den_JK11

    # compute the real part of the response coefficients
    chi_211_mean = np.mean(np.real(chi_211_JK))
    chi_211_err = np.sqrt((nev - 1.)/nev*np.sum(
        (np.real(chi_211_JK) - chi_211_mean)**2))
    results = [chi_211_mean, chi_211_err]

    # compute sqrt<V4L^2>
    V2L_1 = V2_1 - chi_211_mean*(V1_1**2)
    V2L_2 = V2_2 - chi_211_mean*(V1_2**2)
    v2L_rms = np.sqrt(np.mean(np.real(V2L_1*np.conj(V2L_2))))
    v2L_err = np.std(np.real(V2L_1*np.conj(V2L_2)))/np.sqrt(nev)/(2*v2L_rms)
    results += [v2L_rms, v2L_err]

    # compute the imaginary part of the response coefficients
    chi_211_mean = np.mean(np.imag(chi_211_JK))
    chi_211_err = np.sqrt((nev - 1.)/nev*np.sum(
        (np.imag(chi_211_JK) - chi_211_mean)**2))
    results += [chi_211_mean, chi_211_err]

    dN_mean = np.real(np.mean(vn_data_array1[:, 0] + vn_data_array2[:, 0]))
    dN_err = (np.std(vn_data_array1[:, 0] + vn_data_array2[:, 0])/np.sqrt(nev))
    if path.isfile(outputFileNameV2):
        f = open(outputFileNameV2, 'a')
    else:
        f = open(outputFileNameV2, 'w')
        f.write("# cen  Nch  Re{chi_211}  v2L_rms  Im{chi_211}\n")
    f.write("{:.3f}  {:.5e}  {:.5e}".format(cenLabel, dN_mean, dN_err))
    for ires in results:
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
print("Number of good events: {}".format(len(dNdyList)))

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

    cenLabel = (centralityCutList[icen] + centralityCutList[icen + 1])/2.
    print("analysis {}%-{}% nev = {}...".format(centralityCutList[icen],
                                                centralityCutList[icen + 1],
                                                nev))
    print("dNdy: {0:.2f} - {1:.2f}".format(dN_dy_cut_low, dN_dy_cut_high))

    QnArr1 = []
    QnArr2 = []
    QnArr3 = []
    for event_name in selected_events_list:
        QnArr1.append(data[event_name]['ALICE_eta_-0p4_0p4'])
        QnArr2.append(data[event_name]['ALICE_eta_-3p2_-0p4'])
        QnArr3.append(data[event_name]['ALICE_eta_0p4_3p2'])
    QnArr1 = np.array(QnArr1)
    QnArr2 = np.array(QnArr2)
    QnArr3 = np.array(QnArr3)

    calculateNonLinearResponseV2_2sub(QnArr2, QnArr3,
                                      "nonLinearV2_2sub_full.dat", cenLabel)
