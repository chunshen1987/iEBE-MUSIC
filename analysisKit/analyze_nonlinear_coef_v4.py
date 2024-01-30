#!/usr/bin/env python3

import pickle
import numpy as np
from os import path
import sys


def help_message():
    print("Usage: {0} database_file".format(sys.argv[0]))
    exit(0)


Reg_centrality_cut_list = [0., 5., 10., 20., 30., 40., 50.,
                           60., 70., 80., 90., 100.]
centralityCutList = Reg_centrality_cut_list
dNcutList = []    # pre-defined Nch cut if simulation is not minimum bias


def calculateNonLinearResponseV4_2sub(vn_data_array1, vn_data_array2,
                                      outputFileNameV4, flag, cenLabel):
    """
        This function computes the non-linear response coefficients
        for V4 = V6L + chi_422 V2^2 + chi_413 V1*V3
    """
    nev = len(vn_data_array1[:, 0])
    #dN1 = np.real(vn_data_array1[:, -1])
    #dN2 = np.real(vn_data_array2[:, -1])

    V1_1 = vn_data_array1[:, 2]
    V2_1 = vn_data_array1[:, 3]
    V3_1 = vn_data_array1[:, 4]
    V4_1 = vn_data_array1[:, 5]

    V1_2 = vn_data_array2[:, 2]
    V2_2 = vn_data_array2[:, 3]
    V3_2 = vn_data_array2[:, 4]
    V4_2 = vn_data_array2[:, 5]

    chi_422_num = V4_1*(np.conj(V2_2)**2) + V4_2*(np.conj(V2_1)**2)
    V22_V22 = 2*np.real((V2_1*np.conj(V2_2))**2)
    V1V3_V22 = (V1_1*V3_1)*(np.conj(V2_2)**2) + (V1_2*V3_2)*(np.conj(V2_2)**2)

    chi_413_num = V4_1*np.conj(V1_2*V3_2) + V4_2*np.conj(V1_1*V3_1)
    V1V3_V1V3 = 2*np.real(V1_1*V3_1*np.conj(V1_2*V3_2))

    chi_422_JK = np.zeros(nev, dtype=complex)
    chi_413_JK = np.zeros(nev, dtype=complex)
    for iev in range(nev):
        array_idx = [True]*nev
        array_idx[iev] = False
        array_idx = np.array(array_idx)

        num_JK1 = np.mean(chi_422_num[array_idx])
        den_JK11 = np.mean(V22_V22[array_idx])
        den_JK12 = np.mean(V1V3_V22[array_idx])

        num_JK2 = np.mean(chi_413_num[array_idx])
        den_JK21 = np.conj(den_JK12)
        den_JK22 = np.mean(V1V3_V1V3[array_idx])

        array_lhs = np.array([num_JK1, num_JK2], dtype=np.cfloat)

        if flag == 1:
            array_rhs = np.array([[den_JK11, den_JK12],
                                  [den_JK21, den_JK22]],
                                 dtype=np.cfloat)
        else:
            array_rhs = np.array([[den_JK11,        0],
                                  [0,        den_JK22]],
                                 dtype=np.cfloat)

        chi_6 = np.linalg.solve(array_rhs, array_lhs)
        chi_422_JK[iev] = chi_6[0]
        chi_413_JK[iev] = chi_6[1]

    # compute the real part of the response coefficients
    chi_422_mean = np.mean(np.real(chi_422_JK))
    chi_422_err = np.sqrt((nev - 1.)/nev
                          * np.sum((np.real(chi_422_JK) - chi_422_mean)**2))
    chi_413_mean = np.mean(np.real(chi_413_JK))
    chi_413_err = np.sqrt((nev - 1.)/nev
                          * np.sum((np.real(chi_413_JK) - chi_413_mean)**2.))
    results = [chi_422_mean, chi_422_err, chi_413_mean, chi_413_err]

    # compute sqrt<V4L^2>
    V4L_1 = V4_1 - chi_422_mean*(V2_1**2) - chi_413_mean*(V1_1*V3_1)
    V4L_2 = V4_2 - chi_422_mean*(V2_2**2) - chi_413_mean*(V1_2*V3_2)
    v4L_rms = np.sqrt(np.mean(np.real(V4L_1*np.conj(V4L_2))))
    v4L_err = np.std(np.real(V4L_1*np.conj(V4L_2)))/np.sqrt(nev)/(2*v4L_rms)
    results += [v4L_rms, v4L_err]

    # compute the imaginary part of the response coefficients
    chi_422_mean = np.mean(np.imag(chi_422_JK))
    chi_422_err = np.sqrt((nev - 1.)/nev
                          * np.sum((np.imag(chi_422_JK) - chi_422_mean)**2))
    chi_413_mean = np.mean(np.imag(chi_413_JK))
    chi_413_err = np.sqrt((nev - 1.)/nev
                          * np.sum((np.imag(chi_413_JK) - chi_413_mean)**2.))
    results += [chi_422_mean, chi_422_err, chi_413_mean, chi_413_err]

    dN_mean = np.real(np.mean(vn_data_array1[:, 0] + vn_data_array2[:, 0]))
    dN_err = (np.std(vn_data_array1[:, 0] + vn_data_array2[:, 0])
              / np.sqrt(nev))
    if path.isfile(outputFileNameV4):
        f = open(outputFileNameV4, 'a')
    else:
        f = open(outputFileNameV4, 'w')
        f.write("# cen  Nch  Re{chi_422}  Re{chi_413}  "
                + "v4L_rms  Im{chi_422}  Im{chi_413}\n")
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
dNdyList = - np.sort(-np.array(dNdyList))
print("Number of good events: {}".format(len(dNdyList)))

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
                centralityCutList[icen+1])/2.
    print("analysis {}%-{}% nev = {}...".format(
            centralityCutList[icen], centralityCutList[icen+1], nev))
    print("dNdy: {0:.2f} - {1:.2f}".format(dN_dy_cut_low, dN_dy_cut_high))

    QnArr1 = []
    QnArr2 = []
    QnArr3 = []
    for event_name in selected_events_list:
        #QnArr1.append(data[event_name]['STAR_eta_-0p5_0p5'])
        #QnArr2.append(data[event_name]['STAR_eta_-1_-0p5'])
        #QnArr3.append(data[event_name]['STAR_eta_0p5_1'])
        QnArr1.append(data[event_name]['ALICE_eta_-0p4_0p4'])
        QnArr2.append(data[event_name]['ALICE_eta_-0p8_-0p4'])
        QnArr3.append(data[event_name]['ALICE_eta_0p4_0p8'])
    QnArr1 = np.array(QnArr1)
    QnArr2 = np.array(QnArr2)
    QnArr3 = np.array(QnArr3)

    calculateNonLinearResponseV4_2sub(QnArr2, QnArr3,
                                      "nonLinearV4_2sub_diag.dat", 0, cenLabel)
    calculateNonLinearResponseV4_2sub(QnArr2, QnArr3,
                                      "nonLinearV4_2sub_full.dat", 1, cenLabel)
