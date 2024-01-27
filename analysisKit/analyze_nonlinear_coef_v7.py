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


def calculate_chi_422(vn_data_array1, vn_data_array2):
    """
        this function computes the non-linear response coefficients
            chi_422 = Re<V_4*conj(V_2)^2>/(<|V_2|^4>)
        we use one flow vector for V_n and the other for conj(V_n)
    """
    dN1 = np.real(vn_data_array1[:, -1])
    dN2 = np.real(vn_data_array2[:, -1])

    Q2_1 = dN1*vn_data_array1[:, 3]
    Q4_1 = dN1*vn_data_array1[:, 5]
    Q2_2 = dN2*vn_data_array2[:, 3]
    Q4_2 = dN2*vn_data_array2[:, 5]

    # four-particle correlation
    N4_weight = dN1*(dN1 - 1)*dN2*(dN2 - 1)
    Q2_N4 = (np.real(Q2_1*np.conj(Q2_2)*Q2_1*np.conj(Q2_2))
             - np.real(Q4_1*np.conj(Q2_2)*np.conj(Q2_2))
             - np.real(Q2_1*Q2_1*np.conj(Q4_2))
             + np.real(Q4_1*np.conj(Q4_2)))

    # three-particle correlation
    N3_weight = dN1*dN2*(dN2 - 1) + dN1*(dN1 - 1)*dN2
    chi_422_num = (Q4_1*np.conj(Q2_2)*np.conj(Q2_2) - Q4_1*np.conj(Q4_2)
                   + Q4_2*np.conj(Q2_1)*np.conj(Q2_1) - Q4_2*np.conj(Q4_1))

    num = np.real(np.mean(chi_422_num))/np.mean(N3_weight)
    den = np.mean(Q2_N4)/np.mean(N4_weight)
    chi_422_mean = num/den
    return chi_422_mean


def calculate_chi_523(vn_data_array1, vn_data_array2):
    """
        this function computes the non-linear response coefficients
            chi_523 = Re<V_5*conj(V_2*V_3)>/(<|V_2|^2||V_3|^2>)
        we use one flow vector for V_n and the other for conj(V_n)
    """
    dN1 = np.real(vn_data_array1[:, -1])
    dN2 = np.real(vn_data_array2[:, -1])

    Q2_1 = dN1*vn_data_array1[:, 3]
    Q3_1 = dN1*vn_data_array1[:, 4]
    Q5_1 = dN1*vn_data_array1[:, 6]

    Q2_2 = dN2*vn_data_array2[:, 3]
    Q3_2 = dN2*vn_data_array2[:, 4]
    Q5_2 = dN2*vn_data_array2[:, 6]

    N4_weight = dN1*(dN1 - 1)*dN2*(dN2 - 1)
    Q_32 = (np.real(Q2_1*Q3_1*np.conj(Q2_2*Q3_2))
            - np.real(Q5_1*np.conj(Q2_2*Q3_2))
            - np.real(Q2_1*Q3_1*np.conj(Q5_2))
            + np.real(Q5_1*np.conj(Q5_2)))

    N3_weight = dN1*dN2*(dN2 - 1) + dN1*(dN1 - 1)*dN2
    chi_523_num = (Q5_1*np.conj(Q2_2*Q3_2) - Q5_1*np.conj(Q5_2)
                   + Q5_2*np.conj(Q2_1*Q3_1) - Q5_2*np.conj(Q5_1))

    num = np.real(np.mean(chi_523_num))/np.mean(N3_weight)
    den = np.mean(Q_32)/np.mean(N4_weight)
    chi_523_mean = num/den
    return chi_523_mean


def calculateNonLinearResponseV7_2sub(vn_data_array1, vn_data_array2,
                                      outputFileNameV7, flag, cenLabel):
    """
        This function computes the non-linear response coefficients
        for V7 = V7L + chi_7223 V2^2*V3 + chi_734 V3*V4L + chi_725 V2*V5L
    """
    nev = len(vn_data_array1[:, 0])

    V2_1 = vn_data_array1[:, 3]
    V3_1 = vn_data_array1[:, 4]
    V4_1 = vn_data_array1[:, 5]
    V5_1 = vn_data_array1[:, 6]
    V7_1 = vn_data_array1[:, 8]

    V2_2 = vn_data_array2[:, 3]
    V3_2 = vn_data_array2[:, 4]
    V4_2 = vn_data_array2[:, 5]
    V5_2 = vn_data_array2[:, 6]
    V7_2 = vn_data_array2[:, 8]

    chi_422 = calculate_chi_422(vn_data_array1, vn_data_array2)
    V4L_1 = V4_1 - chi_422*(V2_1**2)
    V4L_2 = V4_2 - chi_422*(V2_2**2)

    chi_523 = calculate_chi_523(vn_data_array1, vn_data_array2)
    V5L_1 = V5_1 - chi_523*(V2_1*V3_1)
    V5L_2 = V5_2 - chi_523*(V2_2*V3_2)

    chi_7223_num = V7_1*np.conj((V2_2**2)*V3_2) + V7_2*np.conj((V2_1**2)*V3_1)
    V22V3_V22V3 = 2*np.real((V2_1**2)*V3_1*np.conj((V2_2**2)*V3_2))
    V3V4_V22V3 = (V3_1*V4L_1*np.conj((V2_2**2)*V3_2)
                  + V3_2*V4L_2*np.conj((V2_1**2)*V3_1))
    V2V5_V22V3 = (V2_1*V5L_1*np.conj((V2_2**2)*V3_2)
                  + V2_2*V5L_2*np.conj((V2_1**2)*V3_1))

    chi_734_num = V7_1*np.conj(V3_2*V4L_2) + V7_2*np.conj(V3_1*V4L_1)
    V3V4_V3V4 = 2*np.real(V3_1*V4L_1*np.conj(V3_2*V4L_2))
    V2V5_V3V4 = V2_1*V5L_1*np.conj(V3_2*V4L_2) + V2_2*V5L_2*np.conj(V3_1*V4L_1)

    chi_725_num = V7_1*np.conj(V2_2*V5L_2) + V7_2*np.conj(V2_1*V5L_1)
    V2V5_V2V5 = 2*np.real(V2_1*V5L_1*np.conj(V2_2*V5L_2))

    chi_7223_JK = np.zeros(nev, dtype=complex)
    chi_734_JK = np.zeros(nev, dtype=complex)
    chi_725_JK = np.zeros(nev, dtype=complex)
    for iev in range(nev):
        array_idx = [True]*nev
        array_idx[iev] = False
        array_idx = np.array(array_idx)

        num_JK1 = np.mean(chi_7223_num[array_idx])
        den_JK11 = np.mean(V22V3_V22V3[array_idx])
        den_JK12 = np.mean(V3V4_V22V3[array_idx])
        den_JK13 = np.mean(V2V5_V22V3[array_idx])

        num_JK2 = np.mean(chi_734_num[array_idx])
        den_JK21 = np.conj(den_JK12)
        den_JK22 = np.mean(V3V4_V3V4[array_idx])
        den_JK23 = np.mean(V2V5_V3V4[array_idx])

        num_JK3 = np.mean(chi_725_num[array_idx])
        den_JK31 = np.conj(den_JK13)
        den_JK32 = np.conj(den_JK23)
        den_JK33 = np.mean(V2V5_V2V5[array_idx])

        array_lhs = np.array([num_JK1, num_JK2, num_JK3], dtype=np.cfloat)

        if flag == 1:
            array_rhs = np.array([[den_JK11, den_JK12, den_JK13],
                                  [den_JK21, den_JK22, den_JK23],
                                  [den_JK31, den_JK32, den_JK33]],
                                 dtype=np.cfloat)
        else:
            array_rhs = np.array([[den_JK11,        0,        0],
                                  [0,        den_JK22,        0],
                                  [0,               0, den_JK33]],
                                 dtype=np.cfloat)
        chi_7 = np.linalg.solve(array_rhs, array_lhs)
        chi_7223_JK[iev] = chi_7[0]
        chi_734_JK[iev] = chi_7[1]
        chi_725_JK[iev] = chi_7[2]

    # compute the real part of the response coefficients
    chi_7223_mean = np.mean(np.real(chi_7223_JK))
    chi_7223_err = np.sqrt((nev - 1.)/nev
                           * np.sum((np.real(chi_7223_JK) - chi_7223_mean)**2))
    chi_734_mean = np.mean(np.real(chi_734_JK))
    chi_734_err = np.sqrt((nev - 1.)/nev
                          * np.sum((np.real(chi_734_JK) - chi_734_mean)**2.))
    chi_725_mean = np.mean(np.real(chi_725_JK))
    chi_725_err = np.sqrt((nev - 1.)/nev
                          * np.sum((np.real(chi_725_JK) - chi_725_mean)**2.))
    results = [chi_7223_mean, chi_7223_err, chi_734_mean, chi_734_err,
               chi_725_mean, chi_725_err]

    # compute sqrt<V7L^2>
    V7L_1 = (V7_1 - chi_7223_mean*(V2_1**2*V3_1) - chi_734_mean*(V3_1*V4L_1)
             - chi_725_mean*(V2_1*V5L_1))
    V7L_2 = (V7_2 - chi_7223_mean*(V2_2**2*V3_2) - chi_734_mean*(V3_2*V4L_2)
             - chi_725_mean*(V2_2*V5L_2))
    v7L_rms = np.nan_to_num(np.sqrt(np.mean(np.real(V7L_1*np.conj(V7L_2)))))
    v7L_err = np.nan_to_num(np.std(np.real(V7L_1*np.conj(V7L_2)))/np.sqrt(nev)
                            / (2*v7L_rms))

    results += [v7L_rms, v7L_err]
    # compute the imaginary part of the response coefficients
    chi_7223_mean = np.mean(np.imag(chi_7223_JK))
    chi_7223_err = np.sqrt((nev - 1.)/nev
                           * np.sum((np.imag(chi_7223_JK) - chi_7223_mean)**2))
    chi_734_mean = np.mean(np.imag(chi_734_JK))
    chi_734_err = np.sqrt((nev - 1.)/nev
                          * np.sum((np.imag(chi_734_JK) - chi_734_mean)**2.))
    chi_725_mean = np.mean(np.imag(chi_725_JK))
    chi_725_err = np.sqrt((nev - 1.)/nev
                          * np.sum((np.imag(chi_725_JK) - chi_725_mean)**2.))
    results += [chi_7223_mean, chi_7223_err, chi_734_mean, chi_734_err,
                chi_725_mean, chi_725_err]

    dN_mean = np.real(np.mean(vn_data_array1[:, 0] + vn_data_array2[:, 0]))
    dN_err = (np.std(vn_data_array1[:, 0] + vn_data_array2[:, 0])
              / np.sqrt(nev))
    if path.isfile(outputFileNameV7):
        f = open(outputFileNameV7, 'a')
    else:
        f = open(outputFileNameV7, 'w')
        f.write("# cen  Nch  Re{chi_7223}  Re{chi_734}  Re{chi_725}  "
                + "v7L_rms  Im{chi_7223}  Im{chi_734}  Im{chi_725}\n")
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

    calculateNonLinearResponseV7_2sub(QnArr2, QnArr3,
                                      "nonLinearV7_2sub_diag.dat", 0, cenLabel)
    calculateNonLinearResponseV7_2sub(QnArr2, QnArr3,
                                      "nonLinearV7_2sub_full.dat", 1, cenLabel)
