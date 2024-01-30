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


def calculateNonLinearResponseV6_2sub(vn_data_array1, vn_data_array2,
                                      outputFileNameV6, flag, cenLabel):
    """
        This function computes the non-linear response coefficients
        for V6 = V6L + chi_6222 V2^3 + chi_633 V3^2 + chi_624 V2*V4L
    """
    nev = len(vn_data_array1[:, 0])

    V2_1 = vn_data_array1[:, 3]
    V3_1 = vn_data_array1[:, 4]
    V4_1 = vn_data_array1[:, 5]
    V6_1 = vn_data_array1[:, 7]
    V2_2 = vn_data_array2[:, 3]
    V3_2 = vn_data_array2[:, 4]
    V4_2 = vn_data_array2[:, 5]
    V6_2 = vn_data_array2[:, 7]

    chi_422 = calculate_chi_422(vn_data_array1, vn_data_array2)
    V4L_1 = V4_1 - chi_422*(V2_1**2)
    V4L_2 = V4_2 - chi_422*(V2_2**2)

    chi_6222_num = V6_1*(np.conj(V2_2)**3) + V6_2*(np.conj(V2_1)**3)
    V2_6 = 2*np.real((V2_1*np.conj(V2_2))**3)
    V32_V23 = (V3_1**2)*(np.conj(V2_2)**3) + (V3_2**2)*(np.conj(V2_1)**3)
    V2V4L_V23 = ((V2_1*V4L_1)*(np.conj(V2_2)**3)
                 + (V2_2*V4L_2)*(np.conj(V2_1)**3))

    chi_633_num = V6_1*(np.conj(V3_2)**2) + V6_2*(np.conj(V3_1)**2)
    V3_4 = 2*np.real((V3_1*np.conj(V3_2))**2)
    V2V4L_V32 = ((V2_1*V4L_1)*(np.conj(V3_2)**2)
                 + (V2_2*V4L_2)*(np.conj(V3_1)**2))

    v624_num = V6_1*(np.conj(V2_2*V4L_2)) + V6_2*(np.conj(V2_1*V4L_1))
    V2V4L_2 = 2*np.real(V2_1*np.conj(V2_2)*V4L_1*np.conj(V4L_2))

    chi_6222_JK = np.zeros(nev, dtype=complex)
    chi_633_JK = np.zeros(nev, dtype=complex)
    chi_624_JK = np.zeros(nev, dtype=complex)
    for iev in range(nev):
        array_idx = [True]*nev
        array_idx[iev] = False
        array_idx = np.array(array_idx)

        num_JK1 = np.mean(chi_6222_num[array_idx])
        den_JK11 = np.mean(V2_6[array_idx])
        den_JK12 = np.mean(V32_V23[array_idx])
        den_JK13 = np.mean(V2V4L_V23[array_idx])

        num_JK2 = np.mean(chi_633_num[array_idx])
        den_JK21 = np.conj(den_JK12)
        den_JK22 = np.mean(V3_4[array_idx])
        den_JK23 = np.mean(V2V4L_V32[array_idx])

        num_JK3 = np.mean(v624_num[array_idx])
        den_JK31 = np.conj(den_JK13)
        den_JK32 = np.conj(den_JK23)
        den_JK33 = np.mean(V2V4L_2[array_idx])

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

        chi_6 = np.linalg.solve(array_rhs, array_lhs)
        chi_6222_JK[iev] = chi_6[0]
        chi_633_JK[iev] = chi_6[1]
        chi_624_JK[iev] = chi_6[2]

    # compute the real part of the response coefficients
    chi_6222_mean = np.mean(np.real(chi_6222_JK))
    chi_6222_err = np.sqrt((nev - 1.)/nev
                           * np.sum((np.real(chi_6222_JK) - chi_6222_mean)**2))
    chi_633_mean = np.mean(np.real(chi_633_JK))
    chi_633_err = np.sqrt((nev - 1.)/nev
                          * np.sum((np.real(chi_633_JK) - chi_633_mean)**2.))
    chi_624_mean = np.mean(np.real(chi_624_JK))
    chi_624_err = np.sqrt((nev - 1.)/nev
                          * np.sum((np.real(chi_624_JK) - chi_624_mean)**2.))
    results = [chi_6222_mean, chi_6222_err, chi_633_mean, chi_633_err,
               chi_624_mean, chi_624_err]

    # compute sqrt<V6L^2>
    V6L_1 = (V6_1 - chi_6222_mean*(V2_1**3) - chi_633_mean*(V3_1**2)
             - chi_624_mean*(V2_1*V4L_1))
    V6L_2 = (V6_2 - chi_6222_mean*(V2_2**3) - chi_633_mean*(V3_2**2)
             - chi_624_mean*(V2_2*V4L_2))
    v6L_rms = np.sqrt(np.mean(np.real(V6L_1*np.conj(V6L_2))))
    v6L_err = np.std(np.real(V6L_1*np.conj(V6L_2)))/np.sqrt(nev)/(2*v6L_rms)
    results += [v6L_rms, v6L_err]

    # compute the imaginary part of the response coefficients
    chi_6222_mean = np.mean(np.imag(chi_6222_JK))
    chi_6222_err = np.sqrt((nev - 1.)/nev
                           * np.sum((np.imag(chi_6222_JK) - chi_6222_mean)**2))
    chi_633_mean = np.mean(np.imag(chi_633_JK))
    chi_633_err = np.sqrt((nev - 1.)/nev
                          * np.sum((np.imag(chi_633_JK) - chi_633_mean)**2.))
    chi_624_mean = np.mean(np.imag(chi_624_JK))
    chi_624_err = np.sqrt((nev - 1.)/nev
                          * np.sum((np.imag(chi_624_JK) - chi_624_mean)**2.))
    results += [chi_6222_mean, chi_6222_err, chi_633_mean, chi_633_err,
                chi_624_mean, chi_624_err]

    dN_mean = np.real(np.mean(vn_data_array1[:, 0] + vn_data_array2[:, 0]))
    dN_err = (np.std(vn_data_array1[:, 0] + vn_data_array2[:, 0])
              / np.sqrt(nev))
    if path.isfile(outputFileNameV6):
        f = open(outputFileNameV6, 'a')
    else:
        f = open(outputFileNameV6, 'w')
        f.write("# cen  Nch  Re{chi_6222}  Re{chi_633}  Re{chi_624}  "
                + "v6L_rms  Im{chi_6222}  Im{chi_633}  Im{chi_624}\n")
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

    calculateNonLinearResponseV6_2sub(QnArr2, QnArr3,
                                      "nonLinearV6_2sub_diag.dat", 0, cenLabel)
    calculateNonLinearResponseV6_2sub(QnArr2, QnArr3,
                                      "nonLinearV6_2sub_full.dat", 1, cenLabel)
