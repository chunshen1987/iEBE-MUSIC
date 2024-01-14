#!/usr/bin/env python3

import sys
from os import path
import pickle
import numpy as np


def help_message():
    print("Usage: {0} database_file".format(sys.argv[0]))
    exit(0)


centralityRange = 1.
pidList = ['pi+', 'pi-', 'K+', 'K-', 'p', 'pbar']
Reg_centrality_cut_list = [0., 5., 10., 20., 30., 40., 50.,
                           60., 70., 80., 90., 100.]
centralityCutList = Reg_centrality_cut_list
# centralityCutList = [0, 1, 2, 3, 4, 6, 8, 10, 15, 20, 30, 40, 50, 60,
#                      70, 80, 90, 100]
dNcutList = []    # pre-defined Nch cut if simulation is not minimum bias


def calculate_pid_dN(dN_data_array, outputFilename, cenLabel):
    """
        This function computes the averaged value of the identified particle
        yields for different centralities.
    """
    nev = len(dN_data_array[:, 0])
    dN_mean = np.mean(dN_data_array, axis=0)
    dN_err = np.std(dN_data_array, axis=0)/np.sqrt(nev)

    if path.isfile(outputFilename):
        f = open(outputFilename, 'a')
    else:
        f = open(outputFilename, 'w')
        f.write("# cen  ch  pi+  pi-  K+  K-  p  pbar\n")
    f.write("{:.3f}".format(cenLabel))
    for i in range(len(dN_mean)):
        f.write("  {:.5e}  {:.5e}".format(dN_mean[i], dN_err[i]))
    f.write("\n")
    f.close()
    return


def calculate_pid_meanpT(pT_data_array, outputFilename, cenLabel):
    """
        This function computes the averaged value of the identified particle
        mean pT for different centralities.
    """
    nev = len(pT_data_array[:, 0])
    meanpT_mean = np.mean(pT_data_array, axis=0)
    meanpT_err = np.std(pT_data_array, axis=0)/np.sqrt(nev)

    if path.isfile(outputFilename):
        f = open(outputFilename, 'a')
    else:
        f = open(outputFilename, 'w')
        f.write("# cen  ch  pi+  pi-  K+  K-  p  pbar\n")
    f.write("{:.3f}".format(cenLabel))
    for i in range(len(meanpT_mean)):
        f.write("  {:.5e}  {:.5e}".format(meanpT_mean[i], meanpT_err[i]))
    f.write("\n")
    f.close()
    return


def calculateSymmetricCumulant2sub(vn_data_array1, vn_data_array2,
                                   outputFileName, cenLabel):
    """
        this funciton computes the symmetric cumulant
            SC(m,n) = <v_m*conj(v_m)*v_n*conj(v_n)>
                      - <v_m*conj(v_m)>*<v_n*conj(v_n)>
        we use one flow vecotor for v_{n,m} and another for conj(v_{n,m})
        we use Jackknife resampling method to estimate the statistical error
    """
    nev = len(vn_data_array1[:, 0])
    dN1 = np.real(vn_data_array1[:, -1])
    dN2 = np.real(vn_data_array2[:, -1])

    Q2_1 = dN1*vn_data_array1[:, 3]
    Q3_1 = dN1*vn_data_array1[:, 4]
    Q4_1 = dN1*vn_data_array1[:, 5]
    Q5_1 = dN1*vn_data_array1[:, 6]
    Q6_1 = dN1*vn_data_array1[:, 7]
    Q2_2 = dN2*vn_data_array2[:, 3]
    Q3_2 = dN2*vn_data_array2[:, 4]
    Q4_2 = dN2*vn_data_array2[:, 5]
    Q5_2 = dN2*vn_data_array2[:, 6]
    Q6_2 = dN2*vn_data_array2[:, 7]

    # two-particle correlation
    N2_weight = dN1*dN2
    Q2_N2 = np.real(Q2_1*np.conj(Q2_2))
    Q3_N2 = np.real(Q3_1*np.conj(Q3_2))
    Q4_N2 = np.real(Q4_1*np.conj(Q4_2))

    # four-particle correlation
    N4_weight = dN1*(dN1 - 1)*dN2*(dN2 - 1)
    Q2Q3_N4 = (np.real(Q2_1*np.conj(Q2_2)*Q3_1*np.conj(Q3_2))
               - np.real(Q5_1*np.conj(Q2_2)*np.conj(Q3_2))
               - np.real(Q2_1*Q3_1*np.conj(Q5_2))
               + np.real(Q5_1*np.conj(Q5_2)))
    Q2Q4_N4 = (np.real(Q2_1*np.conj(Q2_2)*Q4_1*np.conj(Q4_2))
               - np.real(Q6_1*np.conj(Q2_2)*np.conj(Q4_2))
               - np.real(Q2_1*Q4_1*np.conj(Q6_2))
               + np.real(Q6_1*np.conj(Q6_2)))

    # calcualte observables with Jackknife resampling method
    SC32_array = np.zeros(nev)
    SC42_array = np.zeros(nev)
    NSC32_array = np.zeros(nev)
    NSC42_array = np.zeros(nev)
    for iev in range(nev):
        array_idx = [True]*nev
        array_idx[iev] = False
        array_idx = np.array(array_idx)

        # SC(3,2)
        v2v3 = ((np.mean(Q3_N2[array_idx])*np.mean(Q2_N2[array_idx]))
                / (np.mean(N2_weight[array_idx])**2.))
        SC32_array[iev] = (np.mean(Q2Q3_N4[array_idx])
                           / np.mean(N4_weight[array_idx]) - v2v3)
        NSC32_array[iev] = SC32_array[iev]/v2v3

        # SC(4,2)
        v2v4 = ((np.mean(Q4_N2[array_idx])*np.mean(Q2_N2[array_idx]))
                / (np.mean(N2_weight[array_idx])**2.))
        SC42_array[iev] = (np.mean(Q2Q4_N4[array_idx])
                           / np.mean(N4_weight[array_idx]) - v2v4)
        NSC42_array[iev] = SC42_array[iev]/v2v4

    SC32_mean = np.mean(SC32_array)
    SC32_err = np.sqrt((nev - 1.)/nev*np.sum((SC32_array - SC32_mean)**2.))
    NSC32_mean = np.mean(NSC32_array)
    NSC32_err = np.sqrt((nev - 1.)/nev*np.sum((NSC32_array - NSC32_mean)**2.))

    SC42_mean = np.mean(SC42_array)
    SC42_err = np.sqrt((nev - 1.)/nev*np.sum((SC42_array - SC42_mean)**2.))
    NSC42_mean = np.mean(NSC42_array)
    NSC42_err = np.sqrt((nev - 1.)/nev*np.sum((NSC42_array - NSC42_mean)**2.))

    results = [SC32_mean, SC32_err, NSC32_mean, NSC32_err,
               SC42_mean, SC42_err, NSC42_mean, NSC42_err]
    dN_mean = np.real(np.mean(vn_data_array1[:, 0] + vn_data_array2[:, 0]))
    dN_err = (np.std(vn_data_array1[:, 0] + vn_data_array2[:, 0])
              / np.sqrt(nev))
    if path.isfile(outputFileName):
        f = open(outputFileName, 'a')
    else:
        f = open(outputFileName, 'w')
        f.write("# cen  Nch  SC(23)  NSC(23)  SC(24)  NSC(24)\n")
    f.write("{:.3f}  {:.5e}  {:.5e}".format(cenLabel, dN_mean, dN_err))
    for ires in results:
        f.write("  {:.5e}".format(ires))
    f.write("\n")
    f.close()


def calculateNonLinearResponseV2_2sub(vn_data_array1, vn_data_array2,
                                      outputFileName, cenLabel):
    """
        this function computes the non-linear response coefficients

        chi_211 = Re<V_2*conj(V_1)^2>/(<|V_1|^4>)
        v_211 = Re<V_2*conj(V_1)^2>/sqrt(<|V_1|^4>)
        rho_211 = Re<V_2*conj(V_1)^2>/sqrt(<|V_2|^2*<|V_1|^4>)
        v2_L = sqrt(|V_2|^2 - Re<V_2*conj(V_1)^2>)

        we use one flow vector for V_n and the other for conj(V_n)
        we will use Jackknife resampling method to estimate
        the statistical error
    """
    nev = len(vn_data_array1[:, 0])
    dN1 = np.real(vn_data_array1[:, -1])
    dN2 = np.real(vn_data_array2[:, -1])

    Q1_1 = dN1*vn_data_array1[:, 2]
    Q2_1 = dN1*vn_data_array1[:, 3]

    Q1_2 = dN2*vn_data_array2[:, 2]
    Q2_2 = dN2*vn_data_array2[:, 3]

    # two-particle correlation
    N2_weight = dN1*dN2
    Q2_N2 = np.real(Q2_1*np.conj(Q2_2))

    # four-particle correlation
    N4_weight = dN1*(dN1 - 1)*dN2*(dN2 - 1)
    Q1_N4 = (np.real(Q1_1*np.conj(Q1_2)*Q1_1*np.conj(Q1_2))
             - np.real(Q2_1*np.conj(Q1_2)*np.conj(Q1_2))
             - np.real(Q1_1*Q1_1*np.conj(Q2_2))
             + np.real(Q2_1*np.conj(Q2_2)))

    # three-particle correlation
    N3_weight = dN1*dN2*(dN2 - 1) + dN1*(dN1 - 1)*dN2
    chi_211_num = (Q2_1*np.conj(Q1_2)*np.conj(Q1_2) - Q2_1*np.conj(Q2_2)
                   + Q2_2*np.conj(Q1_1)*np.conj(Q1_1) - Q2_2*np.conj(Q2_1))

    chi_211_JK = np.zeros(nev)
    v211_JK = np.zeros(nev)
    rho211_JK = np.zeros(nev)
    v2L_JK = np.zeros(nev)
    for iev in range(nev):
        array_idx = [True]*nev
        array_idx[iev] = False
        array_idx = np.array(array_idx)

        num_JK = (np.real(np.mean(chi_211_num[array_idx]))
                  / np.mean(N3_weight[array_idx]))
        den_JK = (np.real(np.mean(Q1_N4[array_idx]))
                  / np.mean(N4_weight[array_idx]))

        v2_Psi2 = np.nan_to_num(np.sqrt(np.mean(Q2_N2[array_idx])
                                        / np.mean(N2_weight[array_idx])))

        chi_211_JK[iev] = num_JK/den_JK
        v211_JK[iev] = np.nan_to_num(num_JK/np.sqrt(den_JK))
        rho211_JK[iev] = v211_JK[iev]/v2_Psi2
        v2L_JK[iev] = np.nan_to_num(np.sqrt(v2_Psi2**2 - v211_JK[iev]**2.))

    chi_211_mean = np.mean(chi_211_JK)
    chi_211_err = np.sqrt((nev - 1.)/nev
                          * np.sum((chi_211_JK - chi_211_mean)**2.))
    v211_mean = np.mean(v211_JK)
    v211_err = np.sqrt((nev - 1.)/nev*np.sum((v211_JK - v211_mean)**2.))
    rho211_mean = np.mean(rho211_JK)
    rho211_err = np.sqrt((nev - 1.)/nev*np.sum((rho211_JK - rho211_mean)**2.))
    v2L_mean = np.mean(v2L_JK)
    v2L_err = np.sqrt((nev - 1.)/nev*np.sum((v2L_JK - v2L_mean)**2.))
    results = [v2L_mean, v2L_err, v211_mean, v211_err, rho211_mean, rho211_err,
               chi_211_mean, chi_211_err]

    dN_mean = np.real(np.mean(vn_data_array1[:, 0] + vn_data_array2[:, 0]))
    dN_err = (np.std(vn_data_array1[:, 0] + vn_data_array2[:, 0])
              / np.sqrt(nev))
    if path.isfile(outputFileName):
        f = open(outputFileName, 'a')
    else:
        f = open(outputFileName, 'w')
        f.write("# cen  Nch  v2L  v_211  rho_211  chi211\n")
    f.write("{:.3f}  {:.5e}  {:.5e}".format(cenLabel, dN_mean, dN_err))
    for ires in results:
        f.write("  {:.5e}".format(ires))
    f.write("\n")
    f.close()


def calculateNonLinearResponseV3_2sub(vn_data_array1, vn_data_array2,
                                      outputFileName, cenLabel):
    """
        this function computes the non-linear response coefficients

        chi_312 = Re<V_3*conj(V_1*V_2)>/(<|V_1|^2|V_2|^2>)
        v_312 = Re<V_3*conj(V_1*V_2)>/sqrt(<|V_1|^2|V_2|^2>)
        rho_312 = Re<V_3*conj(V_1*V_2)>/sqrt(<|V_3|^2*<|V_1|^2*|V_2|^2>)
        v3_L = sqrt(|V_3|^2 - Re<V_3*conj(V_1*V_2)>)

        we use one flow vector for V_n and the other for conj(V_n)
        we will use Jackknife resampling method to estimate
        the statistical error
    """
    nev = len(vn_data_array1[:, 0])
    dN1 = np.real(vn_data_array1[:, -1])
    dN2 = np.real(vn_data_array2[:, -1])

    Q1_1 = dN1*vn_data_array1[:, 2]
    Q2_1 = dN1*vn_data_array1[:, 3]
    Q3_1 = dN1*vn_data_array1[:, 4]

    Q1_2 = dN2*vn_data_array2[:, 2]
    Q2_2 = dN2*vn_data_array2[:, 3]
    Q3_2 = dN2*vn_data_array2[:, 4]

    # two-particle correlation
    N2_weight = dN1*dN2
    Q3_N2 = np.real(Q3_1*np.conj(Q3_2))

    # four-particle correlation
    N4_weight = dN1*(dN1 - 1)*dN2*(dN2 - 1)
    Q1Q2_N4 = (np.real(Q1_1*np.conj(Q1_2)*Q2_1*np.conj(Q2_2))
               - np.real(Q3_1*np.conj(Q1_2)*np.conj(Q2_2))
               - np.real(Q1_1*Q2_1*np.conj(Q3_2))
               + np.real(Q3_1*np.conj(Q3_2)))

    # three-particle correlation
    N3_weight = dN1*dN2*(dN2 - 1) + dN1*(dN1 - 1)*dN2
    chi_312_num = (Q3_1*np.conj(Q1_2)*np.conj(Q2_2) - Q3_1*np.conj(Q3_2)
                   + Q3_2*np.conj(Q1_1)*np.conj(Q2_1) - Q3_2*np.conj(Q3_1))

    chi_312_JK = np.zeros(nev)
    v312_JK = np.zeros(nev)
    rho312_JK = np.zeros(nev)
    v3L_JK = np.zeros(nev)
    for iev in range(nev):
        array_idx = [True]*nev
        array_idx[iev] = False
        array_idx = np.array(array_idx)

        num_JK = (np.real(np.mean(chi_312_num[array_idx]))
                  / np.mean(N3_weight[array_idx]))
        den_JK = (np.real(np.mean(Q1Q2_N4[array_idx]))
                  / np.mean(N4_weight[array_idx]))

        v3_Psi3 = np.nan_to_num(np.sqrt(np.mean(Q3_N2[array_idx])
                                        / np.mean(N2_weight[array_idx])))

        chi_312_JK[iev] = num_JK/den_JK
        v312_JK[iev] = np.nan_to_num(num_JK/np.sqrt(den_JK))
        rho312_JK[iev] = v312_JK[iev]/v3_Psi3
        v3L_JK[iev] = np.nan_to_num(np.sqrt(v3_Psi3**2 - v312_JK[iev]**2.))

    chi_312_mean = np.mean(chi_312_JK)
    chi_312_err = np.sqrt((nev - 1.)/nev
                          * np.sum((chi_312_JK - chi_312_mean)**2.))
    v312_mean = np.mean(v312_JK)
    v312_err = np.sqrt((nev - 1.)/nev*np.sum((v312_JK - v312_mean)**2.))
    rho312_mean = np.mean(rho312_JK)
    rho312_err = np.sqrt((nev - 1.)/nev*np.sum((rho312_JK - rho312_mean)**2.))
    v3L_mean = np.mean(v3L_JK)
    v3L_err = np.sqrt((nev - 1.)/nev*np.sum((v3L_JK - v3L_mean)**2.))
    results = [v3L_mean, v3L_err, v312_mean, v312_err, rho312_mean, rho312_err,
               chi_312_mean, chi_312_err]

    dN_mean = np.real(np.mean(vn_data_array1[:, 0] + vn_data_array2[:, 0]))
    dN_err = (np.std(vn_data_array1[:, 0] + vn_data_array2[:, 0])
              / np.sqrt(nev))
    if path.isfile(outputFileName):
        f = open(outputFileName, 'a')
    else:
        f = open(outputFileName, 'w')
        f.write("# cen  Nch  v3L  v_312  rho_312  chi312\n")
    f.write("{:.3f}  {:.5e}  {:.5e}".format(cenLabel, dN_mean, dN_err))
    for ires in results:
        f.write("  {:.5e}".format(ires))
    f.write("\n")
    f.close()


def calculateNonLinearResponseV4_2sub(vn_data_array1, vn_data_array2,
                                      outputFileName, cenLabel):
    """
        this function computes the non-linear response coefficients

        chi_422 = Re<V_4*conj(V_2)^2>/(<|V_2|^4>)
        v_422 = Re<V_4*conj(V_2)^2>/sqrt(<|V_2|^4>)
        rho_422 = Re<V_4*conj(V_2)^2>/sqrt(<|V_4|^2*<|V_2|^4>)
        v4_L = sqrt(|V_4|^2 - Re<V_4*conj(V_2)^2>)

        we use one flow vector for V_n and the other for conj(V_n)
        we will use Jackknife resampling method to estimate
        the statistical error
    """
    nev = len(vn_data_array1[:, 0])
    dN1 = np.real(vn_data_array1[:, -1])
    dN2 = np.real(vn_data_array2[:, -1])

    Q2_1 = dN1*vn_data_array1[:, 3]
    Q4_1 = dN1*vn_data_array1[:, 5]
    Q2_2 = dN2*vn_data_array2[:, 3]
    Q4_2 = dN2*vn_data_array2[:, 5]

    # two-particle correlation
    N2_weight = dN1*dN2
    Q4_N2 = np.real(Q4_1*np.conj(Q4_2))

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

    chi_422_JK = np.zeros(nev)
    v422_JK = np.zeros(nev)
    rho422_JK = np.zeros(nev)
    v4L_JK = np.zeros(nev)
    for iev in range(nev):
        array_idx = [True]*nev
        array_idx[iev] = False
        array_idx = np.array(array_idx)

        num_JK = (np.real(np.mean(chi_422_num[array_idx]))
                  / np.mean(N3_weight[array_idx]))
        den_JK = (np.real(np.mean(Q2_N4[array_idx]))
                  / np.mean(N4_weight[array_idx]))

        v4_Psi4 = np.nan_to_num(np.sqrt(np.mean(Q4_N2[array_idx])
                                        / np.mean(N2_weight[array_idx])))

        chi_422_JK[iev] = num_JK/den_JK
        v422_JK[iev] = np.nan_to_num(num_JK/np.sqrt(den_JK))
        rho422_JK[iev] = v422_JK[iev]/v4_Psi4
        v4L_JK[iev] = np.nan_to_num(np.sqrt(v4_Psi4**2 - v422_JK[iev]**2.))

    chi_422_mean = np.mean(chi_422_JK)
    chi_422_err = np.sqrt((nev - 1.)/nev
                          * np.sum((chi_422_JK - chi_422_mean)**2.))
    v422_mean = np.mean(v422_JK)
    v422_err = np.sqrt((nev - 1.)/nev*np.sum((v422_JK - v422_mean)**2.))
    rho422_mean = np.mean(rho422_JK)
    rho422_err = np.sqrt((nev - 1.)/nev*np.sum((rho422_JK - rho422_mean)**2.))
    v4L_mean = np.mean(v4L_JK)
    v4L_err = np.sqrt((nev - 1.)/nev*np.sum((v4L_JK - v4L_mean)**2.))
    results = [v4L_mean, v4L_err, v422_mean, v422_err, rho422_mean, rho422_err,
               chi_422_mean, chi_422_err]

    dN_mean = np.real(np.mean(vn_data_array1[:, 0] + vn_data_array2[:, 0]))
    dN_err = (np.std(vn_data_array1[:, 0] + vn_data_array2[:, 0])
              / np.sqrt(nev))
    if path.isfile(outputFileName):
        f = open(outputFileName, 'a')
    else:
        f = open(outputFileName, 'w')
        f.write("# cen  Nch  v4L  v_422  rho_422  chi422\n")
    f.write("{:.3f}  {:.5e}  {:.5e}".format(cenLabel, dN_mean, dN_err))
    for ires in results:
        f.write("  {:.5e}".format(ires))
    f.write("\n")
    f.close()


def calculateNonLinearResponseV5_2sub(vn_data_array1, vn_data_array2,
                                      outputFileName, cenLabel):
    """
        this function computes the non-linear response coefficients

        chi_523 = Re<V_5*conj(V_2*V_3)>/(<|V_2|^2*|V_3|^2>)
        v_523 = Re<V_5*conj(V_2*V_3)>/sqrt(<|V_2|^2*|V_3|^2>)
        rho_523 = Re<V_5*conj(V_2*V_3)>/sqrt(<|V_5|^2><|V_2|^2*|V_3|^2>)
        v5_L = sqrt(|V_5|^2 - Re<V_5*conj(V_2*V_3)>)

        we use one flow vector for V_n and the other for conj(V_n)
        we will use Jackknife resampling method to estimate
        the statistical error
    """
    nev = len(vn_data_array1[:, 0])
    dN1 = np.real(vn_data_array1[:, -1])
    dN2 = np.real(vn_data_array2[:, -1])

    Q2_1 = dN1*vn_data_array1[:, 3]
    Q3_1 = dN1*vn_data_array1[:, 4]
    Q5_1 = dN1*vn_data_array1[:, 6]
    Q2_2 = dN2*vn_data_array2[:, 3]
    Q3_2 = dN2*vn_data_array2[:, 4]
    Q5_2 = dN2*vn_data_array2[:, 6]

    # two-particle correlation
    N2_weight = dN1*dN2
    Q5_N2 = np.real(Q5_1*np.conj(Q5_2))

    # four-particle correlation
    N4_weight = dN1*(dN1 - 1)*dN2*(dN2 - 1)
    Q2Q3_N4 = (np.real(Q2_1*np.conj(Q2_2)*Q3_1*np.conj(Q3_2))
               - np.real(Q5_1*np.conj(Q2_2)*np.conj(Q3_2))
               - np.real(Q2_1*Q3_1*np.conj(Q5_2))
               + np.real(Q5_1*np.conj(Q5_2)))

    # three-particle correlation
    N3_weight = dN1*dN2*(dN2 - 1) + dN1*(dN1 - 1)*dN2
    chi_523_num = (Q5_1*np.conj(Q2_2)*np.conj(Q3_2) - Q5_1*np.conj(Q5_2)
                   + Q5_2*np.conj(Q2_1)*np.conj(Q3_1) - Q5_2*np.conj(Q5_1))

    chi_523_JK = np.zeros(nev)
    v523_JK = np.zeros(nev)
    rho523_JK = np.zeros(nev)
    v5L_JK = np.zeros(nev)
    for iev in range(nev):
        array_idx = [True]*nev
        array_idx[iev] = False
        array_idx = np.array(array_idx)

        num_JK = (np.real(np.mean(chi_523_num[array_idx]))
                  / np.mean(N3_weight[array_idx]))
        den_JK = (np.real(np.mean(Q2Q3_N4[array_idx]))
                  / np.mean(N4_weight[array_idx]))

        v5_Psi5 = np.nan_to_num(np.sqrt(np.mean(Q5_N2[array_idx])
                                        / np.mean(N2_weight[array_idx])))

        chi_523_JK[iev] = num_JK/den_JK
        v523_JK[iev] = np.nan_to_num(num_JK/np.sqrt(den_JK))
        rho523_JK[iev] = v523_JK[iev]/v5_Psi5
        v5L_JK[iev] = np.nan_to_num(np.sqrt(v5_Psi5**2 - v523_JK[iev]**2.))

    chi_523_mean = np.mean(chi_523_JK)
    chi_523_err = np.sqrt((nev - 1.)/nev
                          * np.sum((chi_523_JK - chi_523_mean)**2.))
    v523_mean = np.mean(v523_JK)
    v523_err = np.sqrt((nev - 1.)/nev*np.sum((v523_JK - v523_mean)**2.))
    rho523_mean = np.mean(rho523_JK)
    rho523_err = np.sqrt((nev - 1.)/nev*np.sum((rho523_JK - rho523_mean)**2.))
    v5L_mean = np.mean(v5L_JK)
    v5L_err = np.sqrt((nev - 1.)/nev*np.sum((v5L_JK - v5L_mean)**2.))
    results = [v5L_mean, v5L_err, v523_mean, v523_err, rho523_mean, rho523_err,
               chi_523_mean, chi_523_err]

    dN_mean = np.real(np.mean(vn_data_array1[:, 0] + vn_data_array2[:, 0]))
    dN_err = (np.std(vn_data_array1[:, 0] + vn_data_array2[:, 0])
              / np.sqrt(nev))
    if path.isfile(outputFileName):
        f = open(outputFileName, 'a')
    else:
        f = open(outputFileName, 'w')
        f.write("# cen  Nch  v5L  v_523  rho_523  chi523\n")
    f.write("{:.3f}  {:.5e}  {:.5e}".format(cenLabel, dN_mean, dN_err))
    for ires in results:
        f.write("  {:.5e}".format(ires))
    f.write("\n")
    f.close()


def calculate_vn4_vn6(vn_data_array, outputFileName_vn4,
                      outputFileName42, outputFileName64, cenLabel):
    """
        this funciton computes the 4 particle cumulant vn{4}
            vn{4} = (2 <v_n*conj(v_n)>**2 - <(v_n*conj(v_n))**2.>)**(1/4)

        this funciton also computes the ratio of
        the 4-particle cumulant vn{4} over the 2-particle cumulant vn{2}
        and Fn = sqrt((vn{2}^2 - vn{4}^2)/(vn{2}^2 + vn{4}^2))

            vn{4} = (2 <v_n*conj(v_n)>**2 - <(v_n*conj(v_n))**2.>)**(1/4)
            vn{2} = (<v_n*conj(v_n)>)**(1/2)

        this funciton also computes the ratio of
        the 6-particle cumulant vn{6} over the 4-particle cumulant vn{4}
            cn{6} = <<6>> - 9<<2>><<4>> + 12<<2>>^3
            vn{6} = (cn{6}/4)**(1/6)
            vn{4} = (2 <v_n*conj(v_n)>**2 - <(v_n*conj(v_n))**2.>)**(1/4)
        and compute skewness estimator gamma_1
            gamma_1 = -6\\sqrt{2}*vn{4}^2*(vn{4} - vn{6})
                                         /(vn{2}^2 - vn{4}^2)^(3/2)

        we will use Jackknife resampling method to estimate
        the statistical error
    """
    nev = len(vn_data_array[:, 0])
    dN = np.real(vn_data_array[:, -1])
    Q2 = dN*vn_data_array[:, 3]
    Q3 = dN*vn_data_array[:, 4]
    Q4 = dN*vn_data_array[:, 5]
    Q6 = dN*vn_data_array[:, 7]

    # two-particle correlation
    N2_weight = dN*(dN - 1.)
    Q2_N2 = np.real(Q2*np.conj(Q2)) - dN
    Q3_N2 = np.real(Q3*np.conj(Q3)) - dN

    # four-particle correlation
    N4_weight = dN*(dN - 1.)*(dN - 2.)*(dN - 3.)
    Q2_4 = ((np.abs(Q2)**4.) - 2.*np.real(Q4*np.conj(Q2)*np.conj(Q2))
            - 4.*(dN - 2.)*(np.abs(Q2)**2.) + np.abs(Q4)**2.
            + 2*dN*(dN - 3.))
    Q3_4 = ((np.abs(Q3)**4.) - 2.*np.real(Q6*np.conj(Q3)*np.conj(Q3))
            - 4.*(dN - 2.)*(np.abs(Q3)**2.) + np.abs(Q6)**2.
            + 2*dN*(dN - 3.))

    # six-particle correlation
    N6_weight = dN*(dN - 1.)*(dN - 2.)*(dN - 3.)*(dN - 4.)*(dN - 5.)
    Q2_6 = (np.abs(Q2)**6. + 9*(np.abs(Q4)**2.)*(np.abs(Q2)**2.)
            - 6.*np.real(Q4*Q2*np.conj(Q2)*np.conj(Q2)*np.conj(Q2))
            + 4.*np.real(Q6*np.conj(Q2)*np.conj(Q2)*np.conj(Q2))
            - 12.*np.real(Q6*np.conj(Q4)*np.conj(Q2))
            + 18.*(dN - 4.)*np.real(Q4*np.conj(Q2)*np.conj(Q2))
            + 4.*(np.abs(Q6)**2.)
            - 9.*(dN - 4.)*((np.abs(Q2)**4.) + (np.abs(Q4)**2.))
            + 18.*(dN - 5.)*(dN - 2.)*(np.abs(Q2)**2.)
            - 6.*dN*(dN - 4.)*(dN - 5.))

    # calcualte observables with Jackknife resampling method
    C2_4_array = np.zeros(nev)
    C3_4_array = np.zeros(nev)
    r2_array = np.zeros(nev)
    r3_array = np.zeros(nev)
    F2_array = np.zeros(nev)
    F3_array = np.zeros(nev)
    r26_array = np.zeros(nev)
    gamma1_array = np.zeros(nev)
    for iev in range(nev):
        array_idx = [True]*nev
        array_idx[iev] = False
        array_idx = np.array(array_idx)

        # C_n{4}
        C_2_2 = np.mean(Q2_N2[array_idx])/np.mean(N2_weight[array_idx])
        C24_tmp = np.mean(Q2_4[array_idx])/np.mean(N4_weight[array_idx])
        C_2_4 = C24_tmp - 2.*(C_2_2**2.)
        C_2_6 = (np.mean(Q2_6[array_idx])/np.mean(N6_weight[array_idx])
                 - 9.*C_2_2*C24_tmp + 12.*(C_2_2**3.))
        C2_4_array[iev] = C_2_4
        if C_2_4 < 0. and C_2_2 > 0.:
            v2_4 = (-C_2_4)**0.25
            v2_2 = np.sqrt(C_2_2)
            r2_array[iev] = v2_4/v2_2
            F2_array[iev] = np.nan_to_num(
                    np.sqrt((v2_2**2. - v2_4**2.)
                            / (v2_2**2. + v2_4**2. + 1e-15)))
            if C_2_6 > 0.:
                v2_6 = (C_2_6/4.)**(1./6.)
                r26_array[iev] = v2_6/v2_4
                gamma1_array[iev] = (-6.*np.sqrt(2)*(v2_4**2.)*(v2_4 - v2_6)
                                     / (v2_2**2. - v2_4**2.)**(1.5))

        C_3_2 = np.mean(Q3_N2[array_idx])/np.mean(N2_weight[array_idx])
        C34_tmp = np.mean(Q3_4[array_idx])/np.mean(N4_weight[array_idx])
        C_3_4 = C34_tmp - 2.*C_3_2**2.
        C3_4_array[iev] = C_3_4
        if C_3_4 < 0. and C_3_2 > 0.:
            v3_4 = (-C_3_4)**0.25
            v3_2 = np.sqrt(C_3_2)
            r3_array[iev] = v3_4/v3_2
            F3_array[iev] = np.sqrt((v3_2**2. - v3_4**2.)
                                    / (v3_2**2. + v3_4**2. + 1e-15))

    dN_mean = np.real(np.mean(vn_data_array[:, 0]))
    dN_err = np.std(vn_data_array[:, 0])/np.sqrt(nev)
    C2_4_mean = np.mean(C2_4_array)
    C2_4_err = np.sqrt((nev - 1.)/nev*np.sum((C2_4_array - C2_4_mean)**2.))
    C3_4_mean = np.mean(C3_4_array)
    C3_4_err = np.sqrt((nev - 1.)/nev*np.sum((C3_4_array - C3_4_mean)**2.))

    v2_4 = 0.0
    v2_4_err = 0.0
    if C2_4_mean < 0:
        v2_4 = (-C2_4_mean)**0.25
        v2_4_err = 0.25*((-C2_4_mean)**(-0.75))*C2_4_err

    v3_4 = 0.0
    v3_4_err = 0.0
    if C3_4_mean < 0:
        v3_4 = (-C3_4_mean)**0.25
        v3_4_err = 0.25*((-C3_4_mean)**(-0.75))*C3_4_err
    results = [v2_4, v2_4_err, C2_4_mean, C2_4_err,
               v3_4, v3_4_err, C3_4_mean, C3_4_err,]
    if path.isfile(outputFileName_vn4):
        f = open(outputFileName_vn4, 'a')
    else:
        f = open(outputFileName_vn4, 'w')
        f.write("# cen  Nch  vn{4}  vn{4}_err  Cn{4}  Cn{4}_err (n=2-3)\n")
    f.write("{:.3f}  {:.5e}  {:.5e}".format(cenLabel, dN_mean, dN_err))
    for ires in results:
        f.write("  {:.5e}".format(ires))
    f.write("\n")
    f.close()

    # now the ratios
    r2_mean = np.mean(r2_array)
    r2_err = np.sqrt((nev - 1.)/nev*np.sum((r2_array - r2_mean)**2.))
    r3_mean = np.mean(r3_array)
    r3_err = np.sqrt((nev - 1.)/nev*np.sum((r3_array - r3_mean)**2.))

    F2_mean = np.mean(F2_array)
    F2_err = np.sqrt((nev - 1.)/nev*np.sum((F2_array - F2_mean)**2.))
    F3_mean = np.mean(F3_array)
    F3_err = np.sqrt((nev - 1.)/nev*np.sum((F3_array - F3_mean)**2.))
    results = [r2_mean, r2_err, F2_mean, F2_err,
               r3_mean, r3_err, F3_mean, F3_err]
    if path.isfile(outputFileName42):
        f = open(outputFileName42, 'a')
    else:
        f = open(outputFileName42, 'w')
        f.write("# cen  Nch  vn{4}/vn{2}  (vn{4}/vn{2})_err  Fn  Fn_err\n")
        f.write("# Fn = sqrt((vn{2}^2 - vn{4}^2)/(vn{2}^2 + vn{4}^2))\n")
    f.write("{:.3f}  {:.5e}  {:.5e}".format(cenLabel, dN_mean, dN_err))
    for ires in results:
        f.write("  {:.5e}".format(ires))
    f.write("\n")
    f.close()

    r26_mean = np.mean(r26_array)
    r26_err = np.sqrt((nev - 1.)/nev*np.sum((r26_array - r26_mean)**2.))
    gamma1_mean = np.nan_to_num(np.mean(gamma1_array))
    gamma1_err = np.nan_to_num(np.sqrt((nev - 1.)/nev*np.sum(
                                        (gamma1_array - gamma1_mean)**2.)))
    if path.isfile(outputFileName64):
        f = open(outputFileName64, 'a')
    else:
        f = open(outputFileName64, 'w')
        f.write("# cen  Nch  vn{6}/vn{4}  (vn{6}/vn{4})_err  "
                + "gamma_1  gamma_1_err\n")
    f.write("{:.3f}  {:.5e}  {:.5e}  {:.5e}  {:.5e}  {:.5e}  {:.5e}\n".format(
        cenLabel, dN_mean, dN_err, r26_mean, r26_err, gamma1_mean, gamma1_err))
    f.close()
    return


def calculate_vn4_2sub(vn_data_array1, vn_data_array2,
                       outputFileName_vn4, outputFileName42, cenLabel):
    """
        this funciton computes the 4 particle cumulant vn{4}
            vn{4} = (2 <v_n*conj(v_n)>**2 - <(v_n*conj(v_n))**2.>)**(1/4)

        this funciton also computes the ratio of
        the 4-particle cumulant vn{4} over the 2-particle cumulant vn{2}
        and Fn = sqrt((vn{2}^2 - vn{4}^2)/(vn{2}^2 + vn{4}^2))

            vn{4} = (2 <v_n*conj(v_n)>**2 - <(v_n*conj(v_n))**2.>)**(1/4)
            vn{2} = (<v_n*conj(v_n)>)**(1/2)

        we will use Jackknife resampling method to estimate
        the statistical error
    """
    nev = len(vn_data_array1[:, 0])
    dN1 = np.real(vn_data_array1[:, -1])
    dN2 = np.real(vn_data_array2[:, -1])
    Q2_1 = dN1*vn_data_array1[:, 3]
    Q2_2 = dN2*vn_data_array2[:, 3]
    Q3_1 = dN1*vn_data_array1[:, 4]
    Q3_2 = dN2*vn_data_array2[:, 4]
    Q4_1 = dN1*vn_data_array1[:, 5]
    Q4_2 = dN2*vn_data_array2[:, 5]
    Q6_1 = dN1*vn_data_array1[:, 7]
    Q6_2 = dN2*vn_data_array2[:, 7]

    # two-particle correlation
    N2_weight = dN1*dN2
    Q2_N2 = np.real(Q2_1*np.conj(Q2_2))
    Q3_N2 = np.real(Q3_1*np.conj(Q3_2))

    # four-particle correlation
    N4_weight = dN1*(dN1 - 1)*dN2*(dN2 - 1)
    Q2_N4 = (np.real(Q2_1*np.conj(Q2_2)*Q2_1*np.conj(Q2_2))
             - np.real(Q4_1*np.conj(Q2_2)*np.conj(Q2_2))
             - np.real(Q2_1*Q2_1*np.conj(Q4_2))
             + np.real(Q4_1*np.conj(Q4_2)))

    Q3_N4 = (np.real(Q3_1*np.conj(Q3_2)*Q3_1*np.conj(Q3_2))
             - np.real(Q6_1*np.conj(Q3_2)*np.conj(Q3_2))
             - np.real(Q3_1*Q3_1*np.conj(Q6_2))
             + np.real(Q6_1*np.conj(Q6_2)))

    # calcualte observables with Jackknife resampling method
    C2_4_array = np.zeros(nev)
    C3_4_array = np.zeros(nev)
    r2_array = np.zeros(nev)
    r3_array = np.zeros(nev)
    F2_array = np.zeros(nev)
    F3_array = np.zeros(nev)
    for iev in range(nev):
        array_idx = [True]*nev
        array_idx[iev] = False
        array_idx = np.array(array_idx)

        # C_n{4}
        C_2_2 = np.mean(Q2_N2[array_idx])/np.mean(N2_weight[array_idx])
        C24_tmp = np.mean(Q2_N4[array_idx])/np.mean(N4_weight[array_idx])
        C_2_4 = C24_tmp - 2.*C_2_2**2.
        C2_4_array[iev] = C_2_4
        if C_2_4 < 0. and C_2_2 > 0.:
            v2_4 = (-C_2_4)**0.25
            v2_2 = np.sqrt(C_2_2)
            r2_array[iev] = v2_4/v2_2
            F2_array[iev] = np.sqrt((v2_2**2. - v2_4**2.)
                                    / (v2_2**2. + v2_4**2. + 1e-15))

        C_3_2 = np.mean(Q3_N2[array_idx])/np.mean(N2_weight[array_idx])
        C34_tmp = np.mean(Q3_N4[array_idx])/np.mean(N4_weight[array_idx])
        C_3_4 = C34_tmp - 2.*C_3_2**2.
        C3_4_array[iev] = C_3_4
        if C_3_4 < 0. and C_3_2 > 0.:
            v3_4 = (-C_3_4)**0.25
            v3_2 = np.sqrt(C_3_2)
            r3_array[iev] = v3_4/v3_2
            F3_array[iev] = np.sqrt((v3_2**2. - v3_4**2.)
                                    / (v3_2**2. + v3_4**2. + 1e-15))

    C2_4_mean = np.mean(C2_4_array)
    C2_4_err = np.sqrt((nev - 1.)/nev*np.sum((C2_4_array - C2_4_mean)**2.))
    C3_4_mean = np.mean(C3_4_array)
    C3_4_err = np.sqrt((nev - 1.)/nev*np.sum((C3_4_array - C3_4_mean)**2.))

    v2_4 = 0.0
    v2_4_err = 0.0
    if C2_4_mean < 0:
        v2_4 = (-C2_4_mean)**0.25
        v2_4_err = 0.25*((-C2_4_mean)**(-0.75))*C2_4_err

    v3_4 = 0.0
    v3_4_err = 0.0
    if C3_4_mean < 0:
        v3_4 = (-C3_4_mean)**0.25
        v3_4_err = 0.25*((-C3_4_mean)**(-0.75))*C3_4_err
    results = [v2_4, v2_4_err, C2_4_mean, C2_4_err,
               v3_4, v3_4_err, C3_4_mean, C3_4_err,]

    dN_mean = np.real(np.mean(vn_data_array1[:, 0] + vn_data_array2[:, 0]))
    dN_err = (np.std(vn_data_array1[:, 0] + vn_data_array2[:, 0])
              / np.sqrt(nev))
    if path.isfile(outputFileName_vn4):
        f = open(outputFileName_vn4, 'a')
    else:
        f = open(outputFileName_vn4, 'w')
        f.write("# cen  Nch  vn{4}  vn{4}_err  Cn{4}  Cn{4}_err (n=2-3)\n")
    f.write("{:.3f}  {:.5e}  {:.5e}".format(cenLabel, dN_mean, dN_err))
    for ires in results:
        f.write("  {:.5e}".format(ires))
    f.write("\n")
    f.close()

    # now the ratios
    r2_mean = np.mean(r2_array)
    r2_err = np.sqrt((nev - 1.)/nev*np.sum((r2_array - r2_mean)**2.))
    r3_mean = np.mean(r3_array)
    r3_err = np.sqrt((nev - 1.)/nev*np.sum((r3_array - r3_mean)**2.))

    F2_mean = np.mean(F2_array)
    F2_err = np.sqrt((nev - 1.)/nev*np.sum((F2_array - F2_mean)**2.))
    F3_mean = np.mean(F3_array)
    F3_err = np.sqrt((nev - 1.)/nev*np.sum((F3_array - F3_mean)**2.))

    results = [r2_mean, r2_err, F2_mean, F2_err,
               r3_mean, r3_err, F3_mean, F3_err]
    if path.isfile(outputFileName42):
        f = open(outputFileName42, 'a')
    else:
        f = open(outputFileName42, 'w')
        f.write("# cen  Nch  vn{4}/vn{2}  (vn{4}/vn{2})_err  Fn  Fn_err\n")
        f.write("# Fn = sqrt((vn{2}^2 - vn{4}^2)/(vn{2}^2 + vn{4}^2))\n")
    f.write("{:.3f}  {:.5e}  {:.5e}".format(cenLabel, dN_mean, dN_err))
    for ires in results:
        f.write("  {:.5e}".format(ires))
    f.write("\n")
    f.close()
    return


def calcualte_vn_2_with_gap(vn_data_array_sub1, vn_data_array_sub2,
                            outputFileName, cenLabel):
    """
        this function computes vn{2} and its stat. err.
        using two subevents with a eta gap
    """
    nev = len(vn_data_array_sub1[:, 0])
    dN1 = np.real(vn_data_array_sub1[:, -1])
    dN1 = dN1.reshape(len(dN1), 1)
    dN2 = np.real(vn_data_array_sub2[:, -1])
    dN2 = dN1.reshape(len(dN2), 1)
    Qn_array1 = dN1*vn_data_array_sub1[:, 2:-1]
    Qn_array2 = dN2*vn_data_array_sub2[:, 2:-1]

    corr = (Qn_array1*np.conj(Qn_array2))/(dN1*dN2)
    vn_2 = np.nan_to_num(np.sqrt(np.real(np.mean(corr, 0))))
    vn_2_err = np.nan_to_num(np.std(np.real(corr), 0)/np.sqrt(nev)/2./vn_2)
    dN_mean = np.real(np.mean(vn_data_array_sub1[:, 0]
                              + vn_data_array_sub2[:, 0]))
    dN_err = (np.std(vn_data_array_sub1[:, 0] + vn_data_array_sub2[:, 0])
              / np.sqrt(nev))
    if path.isfile(outputFileName):
        f = open(outputFileName, 'a')
    else:
        f = open(outputFileName, 'w')
        f.write("# cen  Nch  vn{2}  vn{2}_err (n = 1-9)\n")
    f.write("{:.3f}  {:.5e}  {:.5e}".format(cenLabel, dN_mean, dN_err))
    for i in range(len(vn_2)):
        f.write("  {:.5e}  {:.5e}".format(vn_2[i], vn_2_err[i]))
    f.write("\n")
    f.close()
    return


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
    piddNArr = []
    pidmeanpTArr = []
    Ngluons = []
    for event_name in selected_events_list:
        Ngluons.append(data[event_name]['NgluonEst'])
        QnArr1.append(data[event_name]['ALICE_eta_-0p4_0p4'])
        QnArr2.append(data[event_name]['ALICE_eta_-0p8_-0p4'])
        QnArr3.append(data[event_name]['ALICE_eta_0p4_0p8'])
        dNtmp = []
        meanpTtmp = []
        dNtmp.append(data[event_name]['Nch'])
        meanpTtmp.append(data[event_name]['mean_pT_ch'])
        for pidName in pidList:
            tmp = data[event_name][f'{pidName}_dNdy_meanpT']
            dNtmp.append(tmp[0])
            meanpTtmp.append(tmp[1])
        piddNArr.append(dNtmp)
        pidmeanpTArr.append(meanpTtmp)
    Ngluons = np.array(Ngluons)
    print(f"Ngluons = {Ngluons.min():.2f} - {Ngluons.max():.2f}")
    QnArr1 = np.array(QnArr1)
    QnArr2 = np.array(QnArr2)
    QnArr3 = np.array(QnArr3)
    piddNArr = np.array(piddNArr)
    pidmeanpTArr = np.array(pidmeanpTArr)

    calcualte_vn_2_with_gap(QnArr2, QnArr3, "vn2_sub.dat", cenLabel)
    calculate_vn4_2sub(QnArr2, QnArr3, "vn4_2sub.dat",
                       "vn4_2sub_over_vn2.dat", cenLabel)
    calculate_vn4_vn6(QnArr1, "vn4.dat", "vn4_over_vn2.dat",
                      "vn6_over_vn4.dat", cenLabel)
    calculate_pid_dN(piddNArr, "pid_dN.dat", cenLabel)
    calculate_pid_meanpT(pidmeanpTArr, "pid_meanpT.dat", cenLabel)
    calculateNonLinearResponseV2_2sub(QnArr2, QnArr3,
                                      "nonLinearV2_2sub.dat", cenLabel)
    calculateNonLinearResponseV3_2sub(QnArr2, QnArr3,
                                      "nonLinearV3_2sub.dat", cenLabel)
    calculateNonLinearResponseV4_2sub(QnArr2, QnArr3,
                                      "nonLinearV4_2sub.dat", cenLabel)
    calculateNonLinearResponseV5_2sub(QnArr2, QnArr3,
                                      "nonLinearV5_2sub.dat", cenLabel)
    calculateSymmetricCumulant2sub(QnArr2, QnArr3,
                                   "symmetricCumulants_2sub.dat", cenLabel)
