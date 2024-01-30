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

    # four-particle correlation
    N4_weight = dN1*(dN1 - 1)*dN2*(dN2 - 1)
    Q2_N4 = (  np.real(Q2_1*np.conj(Q2_2)*Q2_1*np.conj(Q2_2))
             - np.real(Q4_1*np.conj(Q2_2)*np.conj(Q2_2))
             - np.real(Q2_1*Q2_1*np.conj(Q4_2))
             + np.real(Q4_1*np.conj(Q4_2)))

    # three-particle correlation
    N3_weight = dN1*dN2*(dN2 - 1) + dN1*(dN1 - 1)*dN2
    chi_422_num = (  Q4_1*np.conj(Q2_2)*np.conj(Q2_2) - Q4_1*np.conj(Q4_2)
                   + Q4_2*np.conj(Q2_1)*np.conj(Q2_1) - Q4_2*np.conj(Q4_1))

    num = np.real(np.mean(chi_422_num))/np.mean(N3_weight)
    den = np.mean(Q2_N4)/np.mean(N4_weight)
    chi_422_mean = num/den
    return chi_422_mean


def calculate_chi_523(vn_data_array1, vn_data_array2):
    """
        this function computes the non-linear response coefficients

        chi_523 = Re<V_5*conj(V_2*V_3)>/(<|V_2|^2*|V_3|^2>)

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

    # four-particle correlation
    N4_weight = dN1*(dN1 - 1)*dN2*(dN2 - 1)
    Q2Q3_N4 = (np.real(Q2_1*np.conj(Q2_2)*Q3_1*np.conj(Q3_2))
               - np.real(Q5_1*np.conj(Q2_2*Q3_2))
               - np.real(Q2_1*Q3_1*np.conj(Q5_2))
               + np.real(Q5_1*np.conj(Q5_2)))

    # three-particle correlation
    N3_weight = dN1*dN2*(dN2 - 1) + dN1*(dN1 - 1)*dN2
    chi_523_num = (Q5_1*np.conj(Q2_2*Q3_2) - Q5_1*np.conj(Q5_2)
                   + Q5_2*np.conj(Q2_1*Q3_1) - Q5_2*np.conj(Q5_1))

    num = np.real(np.mean(chi_523_num))/np.mean(N3_weight)
    den = np.mean(Q2Q3_N4)/np.mean(N4_weight)
    chi_523_mean = num/den
    return(chi_523_mean)


def calculateNonLinearResponseV67_2sub(vn_data_array1, vn_data_array2,
                                       V4L_1, V4L_2, V5L_1, V5L_2):
    """
        This function computes the linear components
        for V6 = V6L + chi_6222 V2^3 + chi_633 V3^2 + chi_624 V2*V4L
        and V7 = V7L + chi_7223 V2^2*V3 + chi_734 V3*V4L + chi_725 V2*V5L
    """
    nev = len(vn_data_array1[:, 0])

    V2_1 = vn_data_array1[:, 3]
    V3_1 = vn_data_array1[:, 4]
    V4_1 = vn_data_array1[:, 5]
    V5_1 = vn_data_array1[:, 6]
    V6_1 = vn_data_array1[:, 7]
    V7_1 = vn_data_array1[:, 8]

    V2_2 = vn_data_array2[:, 3]
    V3_2 = vn_data_array2[:, 4]
    V4_2 = vn_data_array2[:, 5]
    V5_2 = vn_data_array2[:, 6]
    V6_2 = vn_data_array2[:, 7]
    V7_2 = vn_data_array2[:, 8]

    chi_6222_num = V6_1*(np.conj(V2_2)**3) + V6_2*(np.conj(V2_1)**3)
    V2_6 = 2*np.real((V2_1*np.conj(V2_2))**3)
    V32_V23 = (V3_1**2)*(np.conj(V2_2)**3) + (V3_2**2)*(np.conj(V2_1)**3)
    V2V4L_V23 = (  (V2_1*V4L_1)*(np.conj(V2_2)**3)
                 + (V2_2*V4L_2)*(np.conj(V2_1)**3))

    chi_633_num = V6_1*(np.conj(V3_2)**2) + V6_2*(np.conj(V3_1)**2)
    V3_4 = 2*np.real((V3_1*np.conj(V3_2))**2)
    V2V4L_V32 = (  (V2_1*V4L_1)*(np.conj(V3_2)**2)
                 + (V2_2*V4L_2)*(np.conj(V3_1)**2))

    v624_num = V6_1*(np.conj(V2_2*V4L_2)) + V6_2*(np.conj(V2_1*V4L_1))
    V2V4L_2 = 2*np.real(V2_1*np.conj(V2_2)*V4L_1*np.conj(V4L_2))

    num_JK1 = np.mean(chi_6222_num)
    den_JK11 = np.mean(V2_6)
    den_JK12 = np.mean(V32_V23)
    den_JK13 = np.mean(V2V4L_V23)

    num_JK2 = np.mean(chi_633_num)
    den_JK21 = np.conj(den_JK12)
    den_JK22 = np.mean(V3_4)
    den_JK23 = np.mean(V2V4L_V32)

    num_JK3 = np.mean(v624_num)
    den_JK31 = np.conj(den_JK13)
    den_JK32 = np.conj(den_JK23)
    den_JK33 = np.mean(V2V4L_2)

    array_lhs = np.array([num_JK1, num_JK2, num_JK3], dtype=np.cfloat)

    array_rhs = np.array([[den_JK11, den_JK12, den_JK13],
                          [den_JK21, den_JK22, den_JK23],
                          [den_JK31, den_JK32, den_JK33],
                         ], dtype=np.cfloat)
    chi_6 = np.linalg.solve(array_rhs, array_lhs)
    chi_6222_mean = chi_6[0]
    chi_633_mean = chi_6[1]
    chi_624_mean = chi_6[2]

    V6L_1 = (V6_1 - chi_6222_mean*(V2_1**3) - chi_633_mean*(V3_1**2)
             - chi_624_mean*(V2_1*V4L_1))
    V6L_2 = (V6_2 - chi_6222_mean*(V2_2**3) - chi_633_mean*(V3_2**2)
             - chi_624_mean*(V2_2*V4L_2))

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

    num_JK1 = np.mean(chi_7223_num)
    den_JK11 = np.mean(V22V3_V22V3)
    den_JK12 = np.mean(V3V4_V22V3)
    den_JK13 = np.mean(V2V5_V22V3)

    num_JK2 = np.mean(chi_734_num)
    den_JK21 = np.conj(den_JK12)
    den_JK22 = np.mean(V3V4_V3V4)
    den_JK23 = np.mean(V2V5_V3V4)

    num_JK3 = np.mean(chi_725_num)
    den_JK31 = np.conj(den_JK13)
    den_JK32 = np.conj(den_JK23)
    den_JK33 = np.mean(V2V5_V2V5)

    array_lhs = np.array([num_JK1, num_JK2, num_JK3], dtype=np.cfloat)

    array_rhs = np.array([[den_JK11, den_JK12, den_JK13],
                          [den_JK21, den_JK22, den_JK23],
                          [den_JK31, den_JK32, den_JK33],
                         ], dtype=np.cfloat)
    chi_7 = np.linalg.solve(array_rhs, array_lhs)
    chi_7223_mean = chi_7[0]
    chi_734_mean = chi_7[1]
    chi_725_mean = chi_7[2]

    V7L_1 = (V7_1 - chi_7223_mean*(V2_1**2*V3_1) - chi_734_mean*(V3_1*V4L_1)
             - chi_725_mean*(V2_1*V5L_1))
    V7L_2 = (V7_2 - chi_7223_mean*(V2_2**2*V3_2) - chi_734_mean*(V3_2*V4L_2)
             - chi_725_mean*(V2_2*V5L_2))
    return(V6L_1, V6L_2, V7L_1, V7L_2)


def calculateNonLinearResponseV9_2sub(vn_data_array1, vn_data_array2,
                                      outputFileNameV9, flag, cenLabel):
    """
        This function computes the non-linear response coefficients
        for V9 = V9L
                 + chi_9333  V3^3
                 + chi_92223 V2^3*V3
                 + chi_9234  V2*V3*V4L
                 + chi_9225  V2^2*V5L
                 + chi_945   V4L*V5L
                 + chi_936   V3*V6L
                 + chi_927   V2*V7L
    """
    nev = len(vn_data_array1[:, 0])

    V2_1 = vn_data_array1[:, 3]
    V3_1 = vn_data_array1[:, 4]
    V4_1 = vn_data_array1[:, 5]
    V5_1 = vn_data_array1[:, 6]
    V6_1 = vn_data_array1[:, 7]
    V7_1 = vn_data_array1[:, 8]
    V8_1 = vn_data_array1[:, 9]
    V9_1 = vn_data_array1[:, 10]

    V2_2 = vn_data_array2[:, 3]
    V3_2 = vn_data_array2[:, 4]
    V4_2 = vn_data_array2[:, 5]
    V5_2 = vn_data_array2[:, 6]
    V6_2 = vn_data_array2[:, 7]
    V7_2 = vn_data_array2[:, 8]
    V8_2 = vn_data_array2[:, 9]
    V9_2 = vn_data_array2[:, 10]

    chi_422 = calculate_chi_422(vn_data_array1, vn_data_array2)
    V4L_1 = V4_1 - chi_422*(V2_1**2)
    V4L_2 = V4_2 - chi_422*(V2_2**2)

    chi_523 = calculate_chi_523(vn_data_array1, vn_data_array2)
    V5L_1 = V5_1 - chi_523*(V2_1*V3_1)
    V5L_2 = V5_2 - chi_523*(V2_2*V3_2)

    V6L_1, V6L_2, V7L_1, V7L_2 = calculateNonLinearResponseV67_2sub(
                vn_data_array1, vn_data_array2, V4L_1, V4L_2, V5L_1, V5L_2)

    conjV333_1 = np.conj(V3_1**3)
    conjV333_2 = np.conj(V3_2**3)
    chi_9333_num = V9_1*conjV333_2 + V9_2*conjV333_1
    V333_V333 = 2*np.real((V3_1**3)*conjV333_2)
    V2223_V333 = (V2_1**3)*V3_1*conjV333_2 + (V2_2**3)*V3_2*conjV333_1
    V234_V333 = V2_1*V3_1*V4L_1*conjV333_2 + V2_2*V3_2*V4L_2*conjV333_1
    V225_V333 = V2_1*V2_1*V5L_1*conjV333_2 + V2_2*V2_2*V5L_2*conjV333_1
    V45_V333 = V4L_1*V5L_1*conjV333_2 + V4L_2*V5L_2*conjV333_1
    V36_V333 = V3_1*V6L_1*conjV333_2 + V3_2*V6L_2*conjV333_1
    V27_V333 = V2_1*V7L_1*conjV333_2 + V2_2*V7L_2*conjV333_1

    conjV2223_1 = np.conj((V2_1**3)*V3_1)
    conjV2223_2 = np.conj((V2_2**3)*V3_2)
    chi_92223_num = V9_1*conjV2223_2 + V9_2*conjV2223_1
    V2223_V2223 = 2*np.real((V2_1**3)*V3_1*conjV2223_2)
    V234_V2223 = V2_1*V3_1*V4L_1*conjV2223_2 + V2_2*V3_2*V4L_2*conjV2223_1
    V225_V2223 = V2_1*V2_1*V5L_1*conjV2223_2 + V2_2*V2_2*V5L_2*conjV2223_1
    V45_V2223 = V4L_1*V5L_1*conjV2223_2 + V4L_2*V5L_2*conjV2223_1
    V36_V2223 = V3_1*V6L_1*conjV2223_2 + V3_2*V6L_2*conjV2223_1
    V27_V2223 = V2_1*V7L_1*conjV2223_2 + V2_2*V7L_2*conjV2223_1

    conjV234_1 = np.conj(V2_1*V3_1*V4L_1)
    conjV234_2 = np.conj(V2_2*V3_2*V4L_2)
    chi_9234_num = V9_1*conjV234_2 + V9_2*conjV234_1
    V234_V234 = 2*np.real(V2_1*V3_1*V4L_1*conjV234_2)
    V225_V234 = V2_1*V2_1*V5L_1*conjV234_2 + V2_2*V2_2*V5L_2*conjV234_1
    V45_V234 = V4L_1*V5L_1*conjV234_2 + V4L_2*V5L_2*conjV234_1
    V36_V234 = V3_1*V6L_1*conjV234_2 + V3_2*V6L_2*conjV234_1
    V27_V234 = V2_1*V7L_1*conjV234_2 + V2_2*V7L_2*conjV234_1

    conjV225_1 = np.conj(V2_1*V2_1*V5L_1)
    conjV225_2 = np.conj(V2_2*V2_2*V5L_2)
    chi_9225_num = V9_1*conjV225_2 + V9_2*conjV225_1
    V225_V225 = V2_1*V2_1*V5L_1*conjV225_2 + V2_2*V2_2*V5L_2*conjV225_1
    V45_V225 = V4L_1*V5L_1*conjV225_2 + V4L_2*V5L_2*conjV225_1
    V36_V225 = V3_1*V6L_1*conjV225_2 + V3_2*V6L_2*conjV225_1
    V27_V225 = V2_1*V7L_1*conjV225_2 + V2_2*V7L_2*conjV225_1

    conjV45_1 = np.conj(V4L_1*V5L_1)
    conjV45_2 = np.conj(V4L_2*V5L_2)
    chi_945_num = V9_1*conjV45_2 + V9_2*conjV45_1
    V45_V45 = V4L_1*V5L_1*conjV45_2 + V4L_2*V5L_2*conjV45_1
    V36_V45 = V3_1*V6L_1*conjV45_2 + V3_2*V6L_2*conjV45_1
    V27_V45 = V2_1*V7L_1*conjV45_2 + V2_2*V7L_2*conjV45_1

    conjV36_1 = np.conj(V3_1*V6L_1)
    conjV36_2 = np.conj(V3_2*V6L_2)
    chi_936_num = V9_1*conjV36_2 + V9_2*conjV36_1
    V36_V36 = V3_1*V6L_1*conjV36_2 + V3_2*V6L_2*conjV36_1
    V27_V36 = V2_1*V7L_1*conjV36_2 + V2_2*V7L_2*conjV36_1

    conjV27_1 = np.conj(V2_1*V7L_1)
    conjV27_2 = np.conj(V2_2*V7L_2)
    chi_927_num = V9_1*conjV27_2 + V9_2*conjV27_1
    V27_V27 = V2_1*V7L_1*conjV27_2 + V2_2*V7L_2*conjV27_1

    chi_9333_JK = np.zeros(nev, dtype=complex)
    chi_92223_JK = np.zeros(nev, dtype=complex)
    chi_9234_JK = np.zeros(nev, dtype=complex)
    chi_9225_JK = np.zeros(nev, dtype=complex)
    chi_945_JK = np.zeros(nev, dtype=complex)
    chi_936_JK = np.zeros(nev, dtype=complex)
    chi_927_JK = np.zeros(nev, dtype=complex)
    for iev in range(nev):
        array_idx = [True]*nev
        array_idx[iev] = False
        array_idx = np.array(array_idx)

        b1 = np.mean(chi_9333_num[array_idx])
        A11 = np.mean(V333_V333[array_idx])
        A12 = np.mean(V2223_V333[array_idx])
        A13 = np.mean(V234_V333[array_idx])
        A14 = np.mean(V225_V333[array_idx])
        A15 = np.mean(V45_V333[array_idx])
        A16 = np.mean(V36_V333[array_idx])
        A17 = np.mean(V27_V333[array_idx])

        b2 = np.mean(chi_92223_num[array_idx])
        A21 = np.conj(A12)
        A22 = np.mean(V2223_V2223[array_idx])
        A23 = np.mean(V234_V2223[array_idx])
        A24 = np.mean(V225_V2223[array_idx])
        A25 = np.mean(V45_V2223[array_idx])
        A26 = np.mean(V36_V2223[array_idx])
        A27 = np.mean(V27_V2223[array_idx])

        b3 = np.mean(chi_9234_num[array_idx])
        A31 = np.conj(A13)
        A32 = np.conj(A23)
        A33 = np.mean(V234_V234[array_idx])
        A34 = np.mean(V225_V234[array_idx])
        A35 = np.mean(V45_V234[array_idx])
        A36 = np.mean(V36_V234[array_idx])
        A37 = np.mean(V27_V234[array_idx])

        b4 = np.mean(chi_9225_num[array_idx])
        A41 = np.conj(A14)
        A42 = np.conj(A24)
        A43 = np.conj(A34)
        A44 = np.mean(V225_V225[array_idx])
        A45 = np.mean(V45_V225[array_idx])
        A46 = np.mean(V36_V225[array_idx])
        A47 = np.mean(V27_V225[array_idx])

        b5 = np.mean(chi_945_num[array_idx])
        A51 = np.conj(A15)
        A52 = np.conj(A25)
        A53 = np.conj(A35)
        A54 = np.conj(A45)
        A55 = np.mean(V45_V45[array_idx])
        A56 = np.mean(V36_V45[array_idx])
        A57 = np.mean(V27_V45[array_idx])

        b6 = np.mean(chi_936_num[array_idx])
        A61 = np.conj(A16)
        A62 = np.conj(A26)
        A63 = np.conj(A36)
        A64 = np.conj(A46)
        A65 = np.conj(A56)
        A66 = np.mean(V36_V36[array_idx])
        A67 = np.mean(V27_V36[array_idx])

        b7 = np.mean(chi_927_num[array_idx])
        A71 = np.conj(A17)
        A72 = np.conj(A27)
        A73 = np.conj(A37)
        A74 = np.conj(A47)
        A75 = np.conj(A57)
        A76 = np.mean(A67)
        A77 = np.mean(V27_V27[array_idx])

        if flag == 0:
            array_lhs = np.array([b1, b2, b3, b4, b5, b6, b7])
            array_rhs = np.array([[A11, 0, 0, 0, 0, 0, 0],
                                  [0, A22, 0, 0, 0, 0, 0],
                                  [0, 0, A33, 0, 0, 0, 0],
                                  [0, 0, 0, A44, 0, 0, 0],
                                  [0, 0, 0, 0, A55, 0, 0],
                                  [0, 0, 0, 0, 0, A66, 0],
                                  [0, 0, 0, 0, 0, 0, A77],
                                 ], dtype=np.cfloat)
        elif flag == 1:
            array_lhs = np.array([b1, b2, b3, b4, b5, b6, b7])
            array_rhs = np.array([[A11, A12, A13, A14, A15, A16, A17],
                                  [A21, A22, A23, A24, A25, A26, A27],
                                  [A31, A32, A33, A34, A35, A36, A37],
                                  [A41, A42, A43, A44, A45, A46, A47],
                                  [A51, A52, A53, A54, A55, A56, A57],
                                  [A61, A62, A63, A64, A65, A66, A67],
                                  [A71, A72, A73, A74, A75, A76, A77],
                                 ], dtype=np.cfloat)
        elif flag == 2:
            array_lhs = np.array([b1, b2, b3, b4, b5])
            array_rhs = np.array([[A11, A12, A13, A14, A15],
                                  [A21, A22, A23, A24, A25],
                                  [A31, A32, A33, A34, A35],
                                  [A41, A42, A43, A44, A45],
                                  [A51, A52, A53, A54, A55],
                                 ], dtype=np.cfloat)
        chi_9 = np.linalg.solve(array_rhs, array_lhs)
        chi_9333_JK[iev] = chi_9[0]
        chi_92223_JK[iev] = chi_9[1]
        chi_9234_JK[iev] = chi_9[2]
        chi_9225_JK[iev] = chi_9[3]
        chi_945_JK[iev] = chi_9[4]
        if flag < 2:
            chi_936_JK[iev] = chi_9[5]
            chi_927_JK[iev] = chi_9[6]

    # compute the real part of the response coefficients
    chi_9333_mean = np.mean(np.real(chi_9333_JK))
    chi_9333_err = np.sqrt((nev - 1.)/nev
                           *np.sum((np.real(chi_9333_JK) - chi_9333_mean)**2))
    chi_92223_mean = np.mean(np.real(chi_92223_JK))
    chi_92223_err = np.sqrt((nev - 1.)/nev
                           *np.sum((np.real(chi_92223_JK) - chi_92223_mean)**2.))
    chi_9234_mean = np.mean(np.real(chi_9234_JK))
    chi_9234_err = np.sqrt((nev - 1.)/nev
                          *np.sum((np.real(chi_9234_JK) - chi_9234_mean)**2.))
    chi_9225_mean = np.mean(np.real(chi_9225_JK))
    chi_9225_err = np.sqrt((nev - 1.)/nev
                          *np.sum((np.real(chi_9225_JK) - chi_9225_mean)**2.))
    chi_945_mean = np.mean(np.real(chi_945_JK))
    chi_945_err = np.sqrt((nev - 1.)/nev
                          *np.sum((np.real(chi_945_JK) - chi_945_mean)**2.))
    chi_936_mean = np.mean(np.real(chi_936_JK))
    chi_936_err = np.sqrt((nev - 1.)/nev
                          *np.sum((np.real(chi_936_JK) - chi_936_mean)**2.))
    chi_927_mean = np.mean(np.real(chi_927_JK))
    chi_927_err = np.sqrt((nev - 1.)/nev
                          *np.sum((np.real(chi_927_JK) - chi_927_mean)**2.))

    results = [chi_9333_mean, chi_9333_err, chi_92223_mean, chi_92223_err,
               chi_9234_mean, chi_9234_err, chi_9225_mean, chi_9225_err,
               chi_945_mean, chi_945_err,
               chi_936_mean, chi_936_err, chi_927_mean, chi_927_err]

    # compute sqrt<V9L^2>
    V9L_1 = (V9_1 - chi_9333_mean*(V3_1**3.)
             - chi_92223_mean*(V2_1**3*V3_1)
             - chi_9234_mean*(V2_1*V3_1*V4L_1)
             - chi_9225_mean*(V2_1*V2_1*V5L_1)
             - chi_945_mean*(V4L_1*V5L_1)
             - chi_936_mean*(V3_1*V6L_1)
             - chi_927_mean*(V2_1*V7L_1))
    V9L_2 = (V9_2 - chi_9333_mean*(V3_2**3.)
             - chi_92223_mean*(V2_2**3*V3_2)
             - chi_9234_mean*(V2_2*V3_2*V4L_2)
             - chi_9225_mean*(V2_2*V2_2*V5L_2)
             - chi_945_mean*(V4L_2*V5L_2)
             - chi_936_mean*(V3_2*V6L_2)
             - chi_927_mean*(V2_2*V7L_2))
    v9L_rms = np.nan_to_num(np.sqrt(np.mean(np.real(V9L_1*np.conj(V9L_2)))))
    v9L_err = np.nan_to_num(np.std(np.real(V9L_1*np.conj(V9L_2)))/np.sqrt(nev)
                            / (2*v9L_rms))

    results += [v9L_rms, v9L_err]

    # compute the imaginary part of the response coefficients
    chi_9333_mean = np.mean(np.imag(chi_9333_JK))
    chi_9333_err = np.sqrt((nev - 1.)/nev
                           *np.sum((np.imag(chi_9333_JK) - chi_9333_mean)**2))
    chi_92223_mean = np.mean(np.imag(chi_92223_JK))
    chi_92223_err = np.sqrt((nev - 1.)/nev
                           *np.sum((np.imag(chi_92223_JK) - chi_92223_mean)**2.))
    chi_9234_mean = np.mean(np.imag(chi_9234_JK))
    chi_9234_err = np.sqrt((nev - 1.)/nev
                          *np.sum((np.imag(chi_9234_JK) - chi_9234_mean)**2.))
    chi_9225_mean = np.mean(np.imag(chi_9225_JK))
    chi_9225_err = np.sqrt((nev - 1.)/nev
                          *np.sum((np.imag(chi_9225_JK) - chi_9225_mean)**2.))
    chi_945_mean = np.mean(np.imag(chi_945_JK))
    chi_945_err = np.sqrt((nev - 1.)/nev
                          *np.sum((np.imag(chi_945_JK) - chi_945_mean)**2.))
    chi_936_mean = np.mean(np.imag(chi_936_JK))
    chi_936_err = np.sqrt((nev - 1.)/nev
                          *np.sum((np.imag(chi_936_JK) - chi_936_mean)**2.))
    chi_927_mean = np.mean(np.imag(chi_927_JK))
    chi_927_err = np.sqrt((nev - 1.)/nev
                          *np.sum((np.imag(chi_927_JK) - chi_927_mean)**2.))

    results += [chi_9333_mean, chi_9333_err, chi_92223_mean, chi_92223_err,
                chi_9234_mean, chi_9234_err, chi_9225_mean, chi_9225_err,
                chi_945_mean, chi_945_err,
                chi_936_mean, chi_936_err, chi_927_mean, chi_927_err]

    dN_mean = np.real(np.mean(vn_data_array1[:, 0] + vn_data_array2[:, 0]))
    dN_err = (np.std(vn_data_array1[:, 0] + vn_data_array2[:, 0])
              / np.sqrt(nev))
    if path.isfile(outputFileNameV9):
        f = open(outputFileNameV9, 'a')
    else:
        f = open(outputFileNameV9, 'w')
        f.write("# cen  Nch  Re{chi_9333}  Re{chi_92223}  Re{chi_9234}  "
                + "Re{chi_9225}  Re{chi_945}  Re{chi_936}  Re{chi_927}  "
                + "v9L_rms  Im{chi_9333}  Im{chi_92223}  Im{chi_9234}  "
                + "Im{chi_9225}  Im{chi_945}  Im{chi_936}  Im{chi_927}\n")
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
    if centralityCutList[icen+1] < centralityCutList[icen]: continue
    selected_events_list = []

    dN_dy_cut_high = dNdyList[
        int(len(dNdyList)*centralityCutList[icen]/100.)
    ]
    dN_dy_cut_low  = dNdyList[
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
    if nev <= 0: continue

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

    calculateNonLinearResponseV9_2sub(QnArr2, QnArr3,
                                      "nonLinearV9_2sub_diag.dat", 0, cenLabel)
    calculateNonLinearResponseV9_2sub(QnArr2, QnArr3,
                                      "nonLinearV9_2sub_full.dat", 1, cenLabel)
    calculateNonLinearResponseV9_2sub(QnArr2, QnArr3,
                                      "nonLinearV9_2sub_sub5.dat", 2, cenLabel)
