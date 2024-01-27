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


def calculateNonLinearResponseV6_2sub(vn_data_array1, vn_data_array2,
                                      V4L_1, V4L_2):
    """
        This function computes the linear coefficients
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
    return(V6L_1, V6L_2)


def calculateNonLinearResponseV8_2sub(vn_data_array1, vn_data_array2,
                                      outputFileNameV8, flag, cenLabel):
    """
        This function computes the non-linear response coefficients
        for V8 = V8L
                 + chi_82222 V2^4
                 + chi_8233 V2*V3^2
                 + chi_8224 V2^2*V4L
                 + chi_835 V3*V5L
                 + chi_826 V2*V6L
                 + chi_844 V4L^2
    """
    nev = len(vn_data_array1[:, 0])

    V2_1 = vn_data_array1[:, 3]
    V3_1 = vn_data_array1[:, 4]
    V4_1 = vn_data_array1[:, 5]
    V5_1 = vn_data_array1[:, 6]
    V6_1 = vn_data_array1[:, 7]
    V8_1 = vn_data_array1[:, 9]

    V2_2 = vn_data_array2[:, 3]
    V3_2 = vn_data_array2[:, 4]
    V4_2 = vn_data_array2[:, 5]
    V5_2 = vn_data_array2[:, 6]
    V6_2 = vn_data_array2[:, 7]
    V8_2 = vn_data_array2[:, 9]

    chi_422 = calculate_chi_422(vn_data_array1, vn_data_array2)
    V4L_1 = V4_1 - chi_422*(V2_1**2)
    V4L_2 = V4_2 - chi_422*(V2_2**2)

    chi_523 = calculate_chi_523(vn_data_array1, vn_data_array2)
    V5L_1 = V5_1 - chi_523*(V2_1*V3_1)
    V5L_2 = V5_2 - chi_523*(V2_2*V3_2)

    V6L_1, V6L_2 = calculateNonLinearResponseV6_2sub(
                        vn_data_array1, vn_data_array2, V4L_1, V4L_2)

    conjV2222_1 = np.conj(V2_1**4)
    conjV2222_2 = np.conj(V2_2**4)
    chi_82222_num = V8_1*conjV2222_2 + V8_2*conjV2222_1
    V2222_V2222 = 2*np.real((V2_1**4)*conjV2222_2)
    V233_V2222 = V2_1*V3_1*V3_1*conjV2222_2 + V2_2*V3_2*V3_2*conjV2222_1
    V224_V2222 = V2_1*V2_1*V4L_1*conjV2222_2 + V2_2*V2_2*V4L_2*conjV2222_1
    V35_V2222 = V3_1*V5L_1*conjV2222_2 + V3_2*V5L_2*conjV2222_1
    V26_V2222 = V2_1*V6L_1*conjV2222_2 + V2_2*V6L_2*conjV2222_1
    V44_V2222 = V4L_1*V4L_1*conjV2222_2 + V4L_2*V4L_2*conjV2222_1

    conjV233_1 = np.conj(V2_1*V3_1*V3_1)
    conjV233_2 = np.conj(V2_2*V3_2*V3_2)
    chi_8233_num = V8_1*conjV233_2 + V8_2*conjV233_1
    V233_V233 = 2*np.real(V2_1*V3_1*V3_1*conjV233_2)
    V224_V233 = V2_1*V2_1*V4L_1*conjV233_2 + V2_2*V2_2*V4L_2*conjV233_1
    V35_V233 = V3_1*V5L_1*conjV233_2 + V3_2*V5L_2*conjV233_1
    V26_V233 = V2_1*V6L_1*conjV233_2 + V2_2*V6L_2*conjV233_1
    V44_V233 = V4L_1*V4L_1*conjV233_2 + V4L_2*V4L_2*conjV233_1

    conjV224_1 = np.conj(V2_1*V2_1*V4L_1)
    conjV224_2 = np.conj(V2_2*V2_2*V4L_2)
    chi_8224_num = V8_1*conjV224_2 + V8_2*conjV224_1
    V224_V224 = 2*np.real(V2_1*V2_1*V4L_1*conjV224_2)
    V35_V224 = V3_1*V5L_1*conjV224_2 + V3_2*V5L_2*conjV224_1
    V26_V224 = V2_1*V6L_1*conjV224_2 + V2_2*V6L_2*conjV224_1
    V44_V224 = V4L_1*V4L_1*conjV224_2 + V4L_2*V4L_2*conjV224_1

    conjV35_1 = np.conj(V3_1*V5L_1)
    conjV35_2 = np.conj(V3_2*V5L_2)
    chi_835_num = V8_1*conjV35_2 + V8_2*conjV35_1
    V35_V35 = 2*np.real(V3_1*V5L_1*conjV35_2)
    V26_V35 = V2_1*V6L_1*conjV35_2 + V2_2*V6L_2*conjV35_1
    V44_V35 = V4L_1*V4L_1*conjV35_2 + V4L_2*V4L_2*conjV35_1

    conjV26_1 = np.conj(V2_1*V6L_1)
    conjV26_2 = np.conj(V2_2*V6L_2)
    chi_826_num = V8_1*conjV26_2 + V8_2*conjV26_1
    V26_V26 = 2*np.real(V2_1*V6L_1*conjV26_2)
    V44_V26 = V4L_1*V4L_1*conjV26_2 + V4L_2*V4L_2*conjV26_1

    conjV44_1 = np.conj(V4L_1**2)
    conjV44_2 = np.conj(V4L_2**2)
    chi_844_num = V8_1*conjV44_2 + V8_2*conjV44_1
    V44_V44 = 2*np.real(V4L_1**2*conjV44_2)

    chi_82222_JK = np.zeros(nev, dtype=complex)
    chi_8233_JK = np.zeros(nev, dtype=complex)
    chi_8224_JK = np.zeros(nev, dtype=complex)
    chi_835_JK = np.zeros(nev, dtype=complex)
    chi_826_JK = np.zeros(nev, dtype=complex)
    chi_844_JK = np.zeros(nev, dtype=complex)
    for iev in range(nev):
        array_idx = [True]*nev
        array_idx[iev] = False
        array_idx = np.array(array_idx)

        b1 = np.mean(chi_82222_num[array_idx])
        A11 = np.mean(V2222_V2222[array_idx])
        A12 = np.mean(V233_V2222[array_idx])
        A13 = np.mean(V224_V2222[array_idx])
        A14 = np.mean(V35_V2222[array_idx])
        A15 = np.mean(V26_V2222[array_idx])
        A16 = np.mean(V44_V2222[array_idx])

        b2 = np.mean(chi_8233_num[array_idx])
        A21 = np.conj(A12)
        A22 = np.mean(V233_V233[array_idx])
        A23 = np.mean(V224_V233[array_idx])
        A24 = np.mean(V35_V233[array_idx])
        A25 = np.mean(V26_V233[array_idx])
        A26 = np.mean(V44_V44[array_idx])

        b3 = np.mean(chi_8224_num[array_idx])
        A31 = np.conj(A13)
        A32 = np.conj(A23)
        A33 = np.mean(V224_V224[array_idx])
        A34 = np.mean(V35_V224[array_idx])
        A35 = np.mean(V26_V224[array_idx])
        A36 = np.mean(V44_V224[array_idx])

        b4 = np.mean(chi_835_num[array_idx])
        A41 = np.conj(A14)
        A42 = np.conj(A24)
        A43 = np.conj(A34)
        A44 = np.mean(V35_V35[array_idx])
        A45 = np.mean(V26_V35[array_idx])
        A46 = np.mean(V44_V35[array_idx])

        b5 = np.mean(chi_826_num[array_idx])
        A51 = np.conj(A15)
        A52 = np.conj(A25)
        A53 = np.conj(A35)
        A54 = np.conj(A45)
        A55 = np.mean(V26_V26[array_idx])
        A56 = np.mean(V44_V26[array_idx])

        b6 = np.mean(chi_844_num[array_idx])
        A61 = np.conj(A16)
        A62 = np.conj(A26)
        A63 = np.conj(A36)
        A64 = np.conj(A46)
        A65 = np.conj(A56)
        A66 = np.mean(V44_V44[array_idx])

        if flag == 0:
            array_lhs = np.array([b1, b2, b3, b4, b5, b6])
            array_rhs = np.array([[A11, 0, 0, 0, 0, 0],
                                  [0, A22, 0, 0, 0, 0],
                                  [0, 0, A33, 0, 0, 0],
                                  [0, 0, 0, A44, 0, 0],
                                  [0, 0, 0, 0, A55, 0],
                                  [0, 0, 0, 0, 0, A66],
                                 ], dtype=np.cfloat)
        elif flag == 1:
            array_lhs = np.array([b1, b2, b3, b4, b5, b6])
            array_rhs = np.array([[A11, A12, A13, A14, A15, A16],
                                  [A21, A22, A23, A24, A25, A26],
                                  [A31, A32, A33, A34, A35, A36],
                                  [A41, A42, A43, A44, A45, A46],
                                  [A51, A52, A53, A54, A55, A56],
                                  [A61, A62, A63, A64, A65, A66],
                                 ], dtype=np.cfloat)
        elif flag == 2:
            array_lhs = np.array([b1, b2, b3])
            array_rhs = np.array([[A11, A12, A13],
                                  [A21, A22, A23],
                                  [A31, A32, A33],
                                 ], dtype=np.cfloat)
        chi_8 = np.linalg.solve(array_rhs, array_lhs)

        chi_82222_JK[iev] = chi_8[0]
        chi_8233_JK[iev] = chi_8[1]
        chi_8224_JK[iev] = chi_8[2]
        if flag < 2:
            chi_835_JK[iev] = chi_8[3]
            chi_826_JK[iev] = chi_8[4]
            chi_844_JK[iev] = chi_8[5]

    # compute the real part of the response coefficients
    chi_82222_mean = np.mean(np.real(chi_82222_JK))
    chi_82222_err = np.sqrt((nev - 1.)/nev
                           *np.sum((np.real(chi_82222_JK) - chi_82222_mean)**2))
    chi_8233_mean = np.mean(np.real(chi_8233_JK))
    chi_8233_err = np.sqrt((nev - 1.)/nev
                          *np.sum((np.real(chi_8233_JK) - chi_8233_mean)**2.))
    chi_8224_mean = np.mean(np.real(chi_8224_JK))
    chi_8224_err = np.sqrt((nev - 1.)/nev
                          *np.sum((np.real(chi_8224_JK) - chi_8224_mean)**2.))
    chi_835_mean = np.mean(np.real(chi_835_JK))
    chi_835_err = np.sqrt((nev - 1.)/nev
                          *np.sum((np.real(chi_835_JK) - chi_835_mean)**2.))
    chi_826_mean = np.mean(np.real(chi_826_JK))
    chi_826_err = np.sqrt((nev - 1.)/nev
                          *np.sum((np.real(chi_826_JK) - chi_826_mean)**2.))
    chi_844_mean = np.mean(np.real(chi_844_JK))
    chi_844_err = np.sqrt((nev - 1.)/nev
                          *np.sum((np.real(chi_844_JK) - chi_844_mean)**2.))

    results = [chi_82222_mean, chi_82222_err, chi_8233_mean, chi_8233_err,
               chi_8224_mean, chi_8224_err, chi_835_mean, chi_835_err,
               chi_826_mean, chi_826_err, chi_844_mean, chi_844_err]

    # compute sqrt<V8L^2>
    V8L_1 = (V8_1 - chi_82222_mean*(V2_1**4.)
             - chi_8233_mean*(V2_1*V3_1**2.)
             - chi_8224_mean*(V2_1**2.*V4L_1)
             - chi_835_mean*(V3_1*V5L_1)
             - chi_826_mean*(V2_1*V6L_1)
             - chi_844_mean*(V4L_1**2.))
    V8L_2 = (V8_2 - chi_82222_mean*(V2_2**4.)
             - chi_8233_mean*(V2_2*V3_2**2.)
             - chi_8224_mean*(V2_2**2.*V4L_2)
             - chi_835_mean*(V3_2*V5L_2)
             - chi_826_mean*(V2_2*V6L_2)
             - chi_844_mean*(V4L_2**2.))
    v8L_rms = np.nan_to_num(np.sqrt(np.mean(np.real(V8L_1*np.conj(V8L_2)))))
    v8L_err = np.nan_to_num(np.std(np.real(V8L_1*np.conj(V8L_2)))/np.sqrt(nev)
                            / (2*v8L_rms))

    results += [v8L_rms, v8L_err]

    # compute the imaginary part of the response coefficients
    chi_82222_mean = np.mean(np.imag(chi_82222_JK))
    chi_82222_err = np.sqrt((nev - 1.)/nev
                           *np.sum((np.imag(chi_82222_JK) - chi_82222_mean)**2))
    chi_8233_mean = np.mean(np.imag(chi_8233_JK))
    chi_8233_err = np.sqrt((nev - 1.)/nev
                          *np.sum((np.imag(chi_8233_JK) - chi_8233_mean)**2.))
    chi_8224_mean = np.mean(np.imag(chi_8224_JK))
    chi_8224_err = np.sqrt((nev - 1.)/nev
                          *np.sum((np.imag(chi_8224_JK) - chi_8224_mean)**2.))
    chi_835_mean = np.mean(np.imag(chi_835_JK))
    chi_835_err = np.sqrt((nev - 1.)/nev
                          *np.sum((np.imag(chi_835_JK) - chi_835_mean)**2.))
    chi_826_mean = np.mean(np.imag(chi_826_JK))
    chi_826_err = np.sqrt((nev - 1.)/nev
                          *np.sum((np.imag(chi_826_JK) - chi_826_mean)**2.))
    chi_844_mean = np.mean(np.imag(chi_844_JK))
    chi_844_err = np.sqrt((nev - 1.)/nev
                          *np.sum((np.imag(chi_844_JK) - chi_844_mean)**2.))

    results += [chi_82222_mean, chi_82222_err, chi_8233_mean, chi_8233_err,
                chi_8224_mean, chi_8224_err, chi_835_mean, chi_835_err,
                chi_826_mean, chi_826_err, chi_844_mean, chi_844_err]

    dN_mean = np.real(np.mean(vn_data_array1[:, 0] + vn_data_array2[:, 0]))
    dN_err = (np.std(vn_data_array1[:, 0] + vn_data_array2[:, 0])
              / np.sqrt(nev))
    if path.isfile(outputFileNameV8):
        f = open(outputFileNameV8, 'a')
    else:
        f = open(outputFileNameV8, 'w')
        f.write("# cen  Nch  Re{chi_82222}  Re{chi_8233}  Re{chi_8224}  "
                + "Re{chi_835}  Re{chi_826}  Re{chi_844}  "
                + "v8L_rms  Im{chi_82222}  Im{chi_8233}  Im{chi_8224}  "
                + "Im{chi_835}  Im{chi_826}  Im{chi_844}\n")
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

    calculateNonLinearResponseV8_2sub(QnArr2, QnArr3,
                                      "nonLinearV8_2sub_diag.dat", 0, cenLabel)
    calculateNonLinearResponseV8_2sub(QnArr2, QnArr3,
                                      "nonLinearV8_2sub_full.dat", 1, cenLabel)
    calculateNonLinearResponseV8_2sub(QnArr2, QnArr3,
                                      "nonLinearV8_2sub_sub3.dat", 2, cenLabel)
