#!/usr/bin/env python3

import pickle
import numpy as np
from os import path


Reg_centrality_cut_list = [0., 5., 10., 20., 30., 40., 50.,
                           60., 70., 80., 90., 100.]
PHOBOS_cen_list = [0., 6., 15., 25., 35., 45., 55.]  # PHOBOS AuAu 200
SPS_cen_list    = [5., 12.5, 23.5, 33.5, 43.5]       # SPS PbPb
PHENIX_cen_list = [0., 20., 40., 60., 88.]           # PHENIX dAu
STAR_cen_list   = [0., 10., 40., 80]                 # STAR v1
ALICE_pp_list   = [0., 100., 0., 1., 5.,
                   0., 5., 10., 15, 20., 30., 40., 50., 70., 100.]
centralityCutList = Reg_centrality_cut_list
dNcutList = []    # pre-defined Nch cut if simulation is not minimum bias


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
            gamma_1 = -6\sqrt{2}*vn{4}^2*(vn{4} - vn{6})
                                         /(vn{2}^2 - vn{4}^2)^(3/2)

        we will use Jackknife resampling method to estimate
        the statistical error
    """
    nev = len(vn_data_array[:, 0])
    dN = np.real(vn_data_array[:, -1])
    Q2 = dN*vn_data_array[:, 3]
    Q3 = dN*vn_data_array[:, 4]
    Q4 = dN*vn_data_array[:, 5]
    Q5 = dN*vn_data_array[:, 6]
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
                            /(v2_2**2. + v2_4**2. + 1e-15)))
            if C_2_6 > 0.:
                v2_6 = (C_2_6/4.)**(1./6.)
                r26_array[iev] = v2_6/v2_4
                gamma1_array[iev] = (-6.*np.sqrt(2)*(v2_4**2.)*(v2_4 - v2_6)
                                     /(v2_2**2. - v2_4**2.)**(1.5))


        C_3_2 = np.mean(Q3_N2[array_idx])/np.mean(N2_weight[array_idx])
        C34_tmp = np.mean(Q3_4[array_idx])/np.mean(N4_weight[array_idx])
        C_3_4 = C34_tmp - 2.*C_3_2**2.
        C3_4_array[iev] = C_3_4
        if C_3_4 < 0. and C_3_2 > 0.:
            v3_4 = (-C_3_4)**0.25
            v3_2 = np.sqrt(C_3_2)
            r3_array[iev] = v3_4/v3_2
            F3_array[iev] = np.sqrt((v3_2**2. - v3_4**2.)
                                    /(v3_2**2. + v3_4**2. + 1e-15))

    dN_mean = np.real(np.mean(vn_data_array[:, 0]))
    dN_err = np.std(vn_data_array[:, 0])/np.sqrt(nev)
    C2_4_mean = np.mean(C2_4_array)
    C2_4_err  = np.sqrt((nev - 1.)/nev*np.sum((C2_4_array - C2_4_mean)**2.))
    C3_4_mean = np.mean(C3_4_array)
    C3_4_err  = np.sqrt((nev - 1.)/nev*np.sum((C3_4_array - C3_4_mean)**2.))

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
        f.write(
        "# cen  Nch  vn{6}/vn{4}  (vn{6}/vn{4})_err  gamma_1  gamma_1_err\n")
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
    Q2_N4 = (np.real(Q2_1*np.conj(Q2_1)*Q2_1*np.conj(Q2_2))
             - np.real(Q4_1*np.conj(Q2_2)*np.conj(Q2_2))
             - np.real(Q2_1*Q2_1*np.conj(Q4_2))
             + np.real(Q4_1*np.conj(Q4_2)))

    Q3_N4 = (np.real(Q3_1*np.conj(Q3_1)*Q3_1*np.conj(Q3_2))
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
                                    /(v2_2**2. + v2_4**2. + 1e-15))

        C_3_2 = np.mean(Q3_N2[array_idx])/np.mean(N2_weight[array_idx])
        C34_tmp = np.mean(Q3_N4[array_idx])/np.mean(N4_weight[array_idx])
        C_3_4 = C34_tmp - 2.*C_3_2**2.
        C3_4_array[iev] = C_3_4
        if C_3_4 < 0. and C_3_2 > 0.:
            v3_4 = (-C_3_4)**0.25
            v3_2 = np.sqrt(C_3_2)
            r3_array[iev] = v3_4/v3_2
            F3_array[iev] = np.sqrt((v3_2**2. - v3_4**2.)
                                    /(v3_2**2. + v3_4**2. + 1e-15))

    C2_4_mean = np.mean(C2_4_array)
    C2_4_err  = np.sqrt((nev - 1.)/nev*np.sum((C2_4_array - C2_4_mean)**2.))
    C3_4_mean = np.mean(C3_4_array)
    C3_4_err  = np.sqrt((nev - 1.)/nev*np.sum((C3_4_array - C3_4_mean)**2.))

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
              /np.sqrt(nev))
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
              /np.sqrt(nev))
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


with open("QnVectors.pickle", "rb") as pf:
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
        QnArr1.append(data[event_name]['STAR_eta_-0p5_0p5'])
        QnArr2.append(data[event_name]['STAR_eta_-1_-0p5'])
        QnArr3.append(data[event_name]['STAR_eta_0p5_1'])
    QnArr1 = np.array(QnArr1)
    QnArr2 = np.array(QnArr2)
    QnArr3 = np.array(QnArr3)

    calcualte_vn_2_with_gap(QnArr2, QnArr3, "vn2_sub.dat", cenLabel)
    calculate_vn4_2sub(QnArr2, QnArr3, "vn4_2sub.dat",
                       "vn4_2sub_over_vn2.dat", cenLabel)
    calculate_vn4_vn6(QnArr1, "vn4.dat", "vn4_over_vn2.dat",
                      "vn6_over_vn4.dat", cenLabel)
