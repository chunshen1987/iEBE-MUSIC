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
Reg_centrality_cut_list = [
    0., 5., 10., 20., 30., 40., 50., 60., 70., 80., 90., 100.
]
#centralityCutList = Reg_centrality_cut_list
centralityCutList = [
    0, 0.1, 0, 1, 2, 4, 6, 8, 10, 15, 20, 25, 30, 40, 50, 60, 70, 80, 90, 100
]
dNcutList = []  # pre-defined Nch cut if simulation is not minimum bias


def calcualte_eccn_2(eccn_data_array, outputFileName: str,
                     cenLabel: str) -> None:
    """
        this function computes ecc_n{2} and its stat. err.
    """
    nev = len(eccn_data_array[:, 0])
    en_array = eccn_data_array[:, ::2] + 1j*eccn_data_array[:, 1::2]

    corr = np.abs(en_array)**2
    en_2 = np.nan_to_num(np.sqrt(np.mean(corr, axis=0)))
    en_2_err = np.nan_to_num(np.std(corr, axis=0)/np.sqrt(nev)/2./en_2)

    if path.isfile(outputFileName):
        f = open(outputFileName, 'a')
    else:
        f = open(outputFileName, 'w')
        f.write("# cen  en{2}  en{2}_err (n = 1-6)\n")
    f.write("{:.3f}".format(cenLabel))
    for i in range(len(en_2)):
        f.write("  {:.5e}  {:.5e}".format(en_2[i], en_2_err[i]))
    f.write("\n")
    f.close()
    return


def calcualte_eccn_4(eccn_data_array, outputFileName: str,
                     cenLabel: str) -> None:
    """
        this function computes ecc_n{4} and its stat. err.
    """
    nev = len(eccn_data_array[:, 0])
    en_array = eccn_data_array[:, ::2] + 1j*eccn_data_array[:, 1::2]

    corr2 = np.abs(en_array)**2
    corr4 = np.abs(en_array)**4
    Cen4 = np.mean(corr4, axis=0) - 2.*np.mean(corr2, axis=0)**2
    Cen4_err = np.std(corr4, axis=0)/np.sqrt(nev)

    en4 = np.nan_to_num((-Cen4)**0.25)
    en4_err = np.nan_to_num(Cen4_err/4./(en4**3.))

    if path.isfile(outputFileName):
        f = open(outputFileName, 'a')
    else:
        f = open(outputFileName, 'w')
        f.write("# cen  en{4}  en{2}_err  Cen{4}  Cen{4}_err (n = 1-6)\n")
    f.write("{:.3f}".format(cenLabel))
    for i in range(len(en4)):
        f.write("  {:.5e}  {:.5e}  {:.5e}  {:.5e}".format(
            en4[i], en4_err[i], Cen4[i], Cen4_err[i]))
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

    QnArr2 = []
    QnArr3 = []
    ecc_n_Arr = []
    for event_name in selected_events_list:
        QnArr2.append(data[event_name]['ALICE_eta_-0p8_-0p4_pT_0p2_3'])
        QnArr3.append(data[event_name]['ALICE_eta_0p4_0p8_pT_0p2_3'])
        ecc_n_Arr.append(data[event_name]['ecc_n'])
    QnArr2 = np.array(QnArr2)
    QnArr3 = np.array(QnArr3)
    ecc_n_Arr = np.array(ecc_n_Arr)

    calcualte_eccn_2(ecc_n_Arr, "eccn2.dat", cenLabel)
    calcualte_eccn_4(ecc_n_Arr, "eccn4.dat", cenLabel)
