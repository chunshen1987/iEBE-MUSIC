#!/usr/bin/env python3

import sys
from os import path
import pickle
import numpy as np


def help_message():
    print("Usage: {0} database_file".format(sys.argv[0]))
    exit(0)


centralityRange = 1.
Reg_centrality_cut_list = [
    0., 5., 10., 20., 30., 40., 50., 60., 70., 80., 90., 100.
]
centralityCutList = Reg_centrality_cut_list
# centralityCutList = [0, 1, 2, 3, 4, 6, 8, 10, 15, 20, 30, 40, 50, 60,
#                      70, 80, 90, 100]
dNcutList = []  # pre-defined Nch cut if simulation is not minimum bias


def calculate_pTSpectra(pTArr, poiSP, outputFileName: str) -> None:
    """
        this function compute the event averaged pT spectra
    """
    nev, npt = poiSP.shape
    pTSp_mean = np.mean(poiSP, axis=0)
    pTSp_err = np.sqrt(pTSp_mean)/np.sqrt(nev)
    results = np.array([pTArr, pTSp_mean, pTSp_err]).transpose()
    np.savetxt(outputFileName,
               results,
               fmt="%.4e",
               delimiter="  ",
               header="pT (GeV)  dN/(2pi pT dpT dy)  dN/(2pi pT dpT dy)err",)


try:
    database_file = str(sys.argv[1])
except IndexError:
    help_message()

with open(database_file, "rb") as pf:
    data = pickle.load(pf)

dNdyDict = {}
for event_name in data.keys():
    if 'global' not in event_name:
        Nch = np.real(
              data[event_name]['ALICE_V0A_eta_2p8_5p1_pT_0_4'][0]
            + data[event_name]['ALICE_V0C_eta_-3p7_-1p7_pT_0_4'][0]
        )
        dNdyDict[event_name] = Nch
dNdyList = -np.sort(-np.array(list(dNdyDict.values())))
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

    for event_name in dNdyDict.keys():
        if (dNdyDict[event_name] > dN_dy_cut_low
                and dNdyDict[event_name] <= dN_dy_cut_high):
            selected_events_list.append(event_name)

    nev = len(selected_events_list)
    if nev <= 0:
        continue

    cenLabel = "{:d}-{:d}".format(
        int(centralityCutList[icen]*centralityRange),
        int(centralityCutList[icen + 1]*centralityRange))
    print("analysis {}%-{}% nev = {}...".format(
        centralityCutList[icen]*centralityRange,
        centralityCutList[icen + 1]*centralityRange, nev))
    print(f"dNdy: {dN_dy_cut_low:.2f} - {dN_dy_cut_high:.2f}")

    pTArr = data['global']['pTArr']
    charged_Sp = []
    pion_Sp = []
    kaon_Sp = []
    proton_Sp = []
    for event_name in selected_events_list:
        charged_Sp.append(data[event_name]['ch_pTArr'][0, :])
        pion_Sp.append(data[event_name]['pi+_pTArr'][0, :]
                       + data[event_name]['pi-_pTArr'][0, :])
        kaon_Sp.append(data[event_name]['K+_pTArr'][0, :]
                       + data[event_name]['K-_pTArr'][0, :])
        proton_Sp.append(data[event_name]['p_pTArr'][0, :]
                         + data[event_name]['pbar_pTArr'][0, :])
    charged_Sp = np.array(charged_Sp)
    pion_Sp = np.array(pion_Sp)
    kaon_Sp = np.array(kaon_Sp)
    proton_Sp = np.array(proton_Sp)

    calculate_pTSpectra(pTArr, charged_Sp,
                        f"pTSp_ChargedHadron_C{cenLabel}.dat")
    calculate_pTSpectra(pTArr, pion_Sp,
                        f"pTSp_pion_C{cenLabel}.dat")
    calculate_pTSpectra(pTArr, kaon_Sp,
                        f"pTSp_kaon_C{cenLabel}.dat")
    calculate_pTSpectra(pTArr, proton_Sp,
                        f"pTSp_proton_C{cenLabel}.dat")
