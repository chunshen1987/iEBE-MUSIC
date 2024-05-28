#!/usr/bin/env python3

import sys
from os import path
import pickle
import numpy as np


def help_message():
    print("Usage: {0} database_file".format(sys.argv[0]))
    exit(0)


centralityRange = 1.
Reg_centrality_cut_list = [0., 5., 10., 20., 30., 40., 50.,
                           60., 70., 80., 90., 100.]
centralityCutList = Reg_centrality_cut_list
# centralityCutList = [0, 1, 2, 3, 4, 6, 8, 10, 15, 20, 30, 40, 50, 60,
#                      70, 80, 90, 100]
dNcutList = []    # pre-defined Nch cut if simulation is not minimum bias


def calculate_v0pT(pTArr, poiSP, refSP, pTmin: float, pTmax: float,
                   outputFileName: str) -> None :
    """
        this function compute the v_0(p_T) according to the work
        arXiv: 2004.00690 [nucl-th]
    """
    pTInterp = np.linspace(pTmin, pTmax, 41)
    dpTInterp = pTInterp[1] - pTInterp[0]
    nev, npT = refSP.shape
    Nref = np.zeros([nev])
    PTref = np.zeros([nev])
    for iev in range(nev):
        refSPInterp = np.exp(np.interp(pTInterp, pTArr,
                                       np.log(refSP[iev, :]+1e-30)))
        Nref[iev] = np.sum(refSPInterp*pTInterp)*dpTInterp*2*np.pi
        PTref[iev] = np.sum(refSPInterp*pTInterp**2.)*dpTInterp*2*np.pi
    deltaNref = Nref - np.mean(Nref)
    deltaPT = PTref - np.mean(PTref)
    delta_poi_SP = poiSP - np.mean(poiSP, axis=0)

    v0pT_arr = np.zeros([nev, npT])
    for iev in range(nev):
        array_idx      = [True]*nev
        array_idx[iev] = False
        array_idx      = np.array(array_idx)

        sigmaNsq = np.mean(deltaNref[array_idx]**2)
        delta_hat_PT = (deltaPT[array_idx]
                        - np.mean(deltaPT[array_idx]*deltaNref[array_idx])
                          /sigmaNsq*deltaNref[array_idx])
        deltaNref_local = deltaNref[array_idx].reshape(-1, 1)
        deltahat_poi_SP = (delta_poi_SP[array_idx, :]
                           - np.mean(delta_poi_SP[array_idx, :]
                                     *deltaNref_local, axis=0)
                             /sigmaNsq*deltaNref_local)
        sigmahatPT = np.sqrt(np.mean(delta_hat_PT**2.))
        delta_hat_PT_local = np.copy(delta_hat_PT).reshape(-1, 1)
        v0pT_arr[iev, :] = (np.mean(deltahat_poi_SP*delta_hat_PT_local, axis=0)
                            /sigmahatPT/np.mean(poiSP[array_idx, :], axis=0))
    v0pT_mean = np.mean(v0pT_arr, axis=0)
    v0pT_err = np.sqrt((nev - 1)/nev
                       *np.sum((v0pT_arr - v0pT_mean)**2., axis=0))
    results = np.array([pTArr, v0pT_mean, v0pT_err]).transpose()
    np.savetxt(outputFileName, results, fmt="%.4e", delimiter="  ",
               header="pT (GeV)  v0(pT)  v0(pT)_err")


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

    cenLabel = "{:d}-{:d}".format(
        int(centralityCutList[icen]*centralityRange),
        int(centralityCutList[icen+1]*centralityRange))
    print("analysis {}%-{}% nev = {}...".format(
            centralityCutList[icen]*centralityRange,
            centralityCutList[icen+1]*centralityRange, nev))
    print(f"dNdy: {dN_dy_cut_low:.2f} - {dN_dy_cut_high:.2f}")

    pTArr = data[selected_events_list[0]]['pT']
    charged_Sp = []
    pion_Sp = []
    kaon_Sp = []
    proton_Sp = []
    for event_name in selected_events_list:
        charged_Sp.append(data[event_name]['ch_sp'])
        pion_Sp.append(data[event_name]['pi+_sp'] + data[event_name]['pi-_sp'])
        kaon_Sp.append(data[event_name]['K+_sp'] + data[event_name]['K-_sp'])
        proton_Sp.append(data[event_name]['p_sp']
                         + data[event_name]['pbar_sp'])
    charged_Sp = np.array(charged_Sp)
    pion_Sp = np.array(pion_Sp)
    kaon_Sp = np.array(kaon_Sp)
    proton_Sp = np.array(proton_Sp)

    calculate_v0pT(pTArr, charged_Sp, charged_Sp, 0., 3.,
                   f"v0pT_ChargedHadron_C{cenLabel}.dat")
    calculate_v0pT(pTArr, pion_Sp, charged_Sp, 0., 3.,
                   f"v0pT_pion_C{cenLabel}.dat")
    calculate_v0pT(pTArr, kaon_Sp, charged_Sp, 0., 3.,
                   f"v0pT_kaon_C{cenLabel}.dat")
    calculate_v0pT(pTArr, proton_Sp, charged_Sp, 0., 3.,
                   f"v0pT_proton_C{cenLabel}.dat")
