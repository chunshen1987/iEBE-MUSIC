#!/usr/bin/env python3

import h5py
import sys
import pickle
import numpy as np


def help_message():
    print("{0} database_file".format(sys.argv[0]))
    exit(0)


def calcualte_yield_and_meanpT(pT_low, pT_high, data):
    """
        this function calculates the pT-integrated particle yield and mean pT
        given pT range (pT_low, pT_high) for every event in the data
    """
    npT = 50
    pT_inte_array = np.linspace(pT_low, pT_high, npT)
    dpT = pT_inte_array[1] - pT_inte_array[0]
    dN_event = data[:, 1]
    pT_event = data[:, 0]
    dN_interp = np.exp(np.interp(pT_inte_array, pT_event,
                                 np.log(dN_event+1e-30)))
    N = 2.*np.pi*np.sum(dN_interp*pT_inte_array)*dpT
    meanpT = np.sum(dN_interp*pT_inte_array**2.)/np.sum(dN_interp*pT_inte_array)
    res_array = [N, meanpT]
    return(res_array)


try:
    database_file = str(sys.argv[1])
except:
    help_message()

h5_data = h5py.File(database_file, "r")
eventList = list(h5_data.keys())

outdata = {}
for ievent, event_i in enumerate(eventList):
    if ievent % 100 == 0:
        print("fetching event: {0} from the database {1} ...".format(
            event_i, database_file))
    eventGroup = h5_data.get(event_i)
    vn_filename = "particle_9999_vndata_diff_eta_-0.5_0.5.dat"
    vn_data = np.nan_to_num(eventGroup.get(vn_filename))
    dN_vector = calcualte_yield_and_meanpT(0.0, 3.0, vn_data)
    outdata[event_i] = {}
    outdata[event_i]["Nch"] = dN_vector[0]
    outdata[event_i]["mean_pT_ch"] = dN_vector[1]

    ini_nB = np.nan_to_num(eventGroup.get("nB_etas_distribution_N_72.dat"))
    FO_nB = np.nan_to_num(eventGroup.get("FO_nBvseta.dat"))
    outdata[event_i]["ini_nB"] = ini_nB
    outdata[event_i]["FO_nB"] = FO_nB
    dNdy_p = np.nan_to_num(eventGroup.get("particle_2212_dNdy_pT_0.2_3.dat"))
    dNdy_pbar = np.nan_to_num(eventGroup.get("particle_-2212_dNdy_pT_0.2_3.dat"))
    outdata[event_i]["proton"] = dNdy_p[:, :2]
    outdata[event_i]["anti-proton"] = dNdy_pbar[:, :2]

with open('nBarrs.pickle', 'wb') as pf:
    pickle.dump(outdata, pf)
