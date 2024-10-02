#!/usr/bin/env python3

import h5py
import sys
import pickle
import numpy as np

NORDER = 4
pidList = [('ch', '9999'), ('pi+', '211'), ('pi-', '-211'), ('K+', '321'),
           ('K-', '-321'), ('p', '2212'), ('pbar', '-2212')]

weakFDFlag = True
weakString = ""
if weakFDFlag:
    weakString = "_weakFD"


def help_message():
    print("Usage: {0} database_file".format(sys.argv[0]))
    exit(0)


def calcualte_yield_and_meanpT_and_vn(pT_low, pT_high, data):
    """
        this function calculates the pT-integrated particle yield and mean pT
        given pT range (pT_low, pT_high) for every event in the data
    """
    npT = 50
    pT_inte_array = np.linspace(pT_low, pT_high, npT)
    dpT = pT_inte_array[1] - pT_inte_array[0]
    dN_event = data[:, 1]
    totalN_event = data[:, -1]
    pT_event = data[:, 0]
    dN_interp = np.exp(
        np.interp(pT_inte_array, pT_event, np.log(dN_event + 1e-30)))
    N = 2.*np.pi*np.sum(dN_interp*pT_inte_array)*dpT
    meanpT = (np.sum(dN_interp*pT_inte_array**2.)
              /np.sum(dN_interp*pT_inte_array))
    res_array = [N, meanpT]
    totalN_interp = np.exp(
        np.interp(pT_inte_array, pT_event, np.log(totalN_event + 1e-30)))
    totalN = np.sum(totalN_interp)*dpT/(pT_event[1] - pT_event[0])
    for iorder in range(1, NORDER + 1):
        vn_real_event = data[:, 2*iorder]
        vn_imag_event = data[:, 2*iorder + 1]
        vn_real_interp = np.interp(pT_inte_array, pT_event, vn_real_event)
        vn_imag_interp = np.interp(pT_inte_array, pT_event, vn_imag_event)
        Vn_real_inte = (np.sum(vn_real_interp*dN_interp*pT_inte_array)
                        /np.sum(dN_interp*pT_inte_array))
        Vn_imag_inte = (np.sum(vn_imag_interp*dN_interp*pT_inte_array)
                        /np.sum(dN_interp*pT_inte_array))
        res_array.append(Vn_real_inte + 1j*Vn_imag_inte)
    res_array.append(totalN)
    return res_array


try:
    database_file = str(sys.argv[1])
except IndexError:
    help_message()

h5_data = h5py.File(database_file, "r")
eventList = list(h5_data.keys())

outdata = {}

for ievent, event_i in enumerate(eventList):
    if ievent % 100 == 0:
        print("fetching event: {0} from the database {1} ...".format(
            event_i, database_file))
    eventGroup = h5_data.get(event_i)
    outdata[event_i] = {}
    vn_filename = f"particle_9999_vndata_diff_eta_-0.5_0.5{weakString}.dat"
    vn_data = np.nan_to_num(eventGroup.get(vn_filename))
    outdata[event_i]["pT"] = vn_data[:, 0]
    dN_vector = calcualte_yield_and_meanpT_and_vn(0.0, 3.0, vn_data)
    outdata[event_i]["Nch"] = dN_vector[0]
    outdata[event_i]["mean_pT_ch"] = dN_vector[1]

    # compute and save the reference Vn
    vn_filename = f"particle_9999_vndata_diff_eta_0.5_1{weakString}.dat"
    vn_data = np.nan_to_num(eventGroup.get(vn_filename))
    vn_vector = calcualte_yield_and_meanpT_and_vn(0.2, 2.0, vn_data)
    outdata[event_i]["Vn_pT_0p2_2_eta_0p5_1"] = vn_vector
    vn_filename = f"particle_9999_vndata_diff_eta_-1_-0.5{weakString}.dat"
    vn_data = np.nan_to_num(eventGroup.get(vn_filename))
    vn_vector = calcualte_yield_and_meanpT_and_vn(0.2, 2.0, vn_data)
    outdata[event_i]["Vn_pT_0p2_2_eta_-1_-0p5"] = vn_vector
    for pidName, pid in pidList:
        if pid == "9999":
            vn_filename = (
                f"particle_9999_vndata_diff_eta_-0.5_0.5{weakString}.dat")
        else:
            vn_filename = (
                f"particle_{pid}_vndata_diff_y_-0.5_0.5{weakString}.dat")
        vn_data = np.nan_to_num(eventGroup.get(vn_filename))
        # dN/(2pi pT dpT dy)
        outdata[event_i]["{}_sp".format(pidName)] = vn_data[:, 1]
        # v2(pT)
        outdata[event_i]["{}_v2pT".format(pidName)] = (vn_data[:, 4]
                                                       + 1j*vn_data[:, 5])
        # v3(pT)
        outdata[event_i]["{}_v3pT".format(pidName)] = (vn_data[:, 6]
                                                       + 1j*vn_data[:, 7])
        # v4(pT)
        outdata[event_i]["{}_v4pT".format(pidName)] = (vn_data[:, 8]
                                                       + 1j*vn_data[:, 9])

print("nev = {}".format(len(eventList)))
with open(f'pTdiffSpVn{weakString}.pickle', 'wb') as pf:
    pickle.dump(outdata, pf)

h5_data.close()
