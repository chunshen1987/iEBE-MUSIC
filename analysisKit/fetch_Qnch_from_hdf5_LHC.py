#!/usr/bin/env python3

import h5py
import sys
import pickle
import numpy as np

NORDER = 9
kinematicCutsDict = {
    "ALICE_eta_-0p4_0p4": {"pTmin": 0.2, "pTmax": 3, "etamin": -0.4, "etamax": 0.4},
    "ALICE_eta_-0p8_-0p4": {"pTmin": 0.2, "pTmax": 3, "etamin": -0.8, "etamax": -0.4},
    "ALICE_eta_0p4_0p8": {"pTmin": 0.2, "pTmax": 3, "etamin": 0.4, "etamax": 0.8},
    "ALICE_eta_-0p8_0p8": {"pTmin": 0.2, "pTmax": 3, "etamin": -0.8, "etamax": 0.8},
}

pidList = [('pi+', '211'), ('pi-', '-211'), ('K+', '321'), ('K-', '-321'),
           ('p', '2212'), ('pbar', '-2212')]


def help_message():
    print("{0} database_file".format(sys.argv[0]))
    exit(0)


def calcualte_inte_Vn_eta(etaMin, etaMax, data):
    """
        this function calculates the eta-integrated vn in a
        given eta range (etaMin, etaMax) for every event in the data
    """
    nEta = 50
    eta_inte_array = np.linspace(etaMin, etaMax, nEta)
    deta = eta_inte_array[1] - eta_inte_array[0]
    dN_event = data[:, 1]
    ET_event = data[:, 2]
    totalN_event = data[:, -1]
    eta_event = data[:, 0]
    dN_interp = np.exp(np.interp(eta_inte_array, eta_event,
                                 np.log(dN_event+1e-30)))
    totalN_interp = np.exp(np.interp(eta_inte_array, eta_event,
                                     np.log(totalN_event+1e-30)))
    ET_interp = np.exp(np.interp(eta_inte_array, eta_event,
                                 np.log(ET_event + 1e-30)))
    N = np.sum(dN_interp)*deta
    totalN = np.sum(totalN_interp)*deta/(eta_event[1] - eta_event[0])
    ET = np.sum(ET_interp)*deta
    temp_vn_array = [N, ET]
    for iorder in range(1, NORDER+1):
        vn_real_event = data[:, 2*iorder+1]
        vn_imag_event = data[:, 2*iorder+2]
        vn_real_interp = np.interp(eta_inte_array, eta_event, vn_real_event)
        vn_imag_interp = np.interp(eta_inte_array, eta_event, vn_imag_event)
        Vn_real_inte = np.sum(vn_real_interp*dN_interp)/np.sum(dN_interp)
        Vn_imag_inte = np.sum(vn_imag_interp*dN_interp)/np.sum(dN_interp)
        temp_vn_array.append(Vn_real_inte + 1j*Vn_imag_inte)
    temp_vn_array.append(totalN)
    return(temp_vn_array)


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
    outdata[event_i] = {}
    vn_filename = "particle_9999_vndata_diff_eta_-0.5_0.5.dat"
    vn_data = np.nan_to_num(eventGroup.get(vn_filename))
    ecc_filename = "eccentricities_evo_ed_tau_0.4.dat"
    eccn_data = np.nan_to_num(eventGroup.get(ecc_filename))
    dN_vector = calcualte_yield_and_meanpT(0.0, 3.0, vn_data)
    outdata[event_i]["Nch"] = dN_vector[0]
    outdata[event_i]["mean_pT_ch"] = dN_vector[1]
    outdata[event_i]["ecc_n"] = eccn_data[2:]
    for pidName, pid in pidList:
        vn_filename = "particle_{}_vndata_diff_y_-0.5_0.5.dat".format(pid)
        vn_data = np.nan_to_num(eventGroup.get(vn_filename))
        dN_vector = calcualte_yield_and_meanpT(0.0, 3.0, vn_data)
        outdata[event_i]["{}_dNdy_meanpT".format(pidName)] = dN_vector
    vn_filename = 'particle_9999_pTeta_distribution.dat'
    for exp_i, expName in enumerate(kinematicCutsDict):
        pTetacut = kinematicCutsDict[expName]
        pTlabel = "{}_{}".format(pTetacut['pTmin'], pTetacut['pTmax'])
        if pTlabel in ['0.2_3', '0.3_3', '0.5_3']:
            vn_filename = 'particle_9999_dNdeta_pT_{}.dat'.format(pTlabel)
        vn_data = np.nan_to_num(eventGroup.get(vn_filename))
        Vn_vector = calcualte_inte_Vn_eta(
                        pTetacut['etamin'], pTetacut['etamax'], vn_data)
        outdata[event_i][expName] = np.array(Vn_vector)

print("nev = {}".format(len(eventList)))
with open('QnVectors.pickle', 'wb') as pf:
    pickle.dump(outdata, pf)

h5_data.close()
