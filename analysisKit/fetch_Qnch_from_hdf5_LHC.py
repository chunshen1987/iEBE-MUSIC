#!/usr/bin/env python3

import h5py
import sys
import pickle
import numpy as np

NORDER = 9
kinematicCutsDict = {
    "ALICE_eta_-0p4_0p4": {"pTmin": 0.2, "pTmax": 3,
                           "etamin": -0.4, "etamax": 0.4},
    "ALICE_eta_-0p8_-0p4": {"pTmin": 0.2, "pTmax": 3,
                            "etamin": -0.8, "etamax": -0.4},
    "ALICE_eta_0p4_0p8": {"pTmin": 0.2, "pTmax": 3,
                          "etamin": 0.4, "etamax": 0.8},
    "ALICE_eta_-0p8_0p8": {"pTmin": 0.2, "pTmax": 3,
                           "etamin": -0.8, "etamax": 0.8},
}

pidList = [('pi+', '211'), ('pi-', '-211'), ('K+', '321'), ('K-', '-321'),
           ('p', '2212'), ('pbar', '-2212')]

LHCetaRangeList = ['-0.4_0.4', '-0.5_0.5', '-0.8_-0.4', '-2.4_-0.5',
                   '-2.5_-0.5', '-3.7_-1.7', '-4.9_-3.1', '-5.1_-2.8',
                   '0.4_0.8', '0.5_2.4', '0.5_2.5', '1.7_3.7', '2.8_5.1',
                   '3.1_4.9', '2.8_5.1']

def help_message():
    print("Usage: {0} database_file".format(sys.argv[0]))
    exit(0)


def calcualte_inte_Vn_pT(pT_low, pT_high, data):
    """
        this function calculates the pT-integrated vn in a
        given pT range (pT_low, pT_high) for every event in the data
    """
    npT = 50
    pT_inte_array = np.linspace(pT_low, pT_high, npT)
    dpT = pT_inte_array[1] - pT_inte_array[0]
    dN_event = data[:, 1]
    totalN_event = data[:, -1]
    pT_event = data[:, 0]
    dN_interp = np.exp(np.interp(pT_inte_array, pT_event,
                                 np.log(dN_event+1e-30)))
    totalN_interp = np.exp(np.interp(pT_inte_array, pT_event,
                                     np.log(totalN_event+1e-30)))
    N = 2.*np.pi*np.sum(dN_interp*pT_inte_array)*dpT
    totalN = np.sum(totalN_interp)*dpT/(pT_event[1] - pT_event[0])
    meanpT = (np.sum(dN_interp*pT_inte_array**2.)
              / np.sum(dN_interp*pT_inte_array))
    temp_vn_array = [N, meanpT]
    for iorder in range(1, NORDER+1):
        vn_real_event = data[:, 2*iorder]
        vn_imag_event = data[:, 2*iorder+1]
        vn_real_interp = np.interp(pT_inte_array, pT_event, vn_real_event)
        vn_imag_interp = np.interp(pT_inte_array, pT_event, vn_imag_event)
        Vn_real_inte = (np.sum(vn_real_interp*dN_interp*pT_inte_array)
                        / np.sum(dN_interp*pT_inte_array))
        Vn_imag_inte = (np.sum(vn_imag_interp*dN_interp*pT_inte_array)
                        / np.sum(dN_interp*pT_inte_array))
        temp_vn_array.append(Vn_real_inte + 1j*Vn_imag_inte)
    temp_vn_array.append(totalN)
    return temp_vn_array


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
    return temp_vn_array


def calcualte_inte_Vn_pTeta(pTMin, pTMax, etaMin, etaMax, data, Nevents):
    """
        this function calculates the pT and eta-integrated vn in a
        given pT range (pTMin, pTMax) and eta range (etaMin, etaMax)
        for every event in the data
    """
    npT = 20; nEta = 71
    pTArr = np.linspace(0, 3.8, npT)
    dpT = pTArr[1] - pTArr[0]
    pTMinIdx = int((pTMin - 0.)/dpT)
    pTMaxIdx = min(npT, int((pTMax - 0.)/dpT) + 1)
    etaArr = np.linspace(-7, 7, nEta)
    deta = etaArr[1] - etaArr[0]
    etaMinIdx = max(0, int((etaMin - etaArr[0])/deta))
    etaMaxIdx = min(nEta, int((etaMax - etaArr[0])/deta) + 1)

    pT_event = data[:, 1].reshape(nEta, npT)
    dN_event = data[:, 2].reshape(nEta, npT)
    N = np.sum(dN_event[etaMinIdx:etaMaxIdx, pTMinIdx:pTMaxIdx])
    meanpT = (np.sum(pT_event[etaMinIdx:etaMaxIdx, pTMinIdx:pTMaxIdx]
                     *dN_event[etaMinIdx:etaMaxIdx, pTMinIdx:pTMaxIdx])/N)
    totalN = N*Nevents
    temp_vn_array = [N, meanpT]
    for iorder in range(1, NORDER+1):
        Qn_real_event = data[:, 2*iorder+2].reshape(nEta, npT)
        Qn_imag_event = data[:, 2*iorder+3].reshape(nEta, npT)
        Vn_real_inte = (
            np.sum(Qn_real_event[etaMinIdx:etaMaxIdx, pTMinIdx:pTMaxIdx])/N)
        Vn_imag_inte = (
            np.sum(Qn_imag_event[etaMinIdx:etaMaxIdx, pTMinIdx:pTMaxIdx])/N)
        temp_vn_array.append(Vn_real_inte + 1j*Vn_imag_inte)
    temp_vn_array.append(totalN)
    return temp_vn_array


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
    meanpT = (np.sum(dN_interp*pT_inte_array**2.)
              / np.sum(dN_interp*pT_inte_array))
    res_array = [N, meanpT]
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
    vn_filename = "particle_9999_vndata_diff_eta_-0.5_0.5.dat"
    vn_data = np.nan_to_num(eventGroup.get(vn_filename))
    ecc_filename = "eccentricities_evo_ed_tau_0.4.dat"
    eccn_data = np.nan_to_num(eventGroup.get(ecc_filename))
    dN_vector = calcualte_yield_and_meanpT(0.0, 3.0, vn_data)
    outdata[event_i]["Nch"] = dN_vector[0]
    outdata[event_i]["mean_pT_ch"] = dN_vector[1]
    outdata[event_i]["ecc_n"] = eccn_data[2:]
    fileList = list(eventGroup.keys())
    for filename in fileList:
        if "NgluonEstimators" in filename:
            data = np.nan_to_num(eventGroup.get(filename))
            outdata[event_i]['NgluonEst'] = data[0]
            continue
    for pidName, pid in pidList:
        vn_filename = "particle_{}_vndata_diff_y_-0.5_0.5.dat".format(pid)
        vn_data = np.nan_to_num(eventGroup.get(vn_filename))
        dN_vector = calcualte_yield_and_meanpT(0.0, 3.0, vn_data)
        outdata[event_i]["{}_dNdy_meanpT".format(pidName)] = dN_vector
    vn_filename = 'particle_9999_pTeta_distribution.dat'
    vnInte_filename = 'particle_9999_vndata_eta_-0.5_0.5.dat'
    for exp_i, expName in enumerate(kinematicCutsDict):
        pTetacut = kinematicCutsDict[expName]
        vn_data = np.nan_to_num(eventGroup.get(vn_filename))
        vnInte_data = np.nan_to_num(eventGroup.get(vnInte_filename))
        N_hadronic_events = vnInte_data[-1, 2]
        Vn_vector = calcualte_inte_Vn_pTeta(
                    pTetacut['pTmin'], pTetacut['pTmax'],
                    pTetacut['etamin'], pTetacut['etamax'],
                    vn_data, N_hadronic_events)
        outdata[event_i][expName] = np.array(Vn_vector)

print("nev = {}".format(len(eventList)))
with open('QnVectors.pickle', 'wb') as pf:
    pickle.dump(outdata, pf)

h5_data.close()
