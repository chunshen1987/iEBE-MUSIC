#!/usr/bin/env python3

import h5py
import sys
import numpy as np

NORDER = 9
kinematicCutsDict = {
    "ALICE": {"pTmin": 0.2, "pTmax": 3.0, "etamin": -0.8, "etamax": 0.8},
    "CMS": {"pTmin": 0.3, "pTmax": 3.0, "etamin": -0.5, "etamax": 0.5},
    "ATLAS": {"pTmin": 0.5, "pTmax": 3.0, "etamin": -0.5, "etamax": 0.5},
}


def help_message():
    print("{0} database_file".format(sys.argv[0]))
    exit(0)


def calcualte_inte_Qn(pT_low, pT_high, data):
    """
        this function calculates the pT-integrated vn in a
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
    temp_vn_array = [N, meanpT]
    for iorder in range(1, NORDER+1):
        vn_real_event = data[:, 2*iorder]
        vn_imag_event = data[:, 2*iorder+1]
        vn_real_interp = np.interp(pT_inte_array, pT_event, vn_real_event)
        vn_imag_interp = np.interp(pT_inte_array, pT_event, vn_imag_event)
        Qn_real_inte = 2.*np.pi*np.sum(
                    vn_real_interp*dN_interp*pT_inte_array)*dpT
        Qn_imag_inte = 2.*np.pi*np.sum(
                    vn_imag_interp*dN_interp*pT_inte_array)*dpT
        temp_vn_array.append(Qn_real_inte)
        temp_vn_array.append(Qn_imag_inte)
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

outFileList = []
for expName in kinematicCutsDict:
    outFile = open("Qn_vectors_{}.txt".format(expName), "w")
    outFile.write("# Nch  <pT>(GeV)  Qn_real  Qn_imag (n=0,{})\n".format(NORDER))
    outFileList.append(outFile)
for ievent, event_i in enumerate(eventList):
    if ievent % 100 == 0:
        print("fetching event: {0} from the database {1} ...".format(
            event_i, database_file))
    eventGroup = h5_data.get(event_i)
    vn_filename = "particle_9999_vndata_diff_eta_-0.5_0.5.dat"
    vn_data = np.nan_to_num(eventGroup.get(vn_filename))
    dN_vector = calcualte_yield_and_meanpT(0.0, 3.0, vn_data)
    Qn_vector = []
    for exp_i, expName in enumerate(kinematicCutsDict):
        pTetacut = kinematicCutsDict[expName]
        vn_filename = 'particle_9999_vndata_diff_eta_{}_{}.dat'.format(
                                        pTetacut["etamin"], pTetacut["etamax"])
        vn_data = np.nan_to_num(eventGroup.get(vn_filename))
        Qn_vector += calcualte_inte_Qn(pTetacut["pTmin"], pTetacut["pTmax"],
                                       vn_data)

    # output Qn vectors
    output = dN_vector + Qn_vector
    for val in output:
        outFileList[exp_i].write("{:.4e}  ".format(val))
    outFileList[exp_i].write("\n")

h5_data.close()
for outFile in outFileList:
    outFile.close()
