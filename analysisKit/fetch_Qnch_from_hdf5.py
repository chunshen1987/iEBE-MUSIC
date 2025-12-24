#!/usr/bin/env python3

import h5py
import sys
import pickle
import numpy as np
from scipy.interpolate import RegularGridInterpolator

EPS = 1e-15

NORDER = 9

expFlag = "ALICE"
initialFlag = False
pTdiffFlag = True
etadiffFlag = True
photonFlag = False

weakFDFlag = False
weakString = ""
if weakFDFlag:
    weakString = "_weakFD"


kinematicCutsDict_STAR = {
    "STAR_eta_-0p5_0p5_pT_0p2_2": {
        "pTmin": 0.2,
        "pTmax": 2,
        "etamin": -0.5,
        "etamax": 0.5
    },
    "STAR_eta_-1_-0p5_pT_0p2_2": {
        "pTmin": 0.2,
        "pTmax": 2,
        "etamin": -1,
        "etamax": -0.5
    },
    "STAR_eta_0p5_1_pT_0p2_2": {
        "pTmin": 0.2,
        "pTmax": 2,
        "etamin": 0.5,
        "etamax": 1
    },
    "STAR_eta_-0p5_0_pT_0p2_2": {
        "pTmin": 0.2,
        "pTmax": 2,
        "etamin": -0.5,
        "etamax": 0
    },
    "STAR_eta_0_0p5_pT_0p2_2": {
        "pTmin": 0.2,
        "pTmax": 2,
        "etamin": 0,
        "etamax": 0.5
    },
}

kinematicCutsDict_ALICE = {
    "ALICE_eta_-0p5_0p5_pT_0_4": {
        "pTmin": 0.0,
        "pTmax": 4,
        "etamin": -0.5,
        "etamax": 0.5
    },
    "ALICE_V0A_eta_2p8_5p1_pT_0_4": {
        "pTmin": 0,
        "pTmax": 4,
        "etamin": 2.8,
        "etamax": 5.1
    },
    "ALICE_V0C_eta_-3p7_-1p7_pT_0_4": {
        "pTmin": 0,
        "pTmax": 4,
        "etamin": -3.7,
        "etamax": -1.7
    },
    "ALICE_eta_-0p5_0p5_pT_0p2_3": {
        "pTmin": 0.2,
        "pTmax": 3,
        "etamin": -0.5,
        "etamax": 0.5
    },
    "ALICE_eta_-0p8_0p8_pT_0p2_4": {
        "pTmin": 0.2,
        "pTmax": 4,
        "etamin": -0.8,
        "etamax": 0.8
    },
    "ALICE_eta_-0p8_0_pT_0p15_2": {
        "pTmin": 0.15,
        "pTmax": 2,
        "etamin": -0.8,
        "etamax": 0
    },
    "ALICE_eta_0_0p8_pT_0p15_2": {
        "pTmin": 0.15,
        "pTmax": 2,
        "etamin": 0,
        "etamax": 0.8
    },
    "ALICE_eta_-0p4_0p4_pT_0p2_3": {
        "pTmin": 0.2,
        "pTmax": 3,
        "etamin": -0.4,
        "etamax": 0.4
    },
    "ALICE_eta_-0p8_-0p5_pT_0p2_4": {
        "pTmin": 0.2,
        "pTmax": 4,
        "etamin": -0.8,
        "etamax": -0.5
    },
    "ALICE_eta_0p5_0p8_pT_0p2_4": {
        "pTmin": 0.2,
        "pTmax": 4,
        "etamin": 0.5,
        "etamax": 0.8
    },
}

kinematicCutsDict = kinematicCutsDict_ALICE
if expFlag == "ALICE":
    kinematicCutsDict = kinematicCutsDict_ALICE
elif expFlag == "STAR":
    kinematicCutsDict = kinematicCutsDict_STAR

pidList = [('ch', '9999'), ('pi+', '211'), ('pi-', '-211'), ('K+', '321'),
           ('K-', '-321'), ('p', '2212'), ('pbar', '-2212')]

photonList = [
    'QGP_2to2_total',
    'QGP_AMYcollinear',
    'HG_2to2_meson_total',
    'HG_omega',
    'HG_rho_spectralfun',
    'HG_pipi_bremsstrahlung',
]

LHCetaRangeList = [
    '-0.4_0.4', '-0.5_0.5', '-0.8_-0.4', '-2.4_-0.5', '-2.5_-0.5', '-3.7_-1.7',
    '-4.9_-3.1', '-5.1_-2.8', '0.4_0.8', '0.5_2.4', '0.5_2.5', '1.7_3.7',
    '2.8_5.1', '3.1_4.9', '2.8_5.1'
]
RHICetaRangeList = ['-0.5_0.5', '-1_-0.5', '-3.9_-3.1', '0.5_1', '3.1_3.9']


def help_message():
    print("Usage: {0} database_file".format(sys.argv[0]))
    exit(0)


def get3DGlauberData(h5Event):
    """
        this function gets the 3D Glauber data from the hdf5 event
    """
    data = []
    for fileName in h5Event.keys():
        if "strings_" in fileName:
            data = h5Event.get(fileName).attrs.get("header")
            break
    data = data.decode("utf-8").split(" ")
    b = float(data[3])
    Npart = int(data[7])
    Ncoll = int(data[10])
    Nstrings = int(data[13])
    Etot = float(data[16])
    Pztot = float(data[20])
    randomSeed = int(data[24])
    return [b, Npart, Ncoll, Nstrings, Etot, Pztot, randomSeed]


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
    dN_interp = np.exp(
        np.interp(pT_inte_array, pT_event, np.log(dN_event + 1e-30)))
    totalN_interp = np.exp(
        np.interp(pT_inte_array, pT_event, np.log(totalN_event + 1e-30)))
    N = 2.*np.pi*np.sum(dN_interp*pT_inte_array)*dpT
    totalN = np.sum(totalN_interp)*dpT/(pT_event[1] - pT_event[0])
    meanpT = (np.sum(dN_interp*pT_inte_array**2.)
              /np.sum(dN_interp*pT_inte_array))
    temp_vn_array = [N, meanpT]
    for iorder in range(1, NORDER + 1):
        vn_real_event = data[:, 2*iorder]
        vn_imag_event = data[:, 2*iorder + 1]
        vn_real_interp = np.interp(pT_inte_array, pT_event, vn_real_event)
        vn_imag_interp = np.interp(pT_inte_array, pT_event, vn_imag_event)
        Vn_real_inte = (np.sum(vn_real_interp*dN_interp*pT_inte_array)
                        /np.sum(dN_interp*pT_inte_array))
        Vn_imag_inte = (np.sum(vn_imag_interp*dN_interp*pT_inte_array)
                        /np.sum(dN_interp*pT_inte_array))
        temp_vn_array.append(Vn_real_inte + 1j*Vn_imag_inte)
    temp_vn_array.append(totalN)
    return temp_vn_array


def calcualte_inte_Vn_eta(etaMin, etaMax, data, vnFlag=True):
    """
        this function calculates the eta-integrated vn in a
        given eta range (etaMin, etaMax) for every event in the data
    """
    nEta = 50
    eta_inte_array = np.linspace(etaMin, etaMax, nEta)
    deta = eta_inte_array[1] - eta_inte_array[0]
    dN_event = data[:, 1]
    ET_event = data[:, -2]
    totalN_event = data[:, -1]
    eta_event = data[:, 0]
    dN_interp = np.exp(
        np.interp(eta_inte_array, eta_event, np.log(dN_event + 1e-30)))
    totalN_interp = np.exp(
        np.interp(eta_inte_array, eta_event, np.log(totalN_event + 1e-30)))
    ET_interp = np.exp(
        np.interp(eta_inte_array, eta_event, np.log(ET_event + 1e-30)))
    N = np.sum(dN_interp)*deta
    totalN = np.sum(totalN_interp)*deta/(eta_event[1] - eta_event[0])
    ET = np.sum(ET_interp)*deta
    temp_vn_array = [N, ET]
    if vnFlag:
        for iorder in range(1, NORDER + 1):
            vn_real_event = data[:, 2*iorder + 1]
            vn_imag_event = data[:, 2*iorder + 2]
            vn_real_interp = np.interp(eta_inte_array, eta_event, vn_real_event)
            vn_imag_interp = np.interp(eta_inte_array, eta_event, vn_imag_event)
            Vn_real_inte = np.sum(vn_real_interp*dN_interp)/np.sum(dN_interp)
            Vn_imag_inte = np.sum(vn_imag_interp*dN_interp)/np.sum(dN_interp)
            temp_vn_array.append(Vn_real_inte + 1j*Vn_imag_inte)
        temp_vn_array.append(totalN)
    return temp_vn_array


def calcualte_inte_Vn_pTeta(pTMin: float, pTMax: float, etaMin: float,
                            etaMax: float, data: np.ndarray,
                            etabin: np.ndarray, pTbin: np.ndarray,
                            Nevents: int, vnOrder:int):
    """
        this function calculates the pT and eta-integrated vn in a
        given pT range (pTMin, pTMax) and eta range (etaMin, etaMax)
        for every event in the data
    """
    npT = len(pTbin); dpT = pTbin[1] - pTbin[0]
    nEta = len(etabin); dEta = etabin[1] - etabin[0]
    weights = np.zeros([nEta, npT])

    pTIdx0 = np.where(pTbin > pTMin)[0][0]
    pTIdx1 = np.where(pTbin < pTMax)[0][-1]
    etaIdx0 = np.where(etabin > etaMin)[0][0]
    etaIdx1 = np.where(etabin < etaMax)[0][-1]
    weights[etaIdx0:etaIdx1+1, pTIdx0:pTIdx1+1] = 1.

    pTminFrac = (pTbin[pTIdx0] - pTMin) / dpT
    pTmaxFrac = (pTMax - pTbin[pTIdx1]) / dpT
    etaminFrac = (etabin[etaIdx0] - etaMin) / dEta
    etamaxFrac = (etaMax - etabin[etaIdx1]) / dEta
    if pTminFrac <= 0.5:
        weights[etaIdx0:etaIdx1+1, pTIdx0] *= (0.5 + pTminFrac)
    elif pTminFrac < 1.0 and pTIdx0 > 0:
        weights[etaIdx0:etaIdx1+1, pTIdx0 - 1] = pTminFrac - 0.5

    if etaminFrac <= 0.5:
        weights[etaIdx0, pTIdx0:pTIdx1+1] *= (0.5 + etaminFrac)
    elif etaminFrac < 1.0 and etaIdx0 > 0:
        weights[etaIdx0 - 1, pTIdx0:pTIdx1+1] = etaminFrac - 0.5
        if pTminFrac < 1.0 and pTIdx0 > 0:
            weights[etaIdx0 - 1, pTIdx0 - 1] = ((pTminFrac - 0.5)
                                                * (etaminFrac - 0.5))

    if pTmaxFrac <= 0.5:
        weights[etaIdx0:etaIdx1+1, pTIdx1] *= (0.5 + pTmaxFrac)
    elif pTmaxFrac < 1.0 and pTIdx1 < npT - 1:
        weights[etaIdx0:etaIdx1+1, pTIdx1 + 1] = pTmaxFrac - 0.5

    if etamaxFrac <= 0.5:
        weights[etaIdx1, pTIdx0:pTIdx1+1] *= (0.5 + etamaxFrac)
    elif etamaxFrac < 1.0 and etaIdx1 < nEta - 1:
        weights[etaIdx1 + 1, pTIdx0:pTIdx1+1] = etamaxFrac - 0.5
        if pTmaxFrac < 1.0 and pTIdx1 < npT - 1:
            weights[etaIdx1 + 1, pTIdx1 + 1] = ((pTmaxFrac - 0.5)
                                                * (etamaxFrac - 0.5))

    dN_event = data[:, 2].reshape(nEta, npT)
    ET_event = data[:, 3].reshape(nEta, npT)
    N = np.sum(dN_event * weights) + EPS
    ET = np.sum(ET_event * weights) + EPS
    pT_event = data[:, 1].reshape(nEta, npT)
    meanpT = np.sum(pT_event * weights * dN_event) / N

    totalN = N*Nevents
    temp_vn_array = [N, meanpT]
    for iorder in range(1, vnOrder + 1):
        Qn_real_event = data[:, 2*iorder + 2].reshape(nEta, npT)
        Qn_imag_event = data[:, 2*iorder + 3].reshape(nEta, npT)
        Vn_real_inte = np.sum(Qn_real_event * weights) / N
        Vn_imag_inte = np.sum(Qn_imag_event * weights) / N
        temp_vn_array.append(Vn_real_inte + 1j*Vn_imag_inte)
    temp_vn_array.append(ET)
    temp_vn_array.append(totalN)
    return temp_vn_array


def calcualte_inte_Vneta_pTeta(pTMin: float, pTMax: float, data: np.ndarray,
                               etabin: np.ndarray, pTbin: np.ndarray,
                               Nevents: int, weightType: int):
    """
        this function calculates the pT-integrated vn(eta) in a
        given pT range (pTMin, pTMax) for every event in the data
    """
    npT = len(pTbin); dpT = pTbin[1] - pTbin[0]
    nEta = len(etabin); dEta = etabin[1] - etabin[0]
    weights = np.zeros([nEta, npT])

    pTIdx0 = np.where(pTbin > pTMin)[0][0]
    pTIdx1 = np.where(pTbin < pTMax)[0][-1]
    weights[:, pTIdx0:pTIdx1+1] = 1.

    pTminFrac = (pTbin[pTIdx0] - pTMin) / dpT
    pTmaxFrac = (pTMax - pTbin[pTIdx1]) / dpT

    if pTminFrac <= 0.5:
        weights[:, pTIdx0] *= (0.5 + pTminFrac)
    elif pTminFrac < 1.0 and pTIdx0 > 0:
        weights[:, pTIdx0 - 1] = pTminFrac - 0.5

    if pTmaxFrac <= 0.5:
        weights[:, pTIdx1] *= (0.5 + pTmaxFrac)
    elif pTmaxFrac < 1.0 and pTIdx1 < npT - 1:
        weights[:, pTIdx1 + 1] = pTmaxFrac - 0.5

    pT_event = data[:, 1].reshape(nEta, npT)
    dN_event = np.real(data[:, 2].reshape(nEta, npT))
    ET_event = np.real(data[:, 3].reshape(nEta, npT))

    N = np.sum(dN_event * weights, axis=1) + EPS
    ET = np.sum(ET_event * weights, axis=1)
    meanpT = np.sum(pT_event * weights * dN_event, axis=1) / N
    totalN = N*Nevents
    temp_vn_array = [N/dEta, meanpT]        # dN/deta, <pT>(eta)
    for iorder in range(1, NORDER + 1):
        Qn_real_event = data[:, 2*iorder + 2].reshape(nEta, npT)
        Qn_imag_event = data[:, 2*iorder + 3].reshape(nEta, npT)
        if weightType == 1:
            Vn_real_inte = (np.sum(Qn_real_event * pT_event * weights, axis=1)
                            / N)
            Vn_imag_inte = (np.sum(Qn_imag_event * pT_event * weights, axis=1)
                            / N)
        else:
            Vn_real_inte = np.sum(Qn_real_event * weights, axis=1) / N
            Vn_imag_inte = np.sum(Qn_imag_event * weights, axis=1) / N
        temp_vn_array.append(Vn_real_inte + 1j*Vn_imag_inte)    # Vn(eta)
    temp_vn_array.append(ET/dEta)            # dET/deta
    temp_vn_array.append(totalN)             # totalN
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
    dN_interp = np.exp(
        np.interp(pT_inte_array, pT_event, np.log(dN_event + 1e-30)))
    N = 2.*np.pi*np.sum(dN_interp*pT_inte_array)*dpT
    meanpT = (np.sum(dN_interp*pT_inte_array**2.)
              /np.sum(dN_interp*pT_inte_array))
    res_array = [N, meanpT]
    return res_array


try:
    database_file = str(sys.argv[1])
except IndexError:
    help_message()

h5_data = h5py.File(database_file, "r")
eventList = list(h5_data.keys())

outdata = {}
outdata["global"] = {}

for ievent, event_i in enumerate(eventList):
    if ievent % 100 == 0:
        print("fetching event: {0} from the database {1} ...".format(
            event_i, database_file))
    eventGroup = h5_data.get(event_i)
    vn_filename = f"particle_9999_vndata_eta_-0.5_0.5{weakString}.dat"
    vn_data = np.nan_to_num(eventGroup.get(vn_filename))
    if vn_data.ndim == 2:
        outdata[event_i] = {}
        outdata[event_i]["Nch"] = vn_data[0, 1]
        outdata[event_i]["mean_pT_ch"] = vn_data[0, 3]
    else:
        print(f"skip {event_i} ...")
        print("vn_data shape: ", vn_data.shape)
        continue

    # identified particle yields and mean pT
    for pidName, pid in pidList[1:]:
        vn_filename = f"particle_{pid}_vndata_y_-0.5_0.5{weakString}.dat"
        vn_data = np.nan_to_num(eventGroup.get(vn_filename))
        outdata[event_i]["{}_dNdy_meanpT".format(pidName)] = (
                                            [vn_data[0, 1], vn_data[0, 3]])

    if initialFlag:
        # initial eccentricity
        ecc_filename = "eccentricities_evo_ed_tau_INITIAL.dat"
        eccn_data = np.nan_to_num(eventGroup.get(ecc_filename))
        if eccn_data.ndim == 1:
            outdata[event_i]["ecc_n"] = eccn_data[2:]
        else:
            # compute ecc_n in the midrapidity
            neta, norder = eccn_data[:, 2:].shape
            eccn = np.zeros([norder])
            etaArr = np.linspace(-0.5, 0.5, 11)
            edInterp = np.interp(etaArr, eccn_data[:, 0], eccn_data[:, 1])
            ednum = np.sum(edInterp)
            for iorder in range(norder):
                eccnInterp = np.interp(etaArr, eccn_data[:, 0],
                                       eccn_data[:, 1]*eccn_data[:, iorder+2])
                eccn[iorder] = np.sum(eccnInterp)/ednum
            outdata[event_i]["ecc_n"] = eccn

    vnInte_filename = f'particle_9999_vndata_eta_-0.5_0.5{weakString}.dat'
    vnInte_data = np.nan_to_num(eventGroup.get(vnInte_filename))
    N_hadronic_events = vnInte_data[-1, 2]
    if ievent == 0:
        vneta_filename = f"particle_99999_dNdeta_pT_0_4{weakString}.dat"
        vneta_data = np.nan_to_num(eventGroup.get(vneta_filename))
        nEta = len(vneta_data[:, 0])
        etaMax = round(vneta_data[-1, 0], 1)
        etaMin = round(vneta_data[0, 0], 1)
        outdata["global"]["etaArr"] = np.linspace(etaMin, etaMax, nEta)
        vnpT_filename = f"particle_99999_vndata_diff_eta_-0.5_0.5{weakString}.dat"
        vnpT_data = np.nan_to_num(eventGroup.get(vnpT_filename))
        npT = len(vnpT_data[:, 0])
        pTMax = round(vnpT_data[-1, 0], 1)
        pTMin = round(vnpT_data[0, 0], 1)
        outdata["global"]["pTArr"] = np.linspace(pTMin, pTMax, npT)

    # compute dET/deta at mid-rapidity
    vn_filename = f'particle_99999_pTeta_distribution{weakString}.dat'
    vn_data = np.nan_to_num(eventGroup.get(vn_filename))
    dN_vector = calcualte_inte_Vn_pTeta(0, 4, -0.5, 0.5, vn_data,
                                        outdata["global"]['etaArr'],
                                        outdata["global"]['pTArr'],
                                        N_hadronic_events, 4)
    outdata[event_i]["ET"] = np.real(dN_vector[-2])


    # charged hadron vn with different kinematic cuts
    vn_filename = f'particle_9999_pTeta_distribution{weakString}.dat'
    vn_data = np.nan_to_num(eventGroup.get(vn_filename))
    for exp_i, expName in enumerate(kinematicCutsDict):
        pTetacut = kinematicCutsDict[expName]
        Vn_vector = calcualte_inte_Vn_pTeta(pTetacut['pTmin'],
                                            pTetacut['pTmax'],
                                            pTetacut['etamin'],
                                            pTetacut['etamax'], vn_data,
                                            outdata["global"]['etaArr'],
                                            outdata["global"]['pTArr'],
                                            N_hadronic_events, NORDER)
        outdata[event_i][expName] = np.array(Vn_vector)

    if pTdiffFlag:
        # pT-differential spectra and vn
        for pidName, pid in pidList:
            if pid == "9999":
                vn_filename = (
                    f"particle_9999_vndata_diff_eta_-0.5_0.5{weakString}.dat")
            else:
                vn_filename = (
                    f"particle_{pid}_vndata_diff_y_-0.5_0.5{weakString}.dat")
            vn_data = np.nan_to_num(eventGroup.get(vn_filename))
            pTdiffData = [vn_data[:, 1]]
            for iOrder in range(1, 5):
                pTdiffData.append(vn_data[:, 2*iOrder]
                                  + 1j*vn_data[:, 2*iOrder + 1])
            outdata[event_i][f"{pidName}_pTArr"] = np.array(pTdiffData)

    if photonFlag:
        eventData = get3DGlauberData(eventGroup)
        outdata[event_i]["Ncoll"] = eventData[2]
        outdata[event_i]["Npart"] = eventData[1]
        outdata[event_i]["b"] = eventData[0]

        photonFullres = []
        for ichannel, channelName in enumerate(photonList):
            vn_filename = f"{channelName}_Spvn_tot_ypTdiff.dat"
            vn_data = np.nan_to_num(eventGroup.get(vn_filename))
            if ichannel == 0:
                photonFullres = vn_data
                photonFullres[:, 3:] = (vn_data[:, 3:]
                                        *vn_data[:, 2].reshape(-1, 1))
            else:
                photonFullres[:, 2] += vn_data[:, 2]
                photonFullres[:, 3:] += (vn_data[:, 3:]
                                         *vn_data[:, 2].reshape(-1, 1))
        photonFullres[:, 3:] /= photonFullres[:, 2].reshape(-1, 1)
        outdata[event_i]["photon_ypTdiff"] = photonFullres[:, 2:]
        outdata[event_i]["photon_pTArr"] = np.unique(photonFullres[:, 1])
        outdata[event_i]["photon_yArr"] = np.unique(photonFullres[:, 0])

        dileptonFileName = "Dilepton_QGPNLO_Spvn_eq_MInv.dat"
        outdata[event_i]["dilepton_MInv"] = (np.nan_to_num(
            eventGroup.get(dileptonFileName)))
        dileptonFileName = "Dilepton_QGPNLO_Spvn_eq_MInvpTdiff.dat"
        outdata[event_i]["dilepton_MInvpTdiff"] = (np.nan_to_num(
            eventGroup.get(dileptonFileName)))

    if etadiffFlag:
        # eta-differential spectra and vn
        vn_filename = f"particle_9999_dNdeta_pT_0.2_3{weakString}.dat"
        vn_data = np.nan_to_num(eventGroup.get(vn_filename))
        outdata[event_i]["dNch/deta"] = vn_data[:, 1]

        vn_filename = f"particle_99999_dNdeta_pT_0_4{weakString}.dat"
        vn_data = np.nan_to_num(eventGroup.get(vn_filename))
        outdata[event_i]["dET/deta"] = vn_data[:, -2]

        vn_filename = f'particle_9999_pTeta_distribution{weakString}.dat'
        vn_data = np.nan_to_num(eventGroup.get(vn_filename))
        if expFlag == "STAR":
            # for longitudinal derrelation
            Vn_vector = calcualte_inte_Vneta_pTeta(0.4, 4.0, vn_data,
                                                   outdata['global']["etaArr"],
                                                   outdata['global']["pTArr"],
                                                   N_hadronic_events, 0)
            outdata[event_i]["chVneta_pT_0p4_4"] = np.array(Vn_vector)
            # for vn(eta)
            Vn_vector = calcualte_inte_Vneta_pTeta(0.1, 4.0, vn_data,
                                                   outdata['global']["etaArr"],
                                                   outdata['global']["pTArr"],
                                                   N_hadronic_events, 0)
            outdata[event_i]["chVneta_pT_0p1_4"] = np.array(Vn_vector)
            Vn_vector = calcualte_inte_Vneta_pTeta(0.15, 2.0, vn_data,
                                                   outdata['global']["etaArr"],
                                                   outdata['global']["pTArr"],
                                                   N_hadronic_events, 0)
            outdata[event_i]["chVneta_pT_0p15_2"] = np.array(Vn_vector)
            Vn_vector = calcualte_inte_Vneta_pTeta(0.15, 2.0, vn_data,
                                                   outdata['global']["etaArr"],
                                                   outdata['global']["pTArr"],
                                                   N_hadronic_events, 1)
            outdata[event_i]["chVneta_pTw_pT_0p15_2"] = np.array(Vn_vector)
            Vn_vector = calcualte_inte_Vneta_pTeta(0.2, 2.0, vn_data,
                                                   outdata['global']["etaArr"],
                                                   outdata['global']["pTArr"],
                                                   N_hadronic_events, 0)
            outdata[event_i]["chVneta_pT_0p2_2"] = np.array(Vn_vector)

        # for vn(eta)
        Vn_vector = calcualte_inte_Vneta_pTeta(0.2, 3.0, vn_data,
                                               outdata['global']["etaArr"],
                                               outdata['global']["pTArr"],
                                               N_hadronic_events, 0)
        outdata[event_i]["chVneta_pT_0p2_3"] = np.array(Vn_vector)

print("nev = {}".format(len(eventList)))
with open(f'QnVectors{weakString}.pickle', 'wb') as pf:
    pickle.dump(outdata, pf)

h5_data.close()
