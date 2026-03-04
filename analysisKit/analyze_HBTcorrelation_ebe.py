#! /usr/bin/env python
"""
     This script performs event averaging for the HBT correlation function
     calculated from event-by-event simulations
"""

from sys import argv, exit
from os import path
import numpy as np
import h5py
import pickle
from scipy.optimize import curve_fit

invFlag = False
HBARC = 0.19733
EPS = 1e-15
Nq = 31

KT_values = ['0.15_0.25', '0.25_0.35', '0.35_0.45', '0.45_0.55']
KTMid = [0.2, 0.3, 0.4, 0.5]
HBTFileList = ['HBT_correlation_function_KT']
if invFlag:
    HBTFileList += ['HBT_correlation_function_inv_KT']

qCutMaxList = [0.06, 0.08, 0.1, 0.12, 0.14]
qCutMinList = [0., 0.01, 0.02, 0.03]

qArr = np.linspace(-0.15, 0.15, Nq)
# 3D ordering is long, out, side
q3DArr = []
for iLong in range(Nq):
    for iOut in range(Nq):
        for iSide in range(Nq):
            qMag = np.sqrt(qArr[iLong]**2 + qArr[iOut]**2 + qArr[iSide]**2)
            q3DArr.append([qArr[iLong], qArr[iOut], qArr[iSide], qMag])
q3DArr = np.array(q3DArr)

header = ("# q_cut[GeV]  lambda  lambda_err  R_out[fm]  R_out_err[fm]  "
          + "R_side[fm]  R_side_err[fm]  R_long [fm]  R_long_err[fm]  "
          + "R_os[fm]  R_os_err[fm]  R_ol[fm]  R_ol_err[fm]")

try:
    databaseFile = path.abspath(argv[1])
except IndexError:
    print("Usage: {} databaseFile.h5".format(argv[0]))
    exit(1)


def check_an_event_is_good(h5_event):
    """This function checks the given event contains all required files"""
    required_files_list = [
        'particle_9999_vndata_eta_-0.5_0.5.dat',
        'HBT_correlation_function_KT_0.15_0.25.dat',
    ]
    if invFlag:
        required_files_list += ['HBT_correlation_function_inv_KT_0.15_0.25.dat']
    event_file_list = list(h5_event.keys())
    for ifile in required_files_list:
        if ifile not in event_file_list:
            print("event {} is bad, missing {} ...".format(
                h5_event.name, ifile))
            return False
    return True


def gaussian_1d(q_arr, lambda_, R_inv):
    R_inv /= HBARC  # [1/GeV]
    gauss = lambda_*np.exp(-(R_inv*q_arr)**2)

    return gauss


def gaussian_3d(q_arr, lambda_, R_out, R_side, R_long, R_os, R_ol):
    """ the fit function is according to arXiv: 1403.4972v1 """
    q_long = q_arr[:, 0]
    q_out = q_arr[:, 1]
    q_side = q_arr[:, 2]

    R_out /= HBARC  # [1/GeV]
    R_side /= HBARC
    R_long /= HBARC
    R_os /= HBARC
    R_ol /= HBARC
    gauss = (lambda_*np.exp(-(
        (R_out*q_out)**2 + (R_side*q_side)**2 + (R_long*q_long)**2
        + 2.*q_out*q_side*R_os**2. + 2.*q_out*q_long*R_ol**2.)))

    return gauss


hf = h5py.File(databaseFile, "r")
event_list = list(hf.keys())
nev = len(event_list)
print(f"total number of events: {nev}")

outdata = {}
for groupName in hf.keys():
    print(f"analysis {groupName}")
    outdata[groupName] = {}

    for qCutMax in qCutMaxList:
        for qCutMin in qCutMinList:
            outdata[groupName][f"R3D_qCut_{qCutMin}_{qCutMax}"] = []
            if invFlag:
                outdata[groupName][f"RInv_qCut_{qCutMin}_{qCutMax}"] = []

    for iK, KTbin in enumerate(KT_values):
        for file_i in HBTFileList:
            fileName = f'{file_i}_{KTbin}'
            event_group = hf.get(groupName)
            temp_data = np.nan_to_num(event_group.get(f"{fileName}.dat"))
            ndim = len(temp_data[:, 0])
            avgCorr = np.zeros([ndim, 2], dtype='float32')

            avgCorr[:, 0] = temp_data[:, 0]/(temp_data[:, 1] + EPS)
            avgCorr[:, 1] = 0.01

            outdata[groupName][fileName] = avgCorr

        # 3D fit
        HBTcorr3D = outdata[groupName][f"{HBTFileList[0]}_{KTbin}"]
        guessVals = [1.0, 5., 5., 5., 0.1, 0.1]
        for qCutMax in qCutMaxList:
            for qCutMin in qCutMinList:
                PosIdx = (q3DArr[:, 3] > qCutMin) & (q3DArr[:, 3] < qCutMax)
                fit_params, cov_mat = curve_fit(gaussian_3d,
                                                q3DArr[PosIdx, :],
                                                HBTcorr3D[PosIdx, 0],
                                                p0=guessVals,
                                                sigma=HBTcorr3D[PosIdx, 1],
                                                absolute_sigma=True,
                                                maxfev=10000)
                fit_errors = np.sqrt(np.diag(cov_mat))
                #fit_corr = gaussian_3d(q3DArr, *fit_params)
                R3Dtemp = [KTMid[iK]]
                for x, y in zip(fit_params, fit_errors):
                    R3Dtemp += [x, y]
                outdata[groupName][f"R3D_qCut_{qCutMin}_{qCutMax}"].append(
                    R3Dtemp)

        if invFlag:
            # perform fits to extract HBT radii
            # 1D fit
            HBTcorr1D = outdata[groupName][f"{HBTFileList[1]}_{KTbin}"]
            guessVals = [1.0, 5.]
            for qCutMax in qCutMaxList:
                for qCutMin in qCutMinList:
                    PosIdx = (qArr > qCutMin) & (qArr < qCutMax)
                    fit_params, cov_mat = curve_fit(gaussian_1d,
                                                    qArr[PosIdx],
                                                    HBTcorr1D[PosIdx, 0],
                                                    p0=guessVals,
                                                    sigma=HBTcorr1D[PosIdx, 1],
                                                    absolute_sigma=True)
                    fit_errors = np.sqrt(np.diag(cov_mat))
                    R1Dtemp = [KTMid[iK]]
                    for x, y in zip(fit_params, fit_errors):
                        R1Dtemp += [x, y]
                    outdata[groupName][f"RInv_qCut_{qCutMin}_{qCutMax}"].append(
                        R1Dtemp)

with open('HBTresults.pickle', 'wb') as pf:
    pickle.dump(outdata, pf)

hf.close()
print("Analysis is done.")
