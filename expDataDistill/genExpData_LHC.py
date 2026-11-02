#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 16 21:41:03 2023

@author: chunshen
"""

#%% setup
import pickle
from os import path

home = path.join(path.expanduser("~"))
desktop = path.join(home, "Desktop")

import numpy as np
import matplotlib.pyplot as plt


def computeCovarianceMatrix(StatErr, SysErr):
    """
        Here we assume the systematic error are fully correlated
    """
    CovMat = np.outer(SysErr, SysErr)
    di = np.diag_indices_from(CovMat)
    CovMat[di] += StatErr**2.
    return CovMat


def packMultipleCov(covList):
    ndim = 0
    for cov_i in covList:
        ndim += cov_i.shape[0]

    outputCov = np.zeros([ndim, ndim])
    offset = 0
    for i, cov_i in enumerate(covList):
        covDim = cov_i.shape[0]
        outputCov[offset:offset + covDim, offset:offset + covDim] = cov_i
        offset = offset + covDim
    return outputCov


def packMultipleObs(obsList, statErrList, sysErrList, covList):
    ndim = 0
    for obs_i in obsList:
        ndim += obs_i.shape[0]

    outputDataArr = np.zeros([ndim, 3])
    outputCov = np.zeros([ndim, ndim])

    offset = 0
    for i, obs_i in enumerate(obsList):
        obsDim = obs_i.shape[0]
        outputDataArr[offset:offset + obsDim, 0] = obs_i
        outputDataArr[offset:offset + obsDim, 1] = statErrList[i]
        outputDataArr[offset:offset + obsDim, 2] = sysErrList[i]
        outputCov[offset:offset + obsDim, offset:offset + obsDim] = covList[i]
        offset = offset + obsDim
    #output = np.concatenate((outputDataArr, outputCov), axis=1)
    outputDict = {"expData": {}}
    outputDict["expData"]["obs"] = outputDataArr.transpose()
    outputDict["expData"]["cov"] = outputCov
    return outputDict


logdNFlag = False  # for future implementation if needed
includev3eta = True
protonFlag = True
includepTVar = False
includepPbdNdeta = True

nameTag = "ALICE"
if logdNFlag:
    nameTag += "_logdN"
if protonFlag:
    nameTag += "_wproton"
if includepTVar:
    nameTag += "_wpTVar"
if includev3eta:
    nameTag += "_wv3eta"
if includepPbdNdeta:
    nameTag += "_wpPbNdeta"

#%% ALICE dNch/deta and pid dN/dy
dNcenCut = 9
vncenCut = 9
v3cenCut = 8
statFrac = 1.0

dNDataArr = []
dNDataStatErr = []
dNDataSysErr = []
dNDataCovList = []

exp_path = path.abspath('../analysisKit/ALICE_run2_dNdy_and_pTSpectra')

# dNch/deta
dNdata = np.loadtxt(path.join(exp_path, "dNch_deta_ALICE.dat"))
dNDataArr.append(dNdata[:dNcenCut, 1])
dNDataStatErr.append(np.sqrt(statFrac)*dNdata[:dNcenCut, 2])
dNDataSysErr.append(np.sqrt(1. - statFrac)*dNdata[:dNcenCut, 2])
dNDataCovList.append(
    computeCovarianceMatrix(
        np.sqrt(statFrac)*dNdata[:dNcenCut, 2],
        np.sqrt(1. - statFrac)*dNdata[:dNcenCut, 2]))

# dN/dy (pi^+ + pi^-)/2
dNdata = np.loadtxt(path.join(exp_path, "PbPb5020_ALICE_pion_dNdy.dat"))
dNDataArr.append(dNdata[:dNcenCut, 1]/2.)
dNDataStatErr.append(dNdata[:dNcenCut, 2]/2.)
dNDataSysErr.append(dNdata[:dNcenCut, 3]/2.)
dNDataCovList.append(
    computeCovarianceMatrix(dNdata[:dNcenCut, 2]/2., dNdata[:dNcenCut, 3]/2.))

# dN/dy (K^+ + K^-)/2
dNdata = np.loadtxt(path.join(exp_path, "PbPb5020_ALICE_kaon_dNdy.dat"))
dNDataArr.append(dNdata[:dNcenCut, 1]/2.)
dNDataStatErr.append(dNdata[:dNcenCut, 2]/2.)
dNDataSysErr.append(dNdata[:dNcenCut, 3]/2.)
dNDataCovList.append(
    computeCovarianceMatrix(dNdata[:dNcenCut, 2]/2., dNdata[:dNcenCut, 3]/2.))

if protonFlag:
    # dN/dy (p + pbar)/2
    dNdata = np.loadtxt(path.join(exp_path, "PbPb5020_ALICE_proton_dNdy.dat"))
    dNDataArr.append(dNdata[:dNcenCut, 1]/2.)
    dNDataStatErr.append(dNdata[:dNcenCut, 2]/2.)
    dNDataSysErr.append(dNdata[:dNcenCut, 3]/2.)
    dNDataCovList.append(
        computeCovarianceMatrix(dNdata[:dNcenCut, 2]/2., dNdata[:dNcenCut, 3]/2.))

    # dN/dy pbar/p ratios
    pbarpRatio = []
    pbarpRatioStatErr = []
    pbarpRatioSysErr = []
    for itable in range(1, 6):
        with open(path.join(exp_path, f"HEPData-ins2725426-v1-Ratios_{itable}.csv"), "r") as f:
            for line in f.readlines():
                if "mathrm{p}" in line:
                    pbarpRatio.append(float(line.split(",")[1]))
                    pbarpRatioStatErr.append(float(line.split(",")[2]))
                    pbarpRatioSysErr.append(float(line.split(",")[4]))
    dNDataArr.append(np.array(pbarpRatio))
    dNDataStatErr.append(np.array(pbarpRatioStatErr))
    dNDataSysErr.append(np.array(pbarpRatioSysErr))
    dNDataCovList.append(
        computeCovarianceMatrix(np.array(pbarpRatioStatErr),
                                np.array(pbarpRatioSysErr)))

dNDataArr = np.concatenate([np.array(arr_i) for arr_i in dNDataArr])
dNDataStatErr = np.concatenate([np.array(arr_i) for arr_i in dNDataStatErr])
dNDataSysErr = np.concatenate([np.array(arr_i) for arr_i in dNDataSysErr])
dNDataCov = packMultipleCov(dNDataCovList)

# mean pT

pTDataArr = []
pTDataStatErr = []
pTDataSysErr = []
pTDataCovList = []

# charged hadron <pT>
pTdata = np.loadtxt(path.join(exp_path, "PbPb_ch_meanpT.dat"))
pTDataArr.append(pTdata[:dNcenCut, 3])
pTDataStatErr.append(np.sqrt(statFrac)*pTdata[:dNcenCut, 4])
pTDataSysErr.append(np.sqrt(1. - statFrac)*pTdata[:dNcenCut, 4])
pTDataCovList.append(
    computeCovarianceMatrix(
        np.sqrt(statFrac)*pTdata[:dNcenCut, 4],
        np.sqrt(1. - statFrac)*pTdata[:dNcenCut, 4]))

# pion <pT>
pTdata = np.loadtxt(path.join(exp_path, "PbPb_5020_meanpT_pion.txt"))
pTDataArr.append(pTdata[:dNcenCut, 1])
pTDataStatErr.append(pTdata[:dNcenCut, 2])
pTDataSysErr.append(pTdata[:dNcenCut, 3])
pTDataCovList.append(
    computeCovarianceMatrix(pTdata[:dNcenCut, 2], pTdata[:dNcenCut, 3]))

# kaon <pT>
pTdata = np.loadtxt(path.join(exp_path, "PbPb_5020_meanpT_kaon.txt"))
pTDataArr.append(pTdata[:dNcenCut, 1])
pTDataStatErr.append(pTdata[:dNcenCut, 2])
pTDataSysErr.append(pTdata[:dNcenCut, 3])
pTDataCovList.append(
    computeCovarianceMatrix(pTdata[:dNcenCut, 2], pTdata[:dNcenCut, 3]))

if protonFlag:
    # proton <pT>
    pTdata = np.loadtxt(path.join(exp_path, "PbPb_5020_meanpT_proton.txt"))
    pTDataArr.append(pTdata[:dNcenCut, 1])
    pTDataStatErr.append(pTdata[:dNcenCut, 2])
    pTDataSysErr.append(pTdata[:dNcenCut, 3])
    pTDataCovList.append(
        computeCovarianceMatrix(pTdata[:dNcenCut, 2], pTdata[:dNcenCut, 3]))

pTDataArr = np.array(pTDataArr).reshape(-1)
pTDataStatErr = np.array(pTDataStatErr).reshape(-1)
pTDataSysErr = np.array(pTDataSysErr).reshape(-1)
pTDataCov = packMultipleCov(pTDataCovList)

# vn{2} at mid-rapidity
exp_path = path.abspath('../analysisKit/ALICE_run2_ch_vn')

vnDataArr = []
vnDataStatErr = []
vnDataSysErr = []
vnDataCovList = []

v2data = np.loadtxt(path.join(exp_path, "HEPData-ins1419244-v2-Table_1.csv"),
                    delimiter=',',
                    skiprows=14,
                    max_rows=9)
vnDataArr.append(v2data[:vncenCut, 3])
vnDataStatErr.append(v2data[:vncenCut, 4])
vnDataSysErr.append(v2data[:vncenCut, 6])
vnDataCovList.append(
    computeCovarianceMatrix(v2data[:vncenCut, 4], v2data[:vncenCut, 6]))

v3data = np.loadtxt(path.join(exp_path, "HEPData-ins1419244-v2-Table_2.csv"),
                    delimiter=',',
                    skiprows=14,
                    max_rows=8)
vnDataArr.append(v3data[:v3cenCut, 3])
vnDataStatErr.append(v3data[:v3cenCut, 4])
vnDataSysErr.append(v3data[:v3cenCut, 6])
vnDataCovList.append(
    computeCovarianceMatrix(v3data[:v3cenCut, 4], v3data[:v3cenCut, 6]))

vnDataArr = np.concatenate([np.array(arr_i) for arr_i in vnDataArr])
vnDataStatErr = np.concatenate([np.array(arr_i) for arr_i in vnDataStatErr])
vnDataSysErr = np.concatenate([np.array(arr_i) for arr_i in vnDataSysErr])
vnDataCov = packMultipleCov(vnDataCovList)

# pT variance at mid-rapidity
exp_path = path.abspath('../analysisKit/ALICE_run2_pTvar')

pTVarDataArr = []
pTVarDataStatErr = []
pTVarDataSysErr = []

pTVarData = np.loadtxt(path.join(exp_path, "HEPData-ins2848476-v1-Table_2.csv"),
                       delimiter=',',
                       skiprows=13)
pTVarDataArr.append(pTVarData[-12:, 3])
pTVarDataStatErr.append(pTVarData[-12:, 4])
pTVarDataSysErr.append(pTVarData[-12:, 6])

pTVarDataArr = np.array(pTVarDataArr).reshape(-1)
pTVarDataStatErr = np.array(pTVarDataStatErr).reshape(-1)
pTVarDataSysErr = np.array(pTVarDataSysErr).reshape(-1)
pTVarDataCov = computeCovarianceMatrix(pTVarDataStatErr, pTVarDataSysErr)

# dNch/deta vs eta
exp_path = path.abspath('../analysisKit/ALICE_run2_dNchdeta_rapdep')

dNchdetaDataArr = []
dNchdetaDataStatErr = []
dNchdetaDataSysErr = []
dNchdetaDataCovList = []

# 0-5%
dNchdetaData = np.loadtxt(path.join(exp_path,
                                    "HEPData-ins1507090-v1-DN_DETARAP.csv"),
                          delimiter=',',
                          skiprows=13,
                          max_rows=34)
dNchdetaDataArr.append(dNchdetaData[:, 3])
dNchdetaDataStatErr.append(dNchdetaData[:, 4])
dNchdetaDataSysErr.append(dNchdetaData[:, 6])
dNchdetaDataCovList.append(
    computeCovarianceMatrix(dNchdetaData[:, 4], dNchdetaData[:, 6]))

# 5-10%
dNchdetaData = np.loadtxt(path.join(exp_path,
                                    "HEPData-ins1507090-v1-DN_DETARAP.csv"),
                          delimiter=',',
                          skiprows=53,
                          max_rows=34)
dNchdetaDataArr.append(dNchdetaData[:, 3])
dNchdetaDataStatErr.append(dNchdetaData[:, 4])
dNchdetaDataSysErr.append(dNchdetaData[:, 6])
dNchdetaDataCovList.append(
    computeCovarianceMatrix(dNchdetaData[:, 4], dNchdetaData[:, 6]))

# 10-20%
dNchdetaData = np.loadtxt(path.join(exp_path,
                                    "HEPData-ins1507090-v1-DN_DETARAP.csv"),
                          delimiter=',',
                          skiprows=93,
                          max_rows=34)
dNchdetaDataArr.append(dNchdetaData[:, 3])
dNchdetaDataStatErr.append(dNchdetaData[:, 4])
dNchdetaDataSysErr.append(dNchdetaData[:, 6])
dNchdetaDataCovList.append(
    computeCovarianceMatrix(dNchdetaData[:, 4], dNchdetaData[:, 6]))

# 20-30%
dNchdetaData = np.loadtxt(path.join(exp_path,
                                    "HEPData-ins1507090-v1-DN_DETARAP.csv"),
                          delimiter=',',
                          skiprows=133,
                          max_rows=34)
dNchdetaDataArr.append(dNchdetaData[:, 3])
dNchdetaDataStatErr.append(dNchdetaData[:, 4])
dNchdetaDataSysErr.append(dNchdetaData[:, 6])
dNchdetaDataCovList.append(
    computeCovarianceMatrix(dNchdetaData[:, 4], dNchdetaData[:, 6]))

# 30-40%
dNchdetaData = np.loadtxt(path.join(exp_path,
                                    "HEPData-ins1507090-v1-DN_DETARAP.csv"),
                          delimiter=',',
                          skiprows=173,
                          max_rows=34)
dNchdetaDataArr.append(dNchdetaData[:, 3])
dNchdetaDataStatErr.append(dNchdetaData[:, 4])
dNchdetaDataSysErr.append(dNchdetaData[:, 6])
dNchdetaDataCovList.append(
    computeCovarianceMatrix(dNchdetaData[:, 4], dNchdetaData[:, 6]))

# 40-50%
dNchdetaData = np.loadtxt(path.join(exp_path,
                                    "HEPData-ins1507090-v1-DN_DETARAP.csv"),
                          delimiter=',',
                          skiprows=213,
                          max_rows=34)
dNchdetaDataArr.append(dNchdetaData[:, 3])
dNchdetaDataStatErr.append(dNchdetaData[:, 4])
dNchdetaDataSysErr.append(dNchdetaData[:, 6])
dNchdetaDataCovList.append(
    computeCovarianceMatrix(dNchdetaData[:, 4], dNchdetaData[:, 6]))

dNchdetaDataArr = np.array(dNchdetaDataArr).reshape(-1)
dNchdetaDataStatErr = np.array(dNchdetaDataStatErr).reshape(-1)
dNchdetaDataSysErr = np.array(dNchdetaDataSysErr).reshape(-1)
dNchdetaDataCov = packMultipleCov(dNchdetaDataCovList)

# vn(eta) vs eta
exp_path = path.abspath('../analysisKit/ALICE_run2_ch_vneta')

vnchdetaDataArr = []
vnchdetaDataStatErr = []
vnchdetaDataSysErr = []
vnchdetaDataCovList = []

# v2(eta) 0-5%
vnchdetaData = np.loadtxt(path.join(exp_path,
                                    "HEPData-ins2679248-v1-Table_1.csv"),
                          delimiter=',',
                          skiprows=14,
                          max_rows=30)
vnchdetaDataArr.append(vnchdetaData[:, 3])
vnchdetaDataStatErr.append(vnchdetaData[:, 4])
vnchdetaDataSysErr.append(vnchdetaData[:, 6])
vnchdetaDataCovList.append(
    computeCovarianceMatrix(vnchdetaData[:, 4], vnchdetaData[:, 6]))

# v2(eta) 5-10%
vnchdetaData = np.loadtxt(path.join(exp_path,
                                    "HEPData-ins2679248-v1-Table_1.csv"),
                          delimiter=',',
                          skiprows=51,
                          max_rows=30)
vnchdetaDataArr.append(vnchdetaData[:, 3])
vnchdetaDataStatErr.append(vnchdetaData[:, 4])
vnchdetaDataSysErr.append(vnchdetaData[:, 6])
vnchdetaDataCovList.append(
    computeCovarianceMatrix(vnchdetaData[:, 4], vnchdetaData[:, 6]))

# v2(eta) 10-20%
vnchdetaData = np.loadtxt(path.join(exp_path,
                                    "HEPData-ins2679248-v1-Table_1.csv"),
                          delimiter=',',
                          skiprows=88,
                          max_rows=30)
vnchdetaDataArr.append(vnchdetaData[:, 3])
vnchdetaDataStatErr.append(vnchdetaData[:, 4])
vnchdetaDataSysErr.append(vnchdetaData[:, 6])
vnchdetaDataCovList.append(
    computeCovarianceMatrix(vnchdetaData[:, 4], vnchdetaData[:, 6]))

# v2(eta) 20-30%
vnchdetaData = np.loadtxt(path.join(exp_path,
                                    "HEPData-ins2679248-v1-Table_1.csv"),
                          delimiter=',',
                          skiprows=125,
                          max_rows=30)
vnchdetaDataArr.append(vnchdetaData[:, 3])
vnchdetaDataStatErr.append(vnchdetaData[:, 4])
vnchdetaDataSysErr.append(vnchdetaData[:, 6])
vnchdetaDataCovList.append(
    computeCovarianceMatrix(vnchdetaData[:, 4], vnchdetaData[:, 6]))

# v2(eta) 30-40%
vnchdetaData = np.loadtxt(path.join(exp_path,
                                    "HEPData-ins2679248-v1-Table_1.csv"),
                          delimiter=',',
                          skiprows=162,
                          max_rows=30)
vnchdetaDataArr.append(vnchdetaData[:, 3])
vnchdetaDataStatErr.append(vnchdetaData[:, 4])
vnchdetaDataSysErr.append(vnchdetaData[:, 6])
vnchdetaDataCovList.append(
    computeCovarianceMatrix(vnchdetaData[:, 4], vnchdetaData[:, 6]))

# v2(eta) 40-50%
vnchdetaData = np.loadtxt(path.join(exp_path,
                                    "HEPData-ins2679248-v1-Table_1.csv"),
                          delimiter=',',
                          skiprows=199,
                          max_rows=30)
vnchdetaDataArr.append(vnchdetaData[:, 3])
vnchdetaDataStatErr.append(vnchdetaData[:, 4])
vnchdetaDataSysErr.append(vnchdetaData[:, 6])
vnchdetaDataCovList.append(
    computeCovarianceMatrix(vnchdetaData[:, 4], vnchdetaData[:, 6]))

if includev3eta:
    # v3(eta) 0-5%
    vnchdetaData = np.loadtxt(path.join(exp_path,
                                        "HEPData-ins2679248-v1-Table_2.csv"),
                              delimiter=',',
                              skiprows=14,
                              max_rows=30)
    vnchdetaDataArr.append(vnchdetaData[:, 3])
    vnchdetaDataStatErr.append(vnchdetaData[:, 4])
    vnchdetaDataSysErr.append(vnchdetaData[:, 6])
    vnchdetaDataCovList.append(
        computeCovarianceMatrix(vnchdetaData[:, 4], vnchdetaData[:, 6]))

    # v3(eta) 5-10%
    vnchdetaData = np.loadtxt(path.join(exp_path,
                                        "HEPData-ins2679248-v1-Table_2.csv"),
                              delimiter=',',
                              skiprows=51,
                              max_rows=30)
    vnchdetaDataArr.append(vnchdetaData[:, 3])
    vnchdetaDataStatErr.append(vnchdetaData[:, 4])
    vnchdetaDataSysErr.append(vnchdetaData[:, 6])
    vnchdetaDataCovList.append(
        computeCovarianceMatrix(vnchdetaData[:, 4], vnchdetaData[:, 6]))

    # v3(eta) 10-20%
    vnchdetaData = np.loadtxt(path.join(exp_path,
                                        "HEPData-ins2679248-v1-Table_2.csv"),
                              delimiter=',',
                              skiprows=88,
                              max_rows=30)
    vnchdetaDataArr.append(vnchdetaData[:, 3])
    vnchdetaDataStatErr.append(vnchdetaData[:, 4])
    vnchdetaDataSysErr.append(vnchdetaData[:, 6])
    vnchdetaDataCovList.append(
        computeCovarianceMatrix(vnchdetaData[:, 4], vnchdetaData[:, 6]))

    # v3(eta) 20-30%
    vnchdetaData = np.loadtxt(path.join(exp_path,
                                        "HEPData-ins2679248-v1-Table_2.csv"),
                              delimiter=',',
                              skiprows=125,
                              max_rows=30)
    vnchdetaDataArr.append(vnchdetaData[:, 3])
    vnchdetaDataStatErr.append(vnchdetaData[:, 4])
    vnchdetaDataSysErr.append(vnchdetaData[:, 6])
    vnchdetaDataCovList.append(
        computeCovarianceMatrix(vnchdetaData[:, 4], vnchdetaData[:, 6]))

    # v3(eta) 30-40%
    vnchdetaData = np.loadtxt(path.join(exp_path,
                                        "HEPData-ins2679248-v1-Table_2.csv"),
                              delimiter=',',
                              skiprows=162,
                              max_rows=30)
    vnchdetaDataArr.append(vnchdetaData[:, 3])
    vnchdetaDataStatErr.append(vnchdetaData[:, 4])
    vnchdetaDataSysErr.append(vnchdetaData[:, 6])
    vnchdetaDataCovList.append(
        computeCovarianceMatrix(vnchdetaData[:, 4], vnchdetaData[:, 6]))

    # v3(eta) 40-50%
    vnchdetaData = np.loadtxt(path.join(exp_path,
                                        "HEPData-ins2679248-v1-Table_2.csv"),
                              delimiter=',',
                              skiprows=199,
                              max_rows=30)
    vnchdetaDataArr.append(vnchdetaData[:, 3])
    vnchdetaDataStatErr.append(vnchdetaData[:, 4])
    vnchdetaDataSysErr.append(vnchdetaData[:, 6])
    vnchdetaDataCovList.append(
        computeCovarianceMatrix(vnchdetaData[:, 4], vnchdetaData[:, 6]))

vnchdetaDataArr = np.array(vnchdetaDataArr).reshape(-1)
vnchdetaDataStatErr = np.array(vnchdetaDataStatErr).reshape(-1)
vnchdetaDataSysErr = np.array(vnchdetaDataSysErr).reshape(-1)
vnchdetaDataCov = packMultipleCov(vnchdetaDataCovList)

if includepPbdNdeta:
    # pPb dNch/deta vs eta
    exp_path = path.abspath('../analysisKit/ALICE_pPb5TeV_dNchdeta_rapdep')

    pPbdNchdetaDataArr = []
    pPbdNchdetaDataStatErr = []
    pPbdNchdetaDataSysErr = []
    pPbdNchdetaDataCovList = []

    # 0-5%
    pPbdNchdetaData = np.loadtxt(path.join(exp_path,
                                           "HEPData-ins2070408-v1-Table_2.csv"),
                              delimiter=',',
                              skiprows=13,
                              max_rows=64)
    pPbdNchdetaDataArr.append(pPbdNchdetaData[:, 3])
    pPbdNchdetaDataStatErr.append(pPbdNchdetaData[:, 4])
    pPbdNchdetaDataSysErr.append(pPbdNchdetaData[:, 6])
    pPbdNchdetaDataCovList.append(
        computeCovarianceMatrix(pPbdNchdetaData[:, 4], pPbdNchdetaData[:, 6]))
    # 5-10%
    pPbdNchdetaData = np.loadtxt(path.join(exp_path,
                                           "HEPData-ins2070408-v1-Table_2.csv"),
                              delimiter=',',
                              skiprows=83,
                              max_rows=64)
    pPbdNchdetaDataArr.append(pPbdNchdetaData[:, 3])
    pPbdNchdetaDataStatErr.append(pPbdNchdetaData[:, 4])
    pPbdNchdetaDataSysErr.append(pPbdNchdetaData[:, 6])
    pPbdNchdetaDataCovList.append(
        computeCovarianceMatrix(pPbdNchdetaData[:, 4], pPbdNchdetaData[:, 6]))
    # 10-20%
    pPbdNchdetaData = np.loadtxt(path.join(exp_path,
                                           "HEPData-ins2070408-v1-Table_2.csv"),
                              delimiter=',',
                              skiprows=153,
                              max_rows=64)
    pPbdNchdetaDataArr.append(pPbdNchdetaData[:, 3])
    pPbdNchdetaDataStatErr.append(pPbdNchdetaData[:, 4])
    pPbdNchdetaDataSysErr.append(pPbdNchdetaData[:, 6])
    pPbdNchdetaDataCovList.append(
        computeCovarianceMatrix(pPbdNchdetaData[:, 4], pPbdNchdetaData[:, 6]))
    # 20-40%
    pPbdNchdetaData = np.loadtxt(path.join(exp_path,
                                           "HEPData-ins2070408-v1-Table_2.csv"),
                              delimiter=',',
                              skiprows=223,
                              max_rows=64)
    pPbdNchdetaDataArr.append(pPbdNchdetaData[:, 3])
    pPbdNchdetaDataStatErr.append(pPbdNchdetaData[:, 4])
    pPbdNchdetaDataSysErr.append(pPbdNchdetaData[:, 6])
    pPbdNchdetaDataCovList.append(
        computeCovarianceMatrix(pPbdNchdetaData[:, 4], pPbdNchdetaData[:, 6]))
    # 40-60%
    pPbdNchdetaData = np.loadtxt(path.join(exp_path,
                                           "HEPData-ins2070408-v1-Table_2.csv"),
                              delimiter=',',
                              skiprows=293,
                              max_rows=64)
    pPbdNchdetaDataArr.append(pPbdNchdetaData[:, 3])
    pPbdNchdetaDataStatErr.append(pPbdNchdetaData[:, 4])
    pPbdNchdetaDataSysErr.append(pPbdNchdetaData[:, 6])
    pPbdNchdetaDataCovList.append(
        computeCovarianceMatrix(pPbdNchdetaData[:, 4], pPbdNchdetaData[:, 6]))
    # 60-80%
    pPbdNchdetaData = np.loadtxt(path.join(exp_path,
                                           "HEPData-ins2070408-v1-Table_2.csv"),
                              delimiter=',',
                              skiprows=363,
                              max_rows=64)
    pPbdNchdetaDataArr.append(pPbdNchdetaData[:, 3])
    pPbdNchdetaDataStatErr.append(pPbdNchdetaData[:, 4])
    pPbdNchdetaDataSysErr.append(pPbdNchdetaData[:, 6])
    pPbdNchdetaDataCovList.append(
        computeCovarianceMatrix(pPbdNchdetaData[:, 4], pPbdNchdetaData[:, 6]))

    pPbdNchdetaDataArr = np.array(pPbdNchdetaDataArr).reshape(-1)
    pPbdNchdetaDataStatErr = np.array(pPbdNchdetaDataStatErr).reshape(-1)
    pPbdNchdetaDataSysErr = np.array(pPbdNchdetaDataSysErr).reshape(-1)
    pPbdNchdetaDataCov = packMultipleCov(pPbdNchdetaDataCovList)

#################
### Save data ###
#################

outputDataList = [
    dNDataArr,
    pTDataArr,
    vnDataArr,
    dNchdetaDataArr,
    vnchdetaDataArr,
]
outputDataStatErrList = [
    dNDataStatErr,
    pTDataStatErr,
    vnDataStatErr,
    dNchdetaDataStatErr,
    vnchdetaDataStatErr,
]
outputDataSysErrList = [
    dNDataSysErr,
    pTDataSysErr,
    vnDataSysErr,
    dNchdetaDataSysErr,
    vnchdetaDataSysErr,
]
outputDataCovList = [
    dNDataCov,
    pTDataCov,
    vnDataCov,
    dNchdetaDataCov,
    vnchdetaDataCov,
]

if includepTVar:
    outputDataList.append(pTVarDataArr)
    outputDataStatErrList.append(pTVarDataStatErr)
    outputDataSysErrList.append(pTVarDataSysErr)
    outputDataCovList.append(pTVarDataCov)

if includepPbdNdeta:
    outputDataList.append(pPbdNchdetaDataArr)
    outputDataStatErrList.append(pPbdNchdetaDataStatErr)
    outputDataSysErrList.append(pPbdNchdetaDataSysErr)
    outputDataCovList.append(pPbdNchdetaDataCov)

outputDict = packMultipleObs(outputDataList, outputDataStatErrList,
                             outputDataSysErrList, outputDataCovList)
with open(f"exp_data_{nameTag}.pkl", "wb") as pf:
    pickle.dump(outputDict, pf)
