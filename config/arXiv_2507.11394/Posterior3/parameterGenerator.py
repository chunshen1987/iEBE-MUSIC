#!/usr/bin/env python3
"""
    This script translates the posterior chain files in the parameters
    can be read in by the iEBE-MUSIC package.
"""

import pickle
import sys
import numpy as np

parameterName = [
    'BG', 'shadowing_factor', 'ylossParam4At2', 'ylossParam4At4',
    'ylossParam4At6', 'ylossParam4var', 'remnant_energy_loss_fraction',
    'lambdaB', 'string_source_sigma_x_200', 'string_source_sigma_x_19p6',
    'string_source_sigma_x_7p7', 'string_source_sigma_eta_200',
    'string_source_sigma_eta_19p6', 'string_source_sigma_eta_7p7',
    'stringTransverseShiftFrac', 'stringPreEqFlowFactor',
    'Shear_to_S_ratio', 'shear_muB_0p2', 'shear_muB_0p4',
    'bulk_viscosity_10_max', 'bulk_viscosity_10_T_peak',
    'bulk_viscosity_10_width_high', 'bulk_viscosity_10_width_low',
    'eps_switch'
]

outputParameterName = [
    'BG', 'shadowing_factor', 'ylossParam4At2', 'ylossParam4At4',
    'ylossParam4At6', 'ylossParam4var', 'remnant_energy_loss_fraction',
    'lambdaB', 'string_source_sigma_x', 'string_source_sigma_eta',
    'stringTransverseShiftFrac', 'stringPreEqFlowFactor', 'Shear_to_S_ratio',
    'shear_muBf0p2', 'shear_muBf0p4', 'bulk_viscosity_10_max',
    'bulk_viscosity_10_T_peak', 'bulk_viscosity_10_width_high',
    'bulk_viscosity_10_width_low', 'eps_switch',
]

setId = int(sys.argv[1])
setFlag = int(sys.argv[2])
paramFile = str(sys.argv[3])
ecm = float(sys.argv[4])

with open("posteriorChain.pkl", 'rb') as f:
    data = pickle.load(f)

setName = 'chain'
if setFlag == 1:
    setName = 'paramClusters'

nParamSets = data[setName].shape[0]
setId = setId % nParamSets
paramSet = data[setName][setId, :]
print(f"Using parameter set: {setId} from {setName}")
paramDict = {}
for i, param_i in enumerate(parameterName):
    paramDict[param_i] = paramSet[i]

paramDict['shear_muBf0p2'] = (paramDict['shear_muB_0p2']
                              /paramDict['Shear_to_S_ratio'])
paramDict['shear_muBf0p4'] = (paramDict['shear_muB_0p4']
                              /paramDict['Shear_to_S_ratio'])
if ecm < 7.7:
    paramDict['string_source_sigma_x'] = paramDict['string_source_sigma_x_7p7']
    paramDict['string_source_sigma_eta'] = (
                            paramDict['string_source_sigma_eta_7p7'])
elif ecm < 19.6:
    frac = (np.log(ecm) - np.log(7.7)) / (np.log(19.6) - np.log(7.7))
    paramDict['string_source_sigma_x'] = (
        (1 - frac) * paramDict['string_source_sigma_x_7p7']
        + frac * paramDict['string_source_sigma_x_19p6']
    )
    paramDict['string_source_sigma_eta'] = (
        (1 - frac) * paramDict['string_source_sigma_eta_7p7']
        + frac * paramDict['string_source_sigma_eta_19p6']
    )
elif ecm < 200.:
    frac = (np.log(ecm) - np.log(19.6)) / (np.log(200.) - np.log(19.6))
    paramDict['string_source_sigma_x'] = (
        (1 - frac) * paramDict['string_source_sigma_x_19p6']
        + frac * paramDict['string_source_sigma_x_200']
    )
    paramDict['string_source_sigma_eta'] = (
        (1 - frac) * paramDict['string_source_sigma_eta_19p6']
        + frac * paramDict['string_source_sigma_eta_200']
    )
else:
    paramDict['string_source_sigma_x'] = paramDict['string_source_sigma_x_200']
    paramDict['string_source_sigma_eta'] = (
                            paramDict['string_source_sigma_eta_200'])

with open(paramFile, "w") as f:
    for param_i in outputParameterName:
        f.write("{}  {}\n".format(param_i, paramDict[param_i]))
