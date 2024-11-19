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
    'ylossParam4At6', 'ylossParam4var', 'remnant_energy_loss_fraction_200',
    'remnant_energy_loss_fraction_19p6', 'remnant_energy_loss_fraction_7p7',
    'lambdaB_200', 'lambdaB_19p6', 'lambdaB_7p7',
    'string_source_sigma_x_200', 'string_source_sigma_x_19p6',
    'string_source_sigma_x_7p7', 'string_source_sigma_eta',
    'stringTransverseShiftFrac', 'stringPreEqFlowFactor', 'Shear_to_S_ratio',
    'shear_muB_0p2', 'shear_muB_0p4', 'bulk_viscosity_10_max',
    'bulk_viscosity_10_T_peak', 'bulk_viscosity_10_width_high',
    'bulk_viscosity_10_width_low', 'eps_switch_200', 'eps_switch_19p6',
    'eps_switch_7p7',
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
paramFile = str(sys.argv[2])
ecm = float(sys.argv[3])

with open("posteriorChain.pkl", 'rb') as f:
    data = pickle.load(f)
nParamSets = data['chain'].shape[0]
setId = setId % nParamSets
print(f"Using parameter set: {setId}")

paramSet = data['chain'][setId, :]
paramDict = {}
for i, param_i in enumerate(parameterName):
    paramDict[param_i] = paramSet[i]

paramDict['shear_muBf0p2'] = (paramDict['shear_muB_0p2']
                              /paramDict['Shear_to_S_ratio'])
paramDict['shear_muBf0p4'] = (paramDict['shear_muB_0p4']
                              /paramDict['Shear_to_S_ratio'])
if ecm < 7.7:
    paramDict['remnant_energy_loss_fraction'] = (
        paramDict['remnant_energy_loss_fraction_7p7'])
    paramDict['lambdaB'] = paramDict['lambdaB_7p7']
    paramDict['string_source_sigma_x'] = paramDict['string_source_sigma_x_7p7']
    paramDict['eps_switch'] = paramDict['eps_switch_7p7']
elif ecm < 19.6:
    frac = (np.log(ecm) - np.log(7.7)) / (np.log(19.6) - np.log(7.7))
    paramDict['remnant_energy_loss_fraction'] = (
        (1 - frac) * paramDict['remnant_energy_loss_fraction_7p7']
        + frac * paramDict['remnant_energy_loss_fraction_19p6']
    )
    paramDict['lambdaB'] = (
        (1 - frac) * paramDict['lambdaB_7p7'] + frac * paramDict['lambdaB_19p6']
    )
    paramDict['string_source_sigma_x'] = (
        (1 - frac) * paramDict['string_source_sigma_x_7p7']
        + frac * paramDict['string_source_sigma_x_19p6']
    )
    paramDict['eps_switch'] = (
        (1 - frac) * paramDict['eps_switch_7p7']
        + frac * paramDict['eps_switch_19p6']
    )
elif ecm < 200.:
    frac = (np.log(ecm) - np.log(19.6)) / (np.log(200.) - np.log(19.6))
    paramDict['remnant_energy_loss_fraction'] = (
        (1 - frac) * paramDict['remnant_energy_loss_fraction_19p6']
        + frac * paramDict['remnant_energy_loss_fraction_200']
    )
    paramDict['lambdaB'] = (
        (1 - frac) * paramDict['lambdaB_19p6'] + frac * paramDict['lambdaB_200']
    )
    paramDict['string_source_sigma_x'] = (
        (1 - frac) * paramDict['string_source_sigma_x_19p6']
        + frac * paramDict['string_source_sigma_x_200']
    )
    paramDict['eps_switch'] = (
        (1 - frac) * paramDict['eps_switch_19p6']
        + frac * paramDict['eps_switch_200']
    )
else:
    paramDict['remnant_energy_loss_fraction'] = (
        paramDict['remnant_energy_loss_fraction_200']
    )
    paramDict['lambdaB'] = paramDict['lambdaB_200']
    paramDict['string_source_sigma_x'] = paramDict['string_source_sigma_x_200']
    paramDict['eps_switch'] = paramDict['eps_switch_200']

with open(paramFile, "w") as f:
    for param_i in outputParameterName:
        f.write("{}  {}\n".format(param_i, paramDict[param_i]))
