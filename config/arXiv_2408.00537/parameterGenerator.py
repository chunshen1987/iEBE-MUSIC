#!/usr/bin/env python3
"""
    This script translates the posterior chain files in the parameters
    can be read in by the iEBE-MUSIC package.
"""

import pickle
import sys

parameterName = [
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

with open("posteriorChain.pkl", 'rb') as f:
    data = pickle.load(f)
nParamSets = data['chain'].shape[0]
setId = setId % nParamSets
print(f"Using parameter set: {setId}")

paramSet = data['chain'][setId, :]

paramSet[13] = paramSet[13]/paramSet[12]
paramSet[14] = paramSet[14]/paramSet[12]

with open(paramFile, "w") as f:
    for i in range(len(parameterName)):
        f.write("{}  {}\n".format(parameterName[i], paramSet[i]))
