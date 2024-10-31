#!/usr/bin/env python3
"""
    This script translates the posterior chain files in the parameters
    can be read in by the iEBE-MUSIC package.
"""

import pickle

with open("posteriorChainRaw.pkl", 'rb') as f:
    data = pickle.load(f)

# add parameter names
data['parameterName'] = [
    'BG', 'shadowing_factor', 'ylossParam4At2', 'ylossParam4At4',
    'ylossParam4At6', 'ylossParam4var', 'remnant_energy_loss_fraction',
    'lambdaB', 'string_source_sigma_x', 'string_source_sigma_eta',
    'stringTransverseShiftFrac', 'stringPreEqFlowFactor', 'Shear_to_S_ratio',
    'shear_muBf0p2', 'shear_muBf0p4', 'bulk_viscosity_10_max',
    'bulk_viscosity_10_T_peak', 'bulk_viscosity_10_width_high',
    'bulk_viscosity_10_width_low', 'eps_switch',
]

data['chain'][:, 13] = data['chain'][:, 13]/data['chain'][:, 12]
data['chain'][:, 14] = data['chain'][:, 14]/data['chain'][:, 12]

with open("posteriorChain.pkl", 'wb') as f:
    pickle.dump(data, f)
