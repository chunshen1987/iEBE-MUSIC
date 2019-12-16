#!/usr/bin/env python3

import h5py
import sys
from numpy import *

n_order = 7

def help_message():
    print("{0} database_file event_id".format(sys.argv[0]))
    exit(0)

def calcualte_inte_Qn(pT_low, pT_high, data):
    """
        this function calculates the pT-integrated vn in a 
        given pT range (pT_low, pT_high) for every event in the data
    """
    npT = 50
    pT_inte_array = linspace(pT_low, pT_high, npT)
    dpT = pT_inte_array[1] - pT_inte_array[0]
    dN_event = data[:, 2]
    pT_event = data[:, 0]
    dN_interp = exp(interp(pT_inte_array, pT_event, log(dN_event+1e-30)))
    N = 2.*pi*sum(dN_interp*pT_inte_array)*dpT
    temp_vn_array = [N + 1j*0.0]
    for iorder in range(1, n_order):
        vn_real_event = data[:, 4*iorder]
        vn_imag_event = data[:, 4*iorder+2]
        vn_real_interp = interp(pT_inte_array, pT_event, vn_real_event)
        vn_imag_interp = interp(pT_inte_array, pT_event, vn_imag_event)
        Qn_real_inte = 2.*pi*sum(vn_real_interp*dN_interp*pT_inte_array)*dpT
        Qn_imag_inte = 2.*pi*sum(vn_imag_interp*dN_interp*pT_inte_array)*dpT
        Qn_inte = Qn_real_inte + 1j*Qn_imag_inte
        temp_vn_array.append(Qn_inte)
    return(temp_vn_array)

try:
    database_file = str(sys.argv[1])
    event_id = str(sys.argv[2])
except:
    help_message()

h5_data = h5py.File(database_file, "r")
h5_group = h5_data.get("spvn_results_{}".format(event_id))
print("fetching event {0} from the database {1} ...".format(
    event_id, database_file))

vn_filename = 'particle_9999_vndata_diff_eta_-0.5_0.5.dat'
vn_data = h5_group.get(vn_filename)
vn_data = nan_to_num(vn_data)

# use ALICE cut
Qn_vector = calcualte_inte_Qn(0.2, 3.0, vn_data)

# output Qn vectors
output = []
for iorder in range(n_order):
    vn_real = real(Qn_vector[iorder]/Qn_vector[0])
    vn_imag = imag(Qn_vector[iorder]/Qn_vector[0])
    temp = [iorder, real(Qn_vector[iorder]), imag(Qn_vector[iorder]),
            vn_real, vn_imag, sqrt(vn_real**2. + vn_imag**2.),
            arctan2(vn_imag, vn_real)/(float(iorder) + 1e-15)]
    output.append(temp)
savetxt("Qn_vectors_{}.dat".format(event_id), output,
        fmt="%.4e", delimiter="  ",
        header="n  Qn_real  Qn_imag  vn_real  vn_imag  vn_mag  psi_n")
