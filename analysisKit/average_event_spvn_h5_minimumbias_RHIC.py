#! /usr/bin/env python3
"""
     This script performs event averaging for particle 
     spectra and anisotropic flow coefficients calculated 
     from event-by-event simulations

     v_n is analyzed up to n = 6

     Format for particle_XXX_vndata.dat file:
     n_order  real_part  real_part_err  imag_part  imag_part_err

     Format for particle_XXX_vndata_diff.dat file:
     pT(GeV)  pT_err(GeV)  dN/(2pi dy pT dpT)(GeV^-2)  dN/(2pi dy pT dpT)_err(GeV^-2)
     vn_real  vn_real_err  vn_imag  vn_imag_err

     All the errors are only statistic errors
"""

from sys import argv, exit
from os import path, mkdir
from glob import glob
from numpy import *
import h5py
import shutil

# define colors
purple = "\033[95m"
green = "\033[92m"
blue = "\033[94m"
yellow = "\033[93m"
red = "\033[91m"
normal = "\033[0m"

kinematicCutsDict = {
        "PHENIX": {"pTmin": 0.20, "pTmax": 2.0},
        "STAR1" : {"pTmin": 0.15, "pTmax": 2.0},
        "STAR2" : {"pTmin": 0.20, "pTmax": 3.0},
}

Reg_centrality_cut_list = [0., 5., 10., 20., 30., 40., 50.,
                           60., 70., 80., 90., 100.]
PHOBOS_cen_list = [0., 6., 15., 25., 35., 45., 55.]  # PHOBOS AuAu 200
SPS_cen_list    = [5., 12.5, 23.5, 33.5, 43.5]       # SPS PbPb
PHENIX_cen_list = [0., 20., 40., 60., 88.]           # PHENIX dAu
STAR_cen_list   = [0., 10., 40., 80]                 # STAR v1
centralityCutList = Reg_centrality_cut_list
dNcutList = []    # pre-defined Nch cut if simulation is not minimum bias

#centralityCutList = Reg_centrality_cut_list + PHOBOS_cen_list
#dNcutList = [1e5, 544.20, 448.98, 309.78, 204.20, 131.71, 83.02, 45.85,
#             24.89, 11.73, 4.58, 0.16,
#             1e5, 517.67, 362.21, 254.18, 161.67, 108.68, 61.31]

CentralityFlag = 1   # 0: use pre-generated centrality label in the database
                     # 1: sort dNch/deta and cut centrality
RapidityTrigger = 0  # 0: mid-rapidity [-0.5, 0.5]
                     # 1: PHENIX BBC trigger [-3.9, -3.1]
                     # 2: ALICE V0A trigger [-5.1, -2.8]
                     # 3: ATLAS forward trigger [-4.9, -3.1]
FastFlag = True      # True: only analyze a subset of charged hadron obs.
                     # False: full analysis

RapTrigLabel = "CL1"
if RapidityTrigger == 1:
    RapTrigLabel = "BBC"

try:
    data_path = path.abspath(argv[1])
    data_name = data_path.split("/")[-1]
    data_path = "/".join(data_path.split("/")[0:-1])
    resultsFolderName = data_name.split(".h5")[0]
    avg_folder_header = path.join(
        path.abspath(argv[2]), "{}_{}".format(resultsFolderName, RapTrigLabel))
    print("output folder: %s" % avg_folder_header)
    if path.isdir(avg_folder_header):
        print("folder %s already exists!" % avg_folder_header)
        var = input("do you want to delete it? [y/N]")
        if 'y' in var.lower():
            shutil.rmtree(avg_folder_header)
            mkdir(avg_folder_header)
        else:
            print("Continue analysis in {} ...".format(avg_folder_header))
    else:
        mkdir(avg_folder_header)
    paraFiles = glob(path.join(data_path, "parameter*"))
    for iFile in paraFiles:
        shutil.copy(iFile, avg_folder_header)
except IndexError:
    print("Usage: {} working_folder results_folder".format(argv[0]))
    exit(1)

particle_list = ['9999', '211', '321', '2212', '-211', '-321', '-2212', 
                 '3122', '-3122', '3312', '-3312', '3334', '-3334',
                 '333']
particle_name_list = ['charged_hadron', 'pion_p', 'kaon_p', 'proton',
                      'pion_m', 'kaon_m', 'anti_proton',
                      'Lambda', 'anti_Lambda', 'Xi_m', 'anti_Xi_p',
                      'Omega', 'anti_Omega', 'phi']

nonlinear_reponse_correlator_name_list = [
                'v4_L', 'v4(Psi2)', 'rho_422', 'chi_422',
                'v5_L', 'v5(Psi23)', 'rho_523', 'chi_523',
                'v6(Psi2)', 'v6(Psi3)', 'v6(Psi24)',
                'rho_6222', 'rho_633', 'chi_6222', 'chi_633', 'chi_624',
                'v7(Psi23)', 'rho_7223', 'chi_7223']
symmetric_cumulant_name_list = ['SC_32', 'SC_42']

n_order = 10
if FastFlag:
    particle_list = particle_list[0:4]


def check_an_event_is_good(h5_event):
    """This function checks the given event contains all required files"""
    required_files_list = [
        'particle_9999_vndata_eta_-0.5_0.5.dat',
        'particle_9999_vndata_eta_-2.5_-0.5.dat',
        'particle_9999_vndata_eta_0.5_2.5.dat',
        'particle_211_vndata_diff_y_-0.5_0.5.dat',
        'particle_321_vndata_diff_y_-0.5_0.5.dat',
        'particle_2212_vndata_diff_y_-0.5_0.5.dat',
        'particle_3122_vndata_diff_y_-0.5_0.5.dat',
        'particle_-3122_vndata_diff_y_-0.5_0.5.dat',
    ]
    event_file_list = list(h5_event.keys())
    for ifile in required_files_list:
        if ifile not in event_file_list:
            print("event {} is bad, missing {} ...".format(h5_event.name,
                                                           ifile))
            return False
    return True


def calcualte_inte_vn(pT_low, pT_high, data):
    """
        this function calculates the pT-integrated vn in a 
        given pT range (pT_low, pT_high) for every event in the data
    """
    npT = 50
    pT_inte_array = linspace(pT_low, pT_high, npT)
    dpT = pT_inte_array[1] - pT_inte_array[0]
    dN_event = data[:, 1]
    pT_event = data[:, 0]
    dN_interp = exp(interp(pT_inte_array, pT_event, log(dN_event+1e-30)))
    N_event = data[:, -1]
    N_interp = exp(interp(pT_inte_array, pT_event, log(N_event+1e-30)))
    N = sum(N_interp)*dpT/0.1
    temp_vn_array = [N,]
    for iorder in range(1, n_order):
        vn_real_event = data[:, 2*iorder]
        vn_imag_event = data[:, 2*iorder+1]
        vn_real_interp = interp(pT_inte_array, pT_event, vn_real_event)
        vn_imag_interp = interp(pT_inte_array, pT_event, vn_imag_event)
        vn_real_inte = (
            sum(vn_real_interp*dN_interp*pT_inte_array)
            /sum(dN_interp*pT_inte_array))
        vn_imag_inte = (
            sum(vn_imag_interp*dN_interp*pT_inte_array)
            /sum(dN_interp*pT_inte_array))
        vn_inte = vn_real_inte + 1j*vn_imag_inte
        temp_vn_array.append(vn_inte)
    return(temp_vn_array)


def calculate_chi_422(vn_array):
    """chi_422 = Re(v4*conj(v2)**2.)/|v2|^4
       v_422 = Re(v4*conj(v2)**2.)/sqrt(|v2|^4)
       rho_422 = v_422/v4(Psi4)
       v4_L = sqrt(v4(Psi4)^2 - v4(Psi2)^2)
    """
    dN = real(vn_array[:, 0])
    Q2 = dN*vn_array[:, 2]
    Q4 = dN*vn_array[:, 4]
    nev = len(dN)

    N2_weight = dN*(dN - 1.)
    Q4_2 = abs(Q4)**2. - dN

    N4_weight = dN*(dN - 1.)*(dN - 2.)*(dN - 3.)
    Q2_4 = ((abs(Q2)**4.) - 2.*real(Q4*conj(Q2)*conj(Q2))
             - 4.*(dN - 2.)*(abs(Q2)**2.) + abs(Q4)**2.
             + 2*dN*(dN - 3.))

    N3_weight = dN*(dN - 1.)*(dN - 2.)
    chi_422_num = Q4*conj(Q2)*conj(Q2) - 2.*Q2*conj(Q2) - Q4*conj(Q4) + 2.*dN

    chi_422_JK = zeros(nev)
    v422_JK = zeros(nev)
    rho422_JK = zeros(nev)
    v4L_JK = zeros(nev)
    for iev in range(nev):
        array_idx = [True]*nev
        array_idx[iev] = False
        array_idx = array(array_idx)

        num_JK = real(mean(chi_422_num[array_idx]))/mean(N3_weight[array_idx])
        den_JK = real(mean(Q2_4[array_idx]))/mean(N4_weight[array_idx])

        v4_Psi4 = nan_to_num(sqrt(mean(Q4_2[array_idx])
                                  /mean(N2_weight[array_idx])))

        chi_422_JK[iev] = num_JK/den_JK
        v422_JK[iev] = nan_to_num(num_JK/sqrt(den_JK))
        rho422_JK[iev] = v422_JK[iev]/v4_Psi4
        v4L_JK[iev] = nan_to_num(sqrt(v4_Psi4**2 - v422_JK[iev]**2.))

    chi_422_mean = mean(chi_422_JK)
    chi_422_err = sqrt((nev - 1.)/nev*sum((chi_422_JK - chi_422_mean)**2.))
    v422_mean = mean(v422_JK)
    v422_err = sqrt((nev - 1.)/nev*sum((v422_JK - v422_mean)**2.))
    rho422_mean = mean(rho422_JK)
    rho422_err = sqrt((nev - 1.)/nev*sum((rho422_JK - rho422_mean)**2.))
    v4L_mean = mean(v4L_JK)
    v4L_err = sqrt((nev - 1.)/nev*sum((v4L_JK - v4L_mean)**2.))
    return(v4L_mean, v4L_err, v422_mean, v422_err, rho422_mean, rho422_err,
           chi_422_mean, chi_422_err)


def calculate_chi_523(vn_array):
    """chi_523 = Re(v5*conj(v2*v3))/|v2|^2/|v3|^2
       v_523 = Re(v5*conj(v2)*conj(v3))/sqrt(|v2|^2*|v3|^2)
       rho_523 = v_523/v5(Psi5)
       v5_L = sqrt(v5(Psi5)^2 - v5(Psi23)^2)
    """
    dN = real(vn_array[:, 0])
    Q1 = dN*vn_array[:, 1]
    Q2 = dN*vn_array[:, 2]
    Q3 = dN*vn_array[:, 3]
    Q5 = dN*vn_array[:, 5]
    nev = len(dN)

    N2_weight = dN*(dN - 1.)
    Q5_2 = abs(Q5)**2. - dN

    N4_weight = dN*(dN - 1.)*(dN - 2.)*(dN - 3.)
    Q_32 = ((abs(Q2)**2.)*(abs(Q3)**2.) - 2.*real(Q5*conj(Q2)*conj(Q3))
        - 2.*real(Q3*conj(Q1)*conj(Q2)) + abs(Q5)**2. + abs(Q1)**2.
        - (dN - 4.)*(abs(Q2)**2. + abs(Q3)**2.) + dN*(dN - 6.)
    )

    N3_weight = dN*(dN - 1.)*(dN - 2.)
    chi_523_num = (Q5*conj(Q2)*conj(Q3) - Q3*conj(Q3) - Q2*conj(Q2)
                   - Q5*conj(Q5) + 2.*dN)

    chi_523_JK = zeros(nev)
    v523_JK = zeros(nev)
    rho523_JK = zeros(nev)
    v5L_JK = zeros(nev)
    for iev in range(nev):
        array_idx = [True]*nev
        array_idx[iev] = False
        array_idx = array(array_idx)

        num_JK = real(mean(chi_523_num[array_idx]))/mean(N3_weight[array_idx])
        den_JK = real(mean(Q_32[array_idx]))/mean(N4_weight[array_idx])

        v5_Psi5 = nan_to_num(sqrt(mean(Q5_2[array_idx])
                             /mean(N2_weight[array_idx])))

        chi_523_JK[iev] = num_JK/den_JK
        v523_JK[iev] = nan_to_num(num_JK/sqrt(den_JK))
        rho523_JK[iev] = v523_JK[iev]/v5_Psi5
        v5L_JK[iev] = nan_to_num(sqrt(v5_Psi5**2. - v523_JK[iev]**2.))

    chi_523_mean = mean(chi_523_JK)
    chi_523_err = sqrt((nev - 1.)/nev*sum((chi_523_JK - chi_523_mean)**2.))
    v523_mean = mean(v523_JK)
    v523_err = sqrt((nev - 1.)/nev*sum((v523_JK - v523_mean)**2.))
    rho523_mean = mean(rho523_JK)
    rho523_err = sqrt((nev - 1.)/nev*sum((rho523_JK - rho523_mean)**2.))
    v5L_mean = mean(v5L_JK)
    v5L_err = sqrt((nev - 1.)/nev*sum((v5L_JK - v5L_mean)**2.))
    return(v5L_mean, v5L_err, v523_mean, v523_err, rho523_mean, rho523_err,
           chi_523_mean, chi_523_err)


def calculate_chi_6222(vn_array):
    """arXiv: 2002.00633
       chi_6222 = Re(v6*conj(v2)**3.)/|v2|^6
       chi_633 = Re(v6*conj(v3)**2.)/|v3|^4
       v6222 = Re(v6*conj(v2)**3.)/sqrt(|v2|^6)
       v633 = Re(v6*conj(v3)**2.)/sqrt(|v3|^4)
       v624 = Re(v6*conj(v2)*conj(v4).)/sqrt(|v2|^2|v4|^2)
       rho_6222 = v6(Psi2)/v6(Psi6)
       rho_633 = v6(Psi3)/v6(Psi6)
    """
    dN = real(vn_array[:, 0])
    Q2 = dN*vn_array[:, 2]
    Q3 = dN*vn_array[:, 3]
    Q4 = dN*vn_array[:, 4]
    Q6 = dN*vn_array[:, 6]
    nev = len(dN)

    N6_weight = dN*(dN - 1.)*(dN - 2.)*(dN - 3.)*(dN - 4.)*(dN - 5.)
    Q2_6 = (abs(Q2)**6. + 9*(abs(Q4)**2.)*(abs(Q2)**2.)
            - 6.*real(Q4*Q2*conj(Q2)*conj(Q2)*conj(Q2))
            + 4.*real(Q6*conj(Q2)*conj(Q2)*conj(Q2))
            - 12.*real(Q6*conj(Q4)*conj(Q2))
            + 18.*(dN - 4.)*real(Q4*conj(Q2)*conj(Q2))
            + 4.*(abs(Q6)**2.)
            - 9.*(dN - 4.)*((abs(Q2)**4.) + (abs(Q4)**2.))
            + 18.*(dN - 5.)*(dN - 2.)*(abs(Q2)**2.)
            - 6.*dN*(dN - 4.)*(dN - 5.))

    N4_weight = dN*(dN - 1.)*(dN - 2.)*(dN - 3.)
    Q2_4 = ((abs(Q2)**4.) - 2.*real(Q4*conj(Q2)*conj(Q2))
             - 4.*(dN - 2.)*(abs(Q2)**2.) + abs(Q4)**2.
             + 2*dN*(dN - 3.))
    Q3_4 = ((abs(Q3)**4.) - 2.*real(Q6*conj(Q3)*conj(Q3))
             - 4.*(dN - 2.)*(abs(Q3)**2.) + abs(Q6)**2.
             + 2*dN*(dN - 3.))
    Q_42 = ((abs(Q2)**2.)*(abs(Q4)**2.) - 2.*real(Q6*conj(Q2)*conj(Q4))
        - 2.*real(Q4*conj(Q2)*conj(Q2)) + abs(Q6)**2. + abs(Q2)**2.
        - (dN - 4.)*(abs(Q2)**2. + abs(Q4)**2.) + dN*(dN - 6.)
    )
    chi_6222_num = (Q6*conj(Q2)*conj(Q2)*conj(Q2) - 3.*Q6*conj(Q4)*conj(Q2)
                    - 3.*Q4*conj(Q2)*conj(Q2) + 2.*Q6*conj(Q6) + 6.*Q2*conj(Q2)
                    + 3.*Q4*conj(Q4) - 6.*dN)

    N3_weight = dN*(dN - 1.)*(dN - 2.)
    chi_633_num = Q6*conj(Q3)*conj(Q3) - 2.*Q3*conj(Q3) - Q6*conj(Q6) + 2.*dN
    v624_num = (Q6*conj(Q2)*conj(Q4) - Q4*conj(Q4) - Q2*conj(Q2)
                - Q6*conj(Q6) + 2.*dN)
    Q422 = Q4*conj(Q2)*conj(Q2) - 2.*Q2*conj(Q2) - Q4*conj(Q4) + 2.*dN

    N2_weight = dN*(dN - 1.)
    Q2_2 = abs(Q2)**2. - dN
    Q4_2 = abs(Q4)**2. - dN
    Q6_2 = abs(Q6)**2. - dN

    chi_6222_JK = zeros(nev)
    chi_633_JK = zeros(nev)
    v6222_JK = zeros(nev)
    v633_JK = zeros(nev)
    rho6222_JK = zeros(nev)
    rho633_JK = zeros(nev)
    v624_JK = zeros(nev)
    chi_624_JK = zeros(nev)
    for iev in range(nev):
        array_idx = [True]*nev
        array_idx[iev] = False
        array_idx = array(array_idx)

        num_JK1 = (real(mean(chi_6222_num[array_idx]))
                   /mean(N4_weight[array_idx]))
        den_JK1 = real(mean(Q2_6[array_idx]))/mean(N6_weight[array_idx])
        num_JK2 = real(mean(chi_633_num[array_idx]))/mean(N3_weight[array_idx])
        den_JK2 = real(mean(Q3_4[array_idx]))/mean(N4_weight[array_idx])
        num_JK3 = real(mean(v624_num[array_idx]))/mean(N3_weight[array_idx])
        v22 = real(mean(Q2_2[array_idx]))/mean(N2_weight[array_idx])
        v24 = real(mean(Q2_4[array_idx]))/mean(N4_weight[array_idx])
        v422 = real(mean(Q422[array_idx]))/mean(N3_weight[array_idx])
        v42 = real(mean(Q4_2[array_idx]))/mean(N2_weight[array_idx])
        v6_Psi6 = nan_to_num(sqrt(mean(Q6_2[array_idx])
                             /mean(N2_weight[array_idx])))

        chi_6222_JK[iev] = num_JK1/den_JK1
        v6222_JK[iev] = nan_to_num(num_JK1/sqrt(den_JK1))
        rho6222_JK[iev] = nan_to_num(v6222_JK[iev]/v6_Psi6)
        chi_633_JK[iev] = num_JK2/den_JK2
        v633_JK[iev] = nan_to_num(num_JK2/sqrt(den_JK2))
        rho633_JK[iev] = nan_to_num(v633_JK[iev]/v6_Psi6)
        chi_624_JK[iev] = real((num_JK3*v24 - num_JK1*v422)
                               /((v24*v42 - v422**2)*v22))
        v624_JK[iev] = (
            num_JK3/sqrt(real(mean(Q_42[array_idx]))/mean(N4_weight[array_idx]))
        )

    chi_6222_mean = mean(chi_6222_JK)
    chi_6222_err = sqrt((nev - 1.)/nev*sum((chi_6222_JK - chi_6222_mean)**2.))
    chi_633_mean = mean(chi_633_JK)
    chi_633_err = sqrt((nev - 1.)/nev*sum((chi_633_JK - chi_633_mean)**2.))
    chi_624_mean = mean(chi_624_JK)
    chi_624_err = sqrt((nev - 1.)/nev*sum((chi_624_JK - chi_624_mean)**2.))
    v6222_mean = mean(v6222_JK)
    v6222_err = sqrt((nev - 1.)/nev*sum((v6222_JK - v6222_mean)**2.))
    v633_mean = mean(v633_JK)
    v633_err = sqrt((nev - 1.)/nev*sum((v633_JK - v633_mean)**2.))
    v624_mean = mean(v624_JK)
    v624_err = sqrt((nev - 1.)/nev*sum((v624_JK - v624_mean)**2.))
    rho6222_mean = mean(rho6222_JK)
    rho6222_err = sqrt((nev - 1.)/nev*sum((rho6222_JK - rho6222_mean)**2.))
    rho633_mean = mean(rho633_JK)
    rho633_err = sqrt((nev - 1.)/nev*sum((rho633_JK - rho633_mean)**2.))
    return(v6222_mean, v6222_err, v633_mean, v633_err, v624_mean, v624_err,
           rho6222_mean, rho6222_err, rho633_mean, rho633_err,
           chi_6222_mean, chi_6222_err, chi_633_mean, chi_633_err,
           chi_624_mean, chi_624_err)


def calculate_chi_7223(vn_array):
    """chi_7223 = Re(v7*conj(v2*v2*v3))/|v2|^4/|v3|^2
    """
    dN = real(vn_array[:, 0])
    Q2 = dN*vn_array[:, 2]
    Q3 = dN*vn_array[:, 3]
    Q4 = dN*vn_array[:, 4]
    Q5 = dN*vn_array[:, 5]
    Q7 = dN*vn_array[:, 7]
    nev = len(dN)

    N4_weight = dN*(dN - 1.)*(dN - 2.)*(dN - 3.)
    Q_7223 = (Q7*conj(Q2)*conj(Q2)*conj(Q3) - 2.*Q5*conj(Q2)*conj(Q3)
              - Q4*conj(Q2)*conj(Q2) - Q7*conj(Q4)*conj(Q3)
              - 2.*Q7*conj(Q2)*conj(Q5)
              + 2*(abs(Q7)**2. + abs(Q5)**2 + abs(Q3)**2) + 4*abs(Q2)**2
              + abs(Q4)**2 - 6*dN
    )
    Q2_4 = ((abs(Q2)**4.) - 2.*real(Q4*conj(Q2)*conj(Q2))
             - 4.*(dN - 2.)*(abs(Q2)**2.) + abs(Q4)**2.
             + 2*dN*(dN - 3.))

    N2_weight = dN*(dN - 1.)
    Q3_2 = abs(Q3)**2. - dN
    Q7_2 = abs(Q7)**2. - dN

    chi_7223_JK = zeros(nev)
    v7223_JK = zeros(nev)
    rho7223_JK = zeros(nev)
    for iev in range(nev):
        array_idx = [True]*nev
        array_idx[iev] = False
        array_idx = array(array_idx)

        num_JK = real(mean(Q_7223[array_idx]))/mean(N4_weight[array_idx])
        den_JK = real(mean(Q2_4[array_idx]*Q3_2[array_idx])
                      /mean(N4_weight[array_idx]*N2_weight[array_idx]))
        v7_Psi7 = nan_to_num(sqrt(mean(Q7_2[array_idx])
                             /mean(N2_weight[array_idx])))

        chi_7223_JK[iev] = num_JK/den_JK
        v7223_JK[iev] = nan_to_num(num_JK/sqrt(den_JK))
        rho7223_JK[iev] = nan_to_num(v7223_JK[iev]/v7_Psi7)

    chi_7223_mean = mean(chi_7223_JK)
    chi_7223_err = sqrt((nev - 1.)/nev*sum((chi_7223_JK - chi_7223_mean)**2.))
    v7223_mean = mean(v7223_JK)
    v7223_err = sqrt((nev - 1.)/nev*sum((v7223_JK - v7223_mean)**2.))
    rho7223_mean = mean(rho7223_JK)
    rho7223_err = sqrt((nev - 1.)/nev*sum((rho7223_JK - rho7223_mean)**2.))
    return(v7223_mean, v7223_err, rho7223_mean, rho7223_err,
           chi_7223_mean, chi_7223_err)


def calculate_v6_L(chi_6222, chi_6222_err, chi_633, chi_633_err, vn_array):
    """
        v6_L = sqrt(v6(Psi6)^2 - chi_6222^2 v2^6
                    - chi_633^2 v3^4 - 2 Re(chi_6222*chi_633*v2^3 v3^{2*}))
    """
    dN = real(vn_array[:, 0])
    v2_array = vn_array[:, 2]
    v3_array = vn_array[:, 3]
    v6_array = vn_array[:, 6]
    nev = len(dN)

    v6_Psi6_sq = mean(abs(v6_array)**2.)
    v6_Psi6_sq_err = std(abs(v6_array)**2.)/sqrt(nev)
    v2_6 = mean(abs(v2_array)**6.)
    v2_6_err = std(abs(v2_array)**6.)/sqrt(nev)
    v3_4 = mean(abs(v3_array)**4.)
    v3_4_err = std(abs(v3_array)**4.)/sqrt(nev)
    v23 = real(mean(v2_array**3.*conj(v3_array)**2.))
    v23_err = real(std(v2_array**3.*conj(v3_array)**2.))/sqrt(nev)
    v6_L = (v6_Psi6_sq - chi_6222**2.*v2_6 - chi_633**2.*v3_4
            - 2.*chi_6222*chi_633*v23)
    v6_L_err = sqrt(
            v6_Psi6_sq_err**2.
            + (2.*chi_6222*chi_6222_err*v2_6)**2. + (chi_6222**2.*v2_6_err)**2.
            + (2.*chi_633*chi_633_err*v3_4)**2. + (chi_633**2.*v3_4_err)**2.
            + (2.*chi_6222_err*chi_633*v23)**2.
            + (2.*chi_6222*chi_633_err*v23)**2.
            + (2.*chi_6222*chi_633*v23_err)**2.)
    return(v6_L, v6_L_err)


def calculate_nonlinear_reponse(vn_array, outputFileName):
    """
        this function computes all the nonlinear response coefficients
        proposed in the paper arXiv: 1502.02502 up to v6
    """
    v4coef = calculate_chi_422(vn_array)
    v5coef = calculate_chi_523(vn_array)
    v6coef= calculate_chi_6222(vn_array)
    #v6_L, v6_L_err = calculate_v6_L(chi_6222, chi_6222_err,
    #                                chi_633, chi_633_err, vn_array)
    v7coef = calculate_chi_7223(vn_array)

    results = list(v4coef) + list(v5coef) + list(v6coef) + list(v7coef)
    f = open(outputFileName, 'w')
    f.write("# type  value  stat. err\n")
    for i in range(len(nonlinear_reponse_correlator_name_list)):
        f.write("%s  %.10e  %.10e\n"
                % (nonlinear_reponse_correlator_name_list[i],
                   results[2*i], results[2*i+1]))
    f.close()
    return


def calcualte_vn_2(vn_data_array):
    """
        this function computes vn{2} and its stat. err.
        self correlation is substracted
    """
    vn_data_array = array(vn_data_array)
    nev = len(vn_data_array[:, 0])
    dN = real(vn_data_array[:, 0])
    dN = dN.reshape(len(dN), 1)
    Qn_array = dN*vn_data_array[:, 1:]
    corr = 1./(dN*(dN - 1.))*(Qn_array*conj(Qn_array) - dN)
    vn_2 = sqrt(real(mean(corr, 0))) + 1e-30
    vn_2_err = std(real(corr), 0)/sqrt(nev)/2./vn_2
    return(nan_to_num(vn_2), nan_to_num(vn_2_err))


def calcualte_vn_2_with_gap(vn_data_array_sub1, vn_data_array_sub2):
    """
        this function computes vn{2} and its stat. err.
        using two subevents with a eta gap
    """
    vn_data_array_sub1 = array(vn_data_array_sub1)
    vn_data_array_sub2 = array(vn_data_array_sub2)
    nev = len(vn_data_array_sub1[:, 0])
    dN1 = real(vn_data_array_sub1[:, 0])
    dN1 = dN1.reshape(len(dN1), 1)
    dN2 = real(vn_data_array_sub2[:, 0])
    dN2 = dN1.reshape(len(dN2), 1)
    Qn_array1 = dN1*vn_data_array_sub1[:, 1:]
    Qn_array2 = dN2*vn_data_array_sub2[:, 1:]

    corr = (Qn_array1*conj(Qn_array2))/(dN1*dN2)
    vn_2 = sqrt(real(mean(corr, 0))) + 1e-30
    vn_2_err = std(real(corr), 0)/sqrt(nev)/2./vn_2
    return(nan_to_num(vn_2), nan_to_num(vn_2_err))


def get_vn_diff_2PC_from_single_event(data):
    """This function computes the 2PC vn for a single event"""
    dN_event = data[:, -1]
    temp_vn_real_array = []
    temp_vn_imag_array = []
    temp_vn_denorm_array = []
    for iorder in range(1, n_order):
        vn_real_event = data[:, 2*iorder]
        vn_imag_event = data[:, 2*iorder+1]
        vn_pt = vn_real_event + 1j*vn_imag_event
        numerator_real = real(dN_event*vn_pt)
        numerator_imag = imag(dN_event*vn_pt)
        denorm = dN_event
        temp_vn_real_array.append(numerator_real)
        temp_vn_imag_array.append(numerator_imag)
    temp_vn_denorm_array.append(denorm)
    return(temp_vn_real_array, temp_vn_imag_array, temp_vn_denorm_array)


def calculate_diff_vn_single_event(pT_ref_low, pT_ref_high, data, data_ref):
    """
        This function computes pT differential vn{4} for a single event
        It returns [Qn_pT_arr, Qn_ref_arr]
    """
    npT = 50
    pT_inte_array = linspace(pT_ref_low, pT_ref_high, npT)
    dpT = pT_inte_array[1] - pT_inte_array[0]
    dN_event = data[:, -1]
    dN_ref_event = data_ref[:, -1]
    pT_ref_event = data_ref[:, 0]
    dN_ref_interp = exp(interp(pT_inte_array, pT_ref_event,
                               log(dN_ref_event + 1e-30)))
    dN_ref = sum(dN_ref_interp)*dpT/0.1
    temp_Qn_pT_array = [dN_event,]
    temp_Qn_ref_array = [dN_ref]
    for iorder in range(1, n_order):
        vn_real_event = data[:, 2*iorder]
        vn_imag_event = data[:, 2*iorder+1]
        vn_ref_real_event = data_ref[:, 2*iorder]
        vn_ref_imag_event = data_ref[:, 2*iorder+1]
        vn_ref_real_interp = interp(pT_inte_array, pT_ref_event,
                                    vn_ref_real_event)
        vn_ref_imag_interp = interp(pT_inte_array, pT_ref_event,
                                    vn_ref_imag_event)
        vn_ref_real_inte = (
            sum(vn_ref_real_interp*dN_ref_interp)/sum(dN_ref_interp))
        vn_ref_imag_inte = (
            sum(vn_ref_imag_interp*dN_ref_interp)/sum(dN_ref_interp))
        Qn_ref = dN_ref*(vn_ref_real_inte + 1j*vn_ref_imag_inte)
        Qn_pt = dN_event*(vn_real_event + 1j*vn_imag_event)
        temp_Qn_pT_array.append(Qn_pt)
        temp_Qn_ref_array.append(Qn_ref)
    return(temp_Qn_pT_array, temp_Qn_ref_array)


def calculate_vn_diff_SP(QnpT_diff, Qnref):
    """
        this funciton calculates the scalar-product vn
        assumption: no overlap between particles of interest
                    and reference flow Qn vectors
        inputs: QnpT_diff[nev, norder, npT], Qnref[nev, norder]
        return: [vn{SP}(pT), vn{SP}(pT)_err]
    """
    QnpT_diff = array(QnpT_diff)
    Qnref = array(Qnref)
    nev, norder, npT = QnpT_diff.shape

    vn_diff_SP = []
    Nref = real(Qnref[:, 0])
    N2refPairs = Nref*(Nref - 1.)
    NpTPOI = real(QnpT_diff[:, 0, :])
    N2POIPairs = NpTPOI*Nref.reshape(nev, 1)
    for iorder in range(1, norder):
        # compute Cn^ref{2}
        QnRef_tmp = Qnref[:, iorder]
        n2ref = abs(QnRef_tmp)**2. - Nref

        # compute vn{SP}(pT)
        QnpT_tmp = QnpT_diff[:, iorder, :]
        n2pT = real(QnpT_tmp*conj(QnRef_tmp.reshape(nev, 1)))

        # calcualte observables with Jackknife resampling method
        vnSPpT_arr = zeros([nev, npT])
        for iev in range(nev):
            array_idx = [True]*nev
            array_idx[iev] = False
            array_idx = array(array_idx)

            Cn2ref_arr = mean(n2ref[array_idx])/mean(N2refPairs[array_idx])
            vnSPpT_arr[iev, :] = (mean(n2pT[array_idx], 0)
                    /mean(N2POIPairs[array_idx], 0)/sqrt(Cn2ref_arr))
        vnSPpT_mean = mean(vnSPpT_arr, 0)
        vnSPpT_err  = sqrt((nev - 1.)/nev
                           *sum((vnSPpT_arr - vnSPpT_mean)**2., 0))
        vn_diff_SP.append(vnSPpT_mean)
        vn_diff_SP.append(vnSPpT_err)
    return vn_diff_SP


def calculate_vn_diff_2PC(vn_diff_real, vn_diff_imag, vn_diff_denorm):
    """this funciton calculates the rms vn[2](pT)"""
    vn_diff_real = array(vn_diff_real)
    vn_diff_imag = array(vn_diff_imag)
    vn_diff_denorm = array(vn_diff_denorm)
    nev = len(vn_diff_denorm[:, 0])
    vn_diff_2PC = sqrt(
        mean((vn_diff_real**2. + vn_diff_imag**2. - vn_diff_denorm)
             /(vn_diff_denorm**2. - vn_diff_denorm + 1e-15), 0))
    vn_diff_2PC_err = (
        std((vn_diff_real**2. + vn_diff_imag**2. - vn_diff_denorm)
            /(vn_diff_denorm**2. - vn_diff_denorm + 1e-15), 0)
        /sqrt(nev)/(2.*vn_diff_2PC + 1e-15))
    return(nan_to_num(vn_diff_2PC), nan_to_num(vn_diff_2PC_err))


def calculate_vn4_diff(QnpT_diff, Qnref):
    """
        this funciton calculates the 4-particle vn(pT)
        assumption: no overlap between particles of interest
                    and reference flow Qn vectors
        inputs: QnpT_diff[nev, norder, npT], Qnref[nev, norder]
        return: [v1{4}pT, v1{4}pT_err, v2{4}pT, v2{4}pT_err,
                 v3{4}pT, v3{4}pT_err]
    """
    QnpT_diff = array(QnpT_diff)
    Qnref = array(Qnref)
    nev, norder, npT = QnpT_diff.shape

    vn4pT_arr = []
    for iorder in range(1, 4):
        # compute Cn^ref{4}
        Nref = real(Qnref[:, 0])
        QnRef_tmp = Qnref[:, iorder]
        Q2nRef_tmp = Qnref[:, 2*iorder]
        N4refPairs = Nref*(Nref - 1.)*(Nref - 2.)*(Nref - 3.)
        n4ref = (abs(QnRef_tmp)**4.
                 - 2.*real(Q2nRef_tmp*conj(QnRef_tmp)*conj(QnRef_tmp))
                 - 4.*(Nref - 2)*abs(QnRef_tmp)**2. + abs(Q2nRef_tmp)**2.
                 + 2.*Nref*(Nref - 3))
        N2refPairs = Nref*(Nref - 1.)
        n2ref = abs(QnRef_tmp)**2. - Nref

        # compute dn{4}(pT)
        NpTPOI = real(QnpT_diff[:, 0, :])
        QnpT_tmp = QnpT_diff[:, iorder, :]
        Nref = Nref.reshape(len(Nref), 1)
        QnRef_tmp = QnRef_tmp.reshape(len(QnRef_tmp), 1)
        Q2nRef_tmp = Q2nRef_tmp.reshape(len(Q2nRef_tmp), 1)
        N4POIPairs = NpTPOI*(Nref - 1.)*(Nref - 2.)*(Nref - 3.) + 1e-30
        n4pT = real(QnpT_tmp*QnRef_tmp*conj(QnRef_tmp)*conj(QnRef_tmp)
                    - 2.*(Nref - 1)*QnpT_tmp*conj(QnRef_tmp)
                    - QnpT_tmp*QnRef_tmp*conj(Q2nRef_tmp))
        N2POIPairs = NpTPOI*Nref + 1e-30
        n2pT = real(QnpT_tmp*conj(QnRef_tmp))

        # calcualte observables with Jackknife resampling method
        Cn2ref_arr = zeros(nev)
        Cn4ref_arr = zeros(nev)
        dn4pT_arr = zeros(npT)
        vn4pT4_arr = zeros([nev, npT])
        for iev in range(nev):
            array_idx = [True]*nev
            array_idx[iev] = False
            array_idx = array(array_idx)

            Cn2ref_arr[iev] = (
                    mean(n2ref[array_idx])/mean(N2refPairs[array_idx]))
            Cn4ref_arr[iev] = (
                mean(n4ref[array_idx])/mean(N4refPairs[array_idx])
                - 2.*(Cn2ref_arr[iev])**2.)

            dn4pT_arr = (
                mean(n4pT[array_idx, :], 0)/mean(N4POIPairs[array_idx, :], 0)
                - 2.*mean(n2pT[array_idx, :], 0)
                    /mean(N2POIPairs[array_idx, :], 0)
                    *Cn2ref_arr[iev]
            )

            vn4pT4_arr[iev, :] = (-dn4pT_arr)**4./((-Cn4ref_arr[iev])**3.)

        vn4pT4_mean = mean(vn4pT4_arr, axis=0)
        vn4pT4_err  = sqrt((nev - 1.)/nev
                            *sum((vn4pT4_arr - vn4pT4_mean)**2., axis=0))

        vn4pT     = zeros(npT)
        vn4pT_err = zeros(npT)
        idx = vn4pT4_mean > 0
        vn4pT[idx] = vn4pT4_mean[idx]**(0.25)
        vn4pT_err[idx] = vn4pT4_err[idx]/(4.*vn4pT4_mean[idx]**(0.75))
        vn4pT_arr.append(vn4pT)
        vn4pT_arr.append(vn4pT_err)
    return(vn4pT_arr)



def calculate_vn_distribution(vn_array, outputFileName):
    """This function computes the vn distribution"""
    nbin = 20
    vn_dim = len(vn_array[0, :])
    output = []
    for vn_order in range(vn_dim):
        vn_mag_array = abs(vn_array[:, vn_order])
        vn_min = min(vn_mag_array)
        vn_max = max(vn_mag_array)*1.0001
        bin_boundaries = linspace(vn_min, vn_max, nbin+1)
        bin_width = bin_boundaries[1] - bin_boundaries[0]
        bin_center = zeros([nbin])
        bin_value = zeros([nbin])
        for vn_elem in vn_mag_array:
            vn_idx = int(floor((vn_elem - vn_min)/bin_width))
            if (vn_idx == 20):
                print(vn_elem, vn_min, bin_width)
            bin_value[vn_idx] += 1.
            bin_center[vn_idx] += vn_elem
        bin_center = bin_center/(bin_value + 1e-15)
        bin_value = bin_value/len(vn_array)
        bin_value_err = sqrt(bin_value/len(vn_array))
        bin_value = bin_value/bin_width
        bin_value_err = bin_value_err/bin_width
        for i in range(nbin):
            if abs(bin_center[i]) < 1e-15:
                bin_center[i] = (bin_boundaries[i] + bin_boundaries[i+1])/2.
        output.append(bin_center)
        output.append(bin_value)
        output.append(bin_value_err)
    output = array(output).transpose()

    f = open(outputFileName, 'w')
    f.write("#vn  dP(vn)/dvn  dP(vn)/dvn_err\n")
    for ipT in range(len(output[:, 0])):
        for iorder in range(1, n_order):
            f.write("%.10e  %.10e  %.10e  "
                    % (output[ipT, 3*(iorder-1)], 
                       output[ipT, 3*(iorder-1)+1],
                       output[ipT, 3*(iorder-1)+2]))
        f.write("\n")
    f.close()
    return


def calcualte_event_plane_correlations_3sub(vn_array, vn_array_sub1,
                                            vn_array_sub2, outputFileName):
    """
        this function compute the three-particle correlations with Qn
        vectors from three different sub-events
    """
    vn_array = array(vn_array)
    vn_array_sub1 = array(vn_array_sub1)
    vn_array_sub2 = array(vn_array_sub2)
    nev = len(vn_array[:, 0])

    dN = real(vn_array[:, 0].reshape(nev, 1))
    dN_sub1 = real(vn_array_sub1[:, 0].reshape(nev, 1))
    dN_sub2 = real(vn_array_sub2[:, 0].reshape(nev, 1))
    Qn_array = dN*vn_array
    Qn_array_sub1 = dN_sub1*vn_array_sub1
    Qn_array_sub2 = dN_sub2*vn_array_sub2

    corr_224_JK = zeros(nev)
    corr_336_JK = zeros(nev)
    corr_235_JK = zeros(nev)
    corr_246_JK = zeros(nev)
    for iev in range(nev):
        array_idx = [True]*nev
        array_idx[iev] = False
        array_idx = array(array_idx)

        v2_2 = (mean(real(Qn_array_sub1[array_idx, 2]
                          *conj(Qn_array_sub2[array_idx, 2])))
                /mean(dN_sub1[array_idx]*dN_sub2[array_idx]))
        v3_2 = (mean(real(Qn_array_sub1[array_idx, 3]
                          *conj(Qn_array_sub2[array_idx, 3])))
                /mean(dN_sub1[array_idx]*dN_sub2[array_idx]))
        v4_2 = (mean(real(Qn_array_sub1[array_idx, 4]
                          *conj(Qn_array_sub2[array_idx, 4])))
                /mean(dN_sub1[array_idx]*dN_sub2[array_idx]))
        v5_2 = (mean(real(Qn_array_sub1[array_idx, 5]
                          *conj(Qn_array_sub2[array_idx, 5])))
                /mean(dN_sub1[array_idx]*dN_sub2[array_idx]))
        v6_2 = (mean(real(Qn_array_sub1[array_idx, 6]
                          *conj(Qn_array_sub2[array_idx, 6])))
                /mean(dN_sub1[array_idx]*dN_sub2[array_idx]))

        # cos(4(Psi_2 - Psi_4))
        corr_224_num = (
            mean(real(  Qn_array[array_idx, 2]*Qn_array_sub1[array_idx, 2]
                        *conj(Qn_array_sub2[array_idx, 4])
                      + Qn_array[array_idx, 2]*Qn_array_sub2[array_idx, 2]
                        *conj(Qn_array_sub1[array_idx, 4])
                      + Qn_array_sub2[array_idx, 2]*Qn_array_sub1[array_idx, 2]
                        *conj(Qn_array[array_idx, 4])))
            /mean(3.*dN[array_idx]*dN_sub1[array_idx]*dN_sub2[array_idx])
        )
        corr_224_JK[iev] = corr_224_num/sqrt(v2_2*v2_2*v4_2)

        # cos(6(Psi_3 - Psi_6))
        corr_336_num = (
            mean(real(  Qn_array[array_idx, 3]*Qn_array_sub1[array_idx, 3]
                        *conj(Qn_array_sub2[array_idx, 6])
                      + Qn_array[array_idx, 3]*Qn_array_sub2[array_idx, 3]
                        *conj(Qn_array_sub1[array_idx, 6])
                      + Qn_array_sub2[array_idx, 3]*Qn_array_sub1[array_idx, 3]
                        *conj(Qn_array[array_idx, 6])))
            /mean(3.*dN[array_idx]*dN_sub1[array_idx]*dN_sub2[array_idx])
        )
        corr_336_JK[iev] = corr_336_num/sqrt((v3_2**2.)*v6_2)

        # cos(2Psi_2 + 3Psi_3 - 5Psi_5)
        corr_235_num = (
            mean(real(  Qn_array[array_idx, 2]*Qn_array_sub1[array_idx, 3]
                        *conj(Qn_array_sub2[array_idx, 5])
                      + Qn_array[array_idx, 3]*Qn_array_sub2[array_idx, 2]
                        *conj(Qn_array_sub1[array_idx, 5])
                      + Qn_array_sub1[array_idx, 2]*Qn_array_sub2[array_idx, 3]
                        *conj(Qn_array[array_idx, 5])
                      + Qn_array[array_idx, 2]*Qn_array_sub2[array_idx, 3]
                        *conj(Qn_array_sub1[array_idx, 5])
                      + Qn_array[array_idx, 3]*Qn_array_sub1[array_idx, 2]
                        *conj(Qn_array_sub2[array_idx, 5])
                      + Qn_array_sub1[array_idx, 3]*Qn_array_sub2[array_idx, 2]
                        *conj(Qn_array[array_idx, 5])
            ))
            /mean(6.*dN[array_idx]*dN_sub1[array_idx]*dN_sub2[array_idx])
        )
        corr_235_JK[iev] = corr_235_num/sqrt(v2_2*v3_2*v5_2)

        # cos(2Psi_2 + 4Psi_4 - 6Psi_6)
        corr_246_num = (
            mean(real(  Qn_array[array_idx, 2]*Qn_array_sub1[array_idx, 4]
                        *conj(Qn_array_sub2[array_idx, 6])
                      + Qn_array[array_idx, 4]*Qn_array_sub2[array_idx, 2]
                        *conj(Qn_array_sub1[array_idx, 6])
                      + Qn_array_sub1[array_idx, 2]*Qn_array_sub2[array_idx, 4]
                        *conj(Qn_array[array_idx, 6])
                      + Qn_array[array_idx, 2]*Qn_array_sub2[array_idx, 4]
                        *conj(Qn_array_sub1[array_idx, 6])
                      + Qn_array[array_idx, 4]*Qn_array_sub1[array_idx, 2]
                        *conj(Qn_array_sub2[array_idx, 6])
                      + Qn_array_sub1[array_idx, 4]*Qn_array_sub2[array_idx, 2]
                        *conj(Qn_array[array_idx, 6])
            ))
            /mean(6.*dN[array_idx]*dN_sub1[array_idx]*dN_sub2[array_idx])
        )
        corr_246_JK[iev] = corr_246_num/sqrt(v2_2*v4_2*v6_2)

    corr_224 = mean(corr_224_JK)
    corr_224_err = sqrt((nev - 1.)/nev*sum((corr_224_JK - corr_224)**2.))
    corr_336 = mean(corr_336_JK)
    corr_336_err = sqrt((nev - 1.)/nev*sum((corr_336_JK - corr_336)**2.))
    corr_235 = mean(corr_235_JK)
    corr_235_err = sqrt((nev - 1.)/nev*sum((corr_235_JK - corr_235)**2.))
    corr_246 = mean(corr_246_JK)
    corr_246_err = sqrt((nev - 1.)/nev*sum((corr_246_JK - corr_246)**2.))

    # output results to a file
    f = open(outputFileName, 'w')
    f.write("#correlator  value  value_err\n")
    f.write("224  %.5e  %.5e\n" % (corr_224, corr_224_err))
    f.write("336  %.5e  %.5e\n" % (corr_336, corr_336_err))
    f.write("235  %.5e  %.5e\n" % (corr_235, corr_235_err))
    f.write("246  %.5e  %.5e\n" % (corr_246, corr_246_err))
    f.close()
    return


def calcualte_event_plane_correlations(vn_array, outputFileName):
    """
        this function compute the scalar-product event plane correlations
        vn_array is a matrix [event_idx, vn_order]
    """
    nev = len(vn_array[:, 0])
    v2_array = vn_array[:, 2]
    v3_array = vn_array[:, 3]
    v4_array = vn_array[:, 4]
    v5_array = vn_array[:, 5]
    v6_array = vn_array[:, 6]

    corr_224_JK = zeros(nev)
    corr_22233_JK = zeros(nev)
    corr_2226_JK = zeros(nev)
    corr_336_JK = zeros(nev)
    corr_235_JK = zeros(nev)
    corr_246_JK = zeros(nev)
    corr_234_JK = zeros(nev)
    for iev in range(nev):
        array_idx = [True]*nev
        array_idx[iev] = False
        array_idx = array(array_idx)

        v2_2 = mean(abs(v2_array[array_idx])**2.)
        v3_2 = mean(abs(v3_array[array_idx])**2.)
        v4_2 = mean(abs(v4_array[array_idx])**2.)
        v5_2 = mean(abs(v5_array[array_idx])**2.)
        v6_2 = mean(abs(v6_array[array_idx])**2.)

        # cos(4(Psi_2 - Psi_4))
        corr_224_num = mean(real((v2_array[array_idx]**2.)
                                 *conj(v4_array[array_idx])))
        corr_224_JK[iev] = corr_224_num/sqrt(v2_2*v2_2*v4_2)

        # cos(6(Psi_2 - Psi_3))
        corr_22233_num = mean(real((v2_array[array_idx]**3.)
                                   *conj(v3_array[array_idx])**2.))
        corr_22233_JK[iev] = corr_22233_num/sqrt(v2_2**3.*v3_2**2.)

        # cos(6(Psi_2 - Psi_6))
        corr_2226_num = mean(real(v2_array[array_idx]**3.
                                  *conj(v6_array[array_idx])))
        corr_2226_JK[iev] = corr_2226_num/sqrt((v2_2**3.)*v6_2)

        # cos(6(Psi_3 - Psi_6))
        corr_336_num = mean(real((v3_array[array_idx]**2.)
                                 *conj(v6_array[array_idx])))
        corr_336_JK[iev] = corr_336_num/sqrt((v3_2**2.)*v6_2)

        # cos(2Psi_2 + 3Psi_3 - 5Psi_5)
        corr_235_num = mean(real(v2_array[array_idx]*v3_array[array_idx]
                                 *conj(v5_array[array_idx])))
        corr_235_JK[iev] = corr_235_num/sqrt(v2_2*v3_2*v5_2)

        # cos(2Psi_2 + 4Psi_4 - 6Psi_6)
        corr_246_num = mean(real(v2_array[array_idx]*v4_array[array_idx]
                                 *conj(v6_array[array_idx])))
        corr_246_JK[iev] = corr_246_num/sqrt(v2_2*v4_2*v6_2)

        # cos(2Psi_2 - 6Psi_3 + 4Psi_4)
        corr_234_num = mean(real(v2_array[array_idx]
                                 *(conj(v3_array[array_idx])**2.)
                                 *v4_array[array_idx]))
        corr_234_JK[iev] = corr_234_num/sqrt(v2_2*(v3_2**2.)*v4_2)

    corr_224 = mean(corr_224_JK)
    corr_224_err = sqrt((nev - 1.)/nev*sum((corr_224_JK - corr_224)**2.))
    corr_22233 = mean(corr_22233_JK)
    corr_22233_err = sqrt((nev - 1.)/nev*sum((corr_22233_JK - corr_22233)**2.))
    corr_2226 = mean(corr_2226_JK)
    corr_2226_err = sqrt((nev - 1.)/nev*sum((corr_2226_JK - corr_2226)**2.))
    corr_336 = mean(corr_336_JK)
    corr_336_err = sqrt((nev - 1.)/nev*sum((corr_336_JK - corr_336)**2.))
    corr_235 = mean(corr_235_JK)
    corr_235_err = sqrt((nev - 1.)/nev*sum((corr_235_JK - corr_235)**2.))
    corr_246 = mean(corr_246_JK)
    corr_246_err = sqrt((nev - 1.)/nev*sum((corr_246_JK - corr_246)**2.))
    corr_234 = mean(corr_234_JK)
    corr_234_err = sqrt((nev - 1.)/nev*sum((corr_234_JK - corr_234)**2.))

    f = open(outputFileName, 'w')
    f.write("#correlator  value  value_err\n")
    f.write("4(24)  %.5e  %.5e\n" % (corr_224, corr_224_err))
    f.write("6(23)  %.5e  %.5e\n" % (corr_22233, corr_22233_err))
    f.write("6(26)  %.5e  %.5e\n" % (corr_2226, corr_2226_err))
    f.write("6(36)  %.5e  %.5e\n" % (corr_336, corr_336_err))
    f.write("(235)  %.5e  %.5e\n" % (corr_235, corr_235_err))
    f.write("(246)  %.5e  %.5e\n" % (corr_246, corr_246_err))
    f.write("(234)  %.5e  %.5e\n" % (corr_234, corr_234_err))
    f.close()
    return


def calculate_vn_arrays_for_rn_ratios(data):
    """
        this function compute the complex pT-integrated Vn vector
        in different pT ranges for a single event
        it returns a 2d matrix vn_arrays[pT_idx, n_order_idx]
    """
    pT_boundaries = [0.3, 0.5, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0]
    npT = 50
    vn_arrays = []
    for ipT in range(len(pT_boundaries)-1):
        pT_low = pT_boundaries[ipT]
        pT_high = pT_boundaries[ipT + 1]
        pT_mid = (pT_low + pT_high)/2.
        vn_array = calcualte_inte_vn(pT_low, pT_high, data)
        vn_array.insert(0, pT_mid)
        vn_arrays.append(vn_array)
    return(vn_arrays)


def calculate_rn_ratios(vn_event_arrays, avg_folder):
    """
        this function compute rn ratio in different pT bins
        according to the CMS measurements
        it reads in a 3d data cube
              vn_event_arrays[event_idx, pT_idx, n_order_idx]
        it returns rn_arrays[iorder, pT_idx, 3]
    """
    vn_event_arrays = array(vn_event_arrays)
    rn_arrays = []
    for iorder in range(3, 6):
        # compute r2, r3, r4
        rn_array = []
        for itrig in range(3, len(vn_event_arrays[0, :, 0])):
            pT_trig = real(vn_event_arrays[0, itrig, 0])
            dN_trig = real(vn_event_arrays[:, itrig, 1])
            Qn_trig_array = dN_trig*vn_event_arrays[:, itrig, iorder]
            nev = len(Qn_trig_array)

            denorm2_dN = dN_trig*(dN_trig - 1.)
            denorm2_array = abs(Qn_trig_array)**2. - dN_trig

            for iasso in range(0, itrig+1):
                pT_asso = real(vn_event_arrays[0, iasso, 0])
                dN_asso = real(vn_event_arrays[:, iasso, 1])
                Qn_asso_array = dN_asso*vn_event_arrays[:, iasso, iorder]

                num_dN = dN_trig*dN_asso
                num_array = real(Qn_asso_array*conj(Qn_trig_array))
                if iasso == itrig:
                    num_dN -= dN_asso
                    num_array = (real(Qn_asso_array*conj(Qn_trig_array))
                                 - dN_asso)

                denorm1_dN = dN_asso*(dN_asso - 1.)
                denorm1_array = abs(Qn_asso_array)**2. - dN_asso

                rn_jackknife = zeros(nev)
                for iev in range(nev):
                    array_idx = [True]*nev
                    array_idx[iev] = False
                    array_idx = array(array_idx)

                    num = mean(num_array[array_idx])/mean(num_dN[array_idx])
                    denorm1 = (mean(denorm1_array[array_idx])
                               /mean(denorm1_dN[array_idx]))
                    denorm2 = (mean(denorm2_array[array_idx])
                               /mean(denorm2_dN[array_idx]))

                    if denorm1 > 0. and denorm2 > 0.:
                        rn_jackknife[iev] = num/sqrt(denorm1*denorm2)

                rn_mean = mean(rn_jackknife)
                rn_err = sqrt((nev - 1.)/nev*sum((rn_jackknife - rn_mean)**2.))
                rn_array.append([pT_trig - pT_asso, rn_mean, rn_err])
        rn_arrays.append(rn_array)
    rn_arrays = array(rn_arrays)

    # output rn ratios
    pT_trig = ['1.0', '1.5', '2.0', '2.5', '3.0']
    ipTtrig = 0
    output_filename = ("%s_rn_ratios_CMS_pTtrig_%s_%s.dat"
                       % (particle_name_list[ipart],
                          pT_trig[ipTtrig], pT_trig[ipTtrig+1]))
    f = open(path.join(avg_folder, output_filename), 'w')
    f.write("#pT_mid  rn  rn_err (n = 2, 3, 4)\n")
    for ipT in range(len(rn_arrays[0, :, 0])):
        for iorder in range(len(rn_arrays[:, 0, 0])):
            f.write("%.5e  %.5e  %.5e  "
                    % (rn_arrays[iorder, ipT, 0],
                       rn_arrays[iorder, ipT, 1],
                       rn_arrays[iorder, ipT, 2]))
        f.write("\n")
        if rn_arrays[0, ipT, 0] == 0.0:
            f.close()
            ipTtrig += 1
            if ipTtrig < (len(pT_trig) - 1):
                output_filename = ("%s_rn_ratios_CMS_pTtrig_%s_%s.dat"
                                   % (particle_name_list[ipart],
                                      pT_trig[ipTtrig],
                                      pT_trig[ipTtrig+1]))
                f = open(path.join(avg_folder, output_filename), 'w')
                f.write("#pT_mid  rn  rn_err (n = 2, 3, 4)\n")
    return


def calculate_symmetric_cumulant(vn_data_array, outputFileName):
    """
        this funciton computes the symmetric cumulant
            SC(m,n) = <v_m*conj(v_m)*v_n*conj(v_n)>
                      - <v_m*conj(v_m)>*<v_n*conj(v_n)>
        we use Jackknife resampling method to estimate the statistical error
    """
    nev = len(vn_data_array[:, 0])
    dN = real(vn_data_array[:, 0])
    Q1 = dN*vn_data_array[:, 1]
    Q2 = dN*vn_data_array[:, 2]
    Q3 = dN*vn_data_array[:, 3]
    Q4 = dN*vn_data_array[:, 4]
    Q5 = dN*vn_data_array[:, 5]
    Q6 = dN*vn_data_array[:, 6]

    # two-particle correlation
    N2_weight = dN*(dN - 1.)
    Q2_2 = abs(Q2)**2. - dN
    Q3_2 = abs(Q3)**2. - dN
    Q4_2 = abs(Q4)**2. - dN

    # four-particle correlation
    N4_weight = dN*(dN - 1.)*(dN - 2.)*(dN - 3.)
    Q_32 = ((abs(Q2)**2.)*(abs(Q3)**2.) - 2.*real(Q5*conj(Q2)*conj(Q3))
        - 2.*real(Q3*conj(Q1)*conj(Q2)) + abs(Q5)**2. + abs(Q1)**2.
        - (dN - 4.)*(abs(Q2)**2. + abs(Q3)**2.) + dN*(dN - 6.)
    )
    Q_42 = ((abs(Q2)**2.)*(abs(Q4)**2.) - 2.*real(Q6*conj(Q2)*conj(Q4))
        - 2.*real(Q4*conj(Q2)*conj(Q2)) + abs(Q6)**2. + abs(Q2)**2.
        - (dN - 4.)*(abs(Q2)**2. + abs(Q4)**2.) + dN*(dN - 6.)
    )

    # calcualte observables with Jackknife resampling method
    SC32_array = zeros(nev)
    SC42_array = zeros(nev)
    NSC32_array = zeros(nev)
    NSC42_array = zeros(nev)
    for iev in range(nev):
        array_idx = [True]*nev
        array_idx[iev] = False
        array_idx = array(array_idx)

        # SC(3,2)
        v2v3 = ((mean(Q3_2[array_idx])*mean(Q2_2[array_idx]))
                /(mean(N2_weight[array_idx])**2.))
        SC32_array[iev] = (
                mean(Q_32[array_idx])/mean(N4_weight[array_idx]) - v2v3)
        NSC32_array[iev] = SC32_array[iev]/v2v3

        # SC(4,2)
        v2v4 = ((mean(Q4_2[array_idx])*mean(Q2_2[array_idx]))
                /(mean(N2_weight[array_idx])**2.))
        SC42_array[iev] = (
                mean(Q_42[array_idx])/mean(N4_weight[array_idx]) - v2v4)
        NSC42_array[iev] = SC42_array[iev]/v2v4

    SC32_mean = mean(SC32_array)
    SC32_err = sqrt((nev - 1.)/nev*sum((SC32_array - SC32_mean)**2.))
    NSC32_mean = mean(NSC32_array)
    NSC32_err = sqrt((nev - 1.)/nev*sum((NSC32_array - NSC32_mean)**2.))

    SC42_mean = mean(SC42_array)
    SC42_err = sqrt((nev - 1.)/nev*sum((SC42_array - SC42_mean)**2.))
    NSC42_mean = mean(NSC42_array)
    NSC42_err = sqrt((nev - 1.)/nev*sum((NSC42_array - NSC42_mean)**2.))

    results = [SC32_mean, SC32_err, NSC32_mean, NSC32_err,
               SC42_mean, SC42_err, NSC42_mean, NSC42_err]
    f = open(outputFileName, 'w')
    f.write("# type  SC{mn}  SC{mn}_err  NSC{mn}  NSC{mn}_err\n")
    for i in range(len(symmetric_cumulant_name_list)):
        f.write("%s  %.10e  %.10e  %.10e  %.10e\n"
                % (symmetric_cumulant_name_list[i],
                   results[4*i],   results[4*i+1],
                   results[4*i+2], results[4*i+3]))
    f.close()
    return


def calculate_vn4_vn6(vn_data_array, outputFileName_vn4,
                      outputFileName42, outputFileName64):
    """
        this funciton computes the 4 particle cumulant vn{4}
            vn{4} = (2 <v_n*conj(v_n)>**2 - <(v_n*conj(v_n))**2.>)**(1/4)

        this funciton also computes the ratio of
        the 4-particle cumulant vn{4} over the 2-particle cumulant vn{2}
        and Fn = sqrt((vn{2}^2 - vn{4}^2)/(vn{2}^2 + vn{4}^2))

            vn{4} = (2 <v_n*conj(v_n)>**2 - <(v_n*conj(v_n))**2.>)**(1/4)
            vn{2} = (<v_n*conj(v_n)>)**(1/2)

        this funciton also computes the ratio of
        the 6-particle cumulant vn{6} over the 4-particle cumulant vn{4}
            cn{6} = <<6>> - 9<<2>><<4>> + 12<<2>>^3
            vn{6} = (cn{6}/4)**(1/6)
            vn{4} = (2 <v_n*conj(v_n)>**2 - <(v_n*conj(v_n))**2.>)**(1/4)
        and compute skewness estimator gamma_1
            gamma_1 = -6\sqrt{2}*vn{4}^2*(vn{4} - vn{6})
                                         /(vn{2}^2 - vn{4}^2)^(3/2)

        we will use Jackknife resampling method to estimate
        the statistical error
    """
    nev = len(vn_data_array[:, 0])
    dN = real(vn_data_array[:, 0])
    Q1 = dN*vn_data_array[:, 1]
    Q2 = dN*vn_data_array[:, 2]
    Q3 = dN*vn_data_array[:, 3]
    Q4 = dN*vn_data_array[:, 4]
    Q5 = dN*vn_data_array[:, 5]
    Q6 = dN*vn_data_array[:, 6]

    # two-particle correlation
    N2_weight = dN*(dN - 1.)
    #Q1_2 = abs(Q1)**2. - dN
    Q2_2 = abs(Q2)**2. - dN
    Q3_2 = abs(Q3)**2. - dN

    # four-particle correlation
    N4_weight = dN*(dN - 1.)*(dN - 2.)*(dN - 3.)
    #Q1_4 = ((abs(Q1)**4.) - 2.*real(Q2*conj(Q1)*conj(Q1))
    #         - 4.*(dN - 2.)*(abs(Q1)**2.) + abs(Q2)**2.
    #         + 2*dN*(dN - 3.))
    Q2_4 = ((abs(Q2)**4.) - 2.*real(Q4*conj(Q2)*conj(Q2))
             - 4.*(dN - 2.)*(abs(Q2)**2.) + abs(Q4)**2.
             + 2*dN*(dN - 3.))
    Q3_4 = ((abs(Q3)**4.) - 2.*real(Q6*conj(Q3)*conj(Q3))
             - 4.*(dN - 2.)*(abs(Q3)**2.) + abs(Q6)**2.
             + 2*dN*(dN - 3.))

    # six-particle correlation
    N6_weight = dN*(dN - 1.)*(dN - 2.)*(dN - 3.)*(dN - 4.)*(dN - 5.)
    Q2_6 = (abs(Q2)**6. + 9*(abs(Q4)**2.)*(abs(Q2)**2.)
            - 6.*real(Q4*Q2*conj(Q2)*conj(Q2)*conj(Q2))
            + 4.*real(Q6*conj(Q2)*conj(Q2)*conj(Q2))
            - 12.*real(Q6*conj(Q4)*conj(Q2))
            + 18.*(dN - 4.)*real(Q4*conj(Q2)*conj(Q2))
            + 4.*(abs(Q6)**2.)
            - 9.*(dN - 4.)*((abs(Q2)**4.) + (abs(Q4)**2.))
            + 18.*(dN - 5.)*(dN - 2.)*(abs(Q2)**2.)
            - 6.*dN*(dN - 4.)*(dN - 5.))

    # calcualte observables with Jackknife resampling method
    C1_4_array = zeros(nev)
    C2_4_array = zeros(nev)
    C3_4_array = zeros(nev)
    r1_array = zeros(nev)
    r2_array = zeros(nev)
    r3_array = zeros(nev)
    F1_array = zeros(nev)
    F2_array = zeros(nev)
    F3_array = zeros(nev)
    r26_array = zeros(nev)
    gamma1_array = zeros(nev)
    for iev in range(nev):
        array_idx = [True]*nev
        array_idx[iev] = False
        array_idx = array(array_idx)

        # C_n{4}
        #C_1_4 = (mean(Q1_4[array_idx])/mean(N4_weight[array_idx])
        #         - 2.*((mean(Q1_2[array_idx])/mean(N2_weight[array_idx]))**2.))
        #C_1_2 = mean(Q1_2[array_idx])/mean(N2_weight[array_idx])
        #if C_1_4 < 0. and C_1_2 > 0.:
        #    v1_4 = (-C_1_4)**0.25
        #    v1_2 = sqrt(C_1_2)
        #    r1_array[iev] = v1_4/(v1_2 + 1e-15)
        #    F1_array[iev] = sqrt((v1_2**2. - v1_4**2.)
        #                         /(v1_2**2. + v1_4**2. + 1e-15))

        C_2_2 = mean(Q2_2[array_idx])/mean(N2_weight[array_idx])
        C24_tmp = mean(Q2_4[array_idx])/mean(N4_weight[array_idx])
        C_2_4 = C24_tmp - 2.*C_2_2**2.
        C_2_6 = (mean(Q2_6[array_idx])/mean(N6_weight[array_idx])
                 - 9.*C_2_2*C24_tmp + 12.*(C_2_2**3.))
        C2_4_array[iev] = C_2_4
        if C_2_4 < 0. and C_2_2 > 0.:
            v2_4 = (-C_2_4)**0.25
            v2_2 = sqrt(C_2_2)
            r2_array[iev] = v2_4/v2_2
            F2_array[iev] = sqrt((v2_2**2. - v2_4**2.)
                                 /(v2_2**2. + v2_4**2. + 1e-15))
            if C_2_6 > 0.:
                v2_6 = (C_2_6/4.)**(1./6.)
                r26_array[iev] = v2_6/v2_4
                gamma1_array[iev] = (-6.*sqrt(2)*(v2_4**2.)*(v2_4 - v2_6)
                                     /(v2_2**2. - v2_4**2.)**(1.5))


        C_3_2 = mean(Q3_2[array_idx])/mean(N2_weight[array_idx])
        C34_tmp = mean(Q3_4[array_idx])/mean(N4_weight[array_idx])
        C_3_4 = C34_tmp - 2.*C_3_2**2.
        C3_4_array[iev] = C_3_4
        if C_3_4 < 0. and C_3_2 > 0.:
            v3_4 = (-C_3_4)**0.25
            v3_2 = sqrt(C_3_2)
            r3_array[iev] = v3_4/v3_2
            F3_array[iev] = sqrt((v3_2**2. - v3_4**2.)
                                 /(v3_2**2. + v3_4**2. + 1e-15))

    C1_4_mean = mean(C1_4_array)
    C1_4_err  = sqrt((nev - 1.)/nev*sum((C1_4_array - C1_4_mean)**2.))
    C2_4_mean = mean(C2_4_array)
    C2_4_err  = sqrt((nev - 1.)/nev*sum((C2_4_array - C2_4_mean)**2.))
    C3_4_mean = mean(C3_4_array)
    C3_4_err  = sqrt((nev - 1.)/nev*sum((C3_4_array - C3_4_mean)**2.))

    v1_4 = 0.0
    v1_4_err = 0.0
    if C1_4_mean < 0:
        v1_4 = (-C1_4_mean)**0.25
        v1_4_err = 0.25*((-C1_4_mean)**(-0.75))*C1_4_err

    v2_4 = 0.0
    v2_4_err = 0.0
    if C2_4_mean < 0:
        v2_4 = (-C2_4_mean)**0.25
        v2_4_err = 0.25*((-C2_4_mean)**(-0.75))*C2_4_err

    v3_4 = 0.0
    v3_4_err = 0.0
    if C3_4_mean < 0:
        v3_4 = (-C3_4_mean)**0.25
        v3_4_err = 0.25*((-C3_4_mean)**(-0.75))*C3_4_err
    results = [v1_4, v1_4_err, C1_4_mean, C1_4_err,
               v2_4, v2_4_err, C2_4_mean, C2_4_err,
               v3_4, v3_4_err, C3_4_mean, C3_4_err,]
    f = open(outputFileName_vn4, 'w')
    f.write("# n  vn{4}  vn{4}_err  Cn{4}  Cn{4}_err\n")
    for i in range(1, 4):
        f.write("%d  %.10e  %.10e  %.10e  %.10e\n"
                % (i, results[4*i-4], results[4*i-3],
                   results[4*i-2], results[4*i-1]))
    f.close()

    # now the ratios
    r1_mean = mean(r1_array)
    r1_err = sqrt((nev - 1.)/nev*sum((r1_array - r1_mean)**2.))
    r2_mean = mean(r2_array)
    r2_err = sqrt((nev - 1.)/nev*sum((r2_array - r2_mean)**2.))
    r3_mean = mean(r3_array)
    r3_err = sqrt((nev - 1.)/nev*sum((r3_array - r3_mean)**2.))

    F1_mean = mean(F1_array)
    F1_err = sqrt((nev - 1.)/nev*sum((F1_array - F1_mean)**2.))
    F2_mean = mean(F2_array)
    F2_err = sqrt((nev - 1.)/nev*sum((F2_array - F2_mean)**2.))
    F3_mean = mean(F3_array)
    F3_err = sqrt((nev - 1.)/nev*sum((F3_array - F3_mean)**2.))

    results = [r1_mean, r1_err, F1_mean, F1_err,
               r2_mean, r2_err, F2_mean, F2_err,
               r3_mean, r3_err, F3_mean, F3_err]
    f = open(outputFileName42, 'w')
    f.write("# n  vn{4}/vn{2}  (vn{4}/vn{2})_err  Fn  Fn_err\n")
    f.write("# Fn = sqrt((vn{2}^2 - vn{4}^2)/(vn{2}^2 + vn{4}^2))\n")
    for i in range(1, 4):
        f.write("%d  %.10e  %.10e  %.10e  %.10e\n"
                % (i, results[4*i - 4], results[4*i - 3],
                   results[4*i - 2], results[4*i-1]))
    f.close()
    r26_mean = mean(r26_array)
    r26_err = sqrt((nev - 1.)/nev*sum((r26_array - r26_mean)**2.))
    gamma1_mean = mean(gamma1_array)
    gamma1_err = sqrt((nev - 1.)/nev*sum((gamma1_array - gamma1_mean)**2.))

    f = open(outputFileName64, 'w')
    f.write(
        "# n  vn{6}/vn{4}  (vn{6}/vn{4})_err  gamma_1  gamma_1_err\n")
    f.write("%d  %.10e  %.10e  %.10e  %.10e\n"
            % (2, r26_mean, r26_err, gamma1_mean, gamma1_err))
    f.close()
    return


def calculate_vn_eta(eta_array, Qn_rap_array, eta_min, eta_max):
    """
        This function computes vn(eta).
        eta_min and eta_max specify the rapidity range of reference flow vector
    """
    nev, nvn, neta = Qn_rap_array[:, 1:, :].shape
    dN_array = Qn_rap_array[:, 0, :].reshape((nev, 1, neta))
    idx = (eta_array > eta_min) & (eta_array < eta_max)
    Qn_ref = sum(Qn_rap_array[:, 1:, idx], axis=2)
    Qn_ref = Qn_ref.reshape((nev, nvn, 1))
    vn_SP_ev    = real(Qn_rap_array[:, 1:, :]*conj(Qn_ref))
    vn_SP_array = zeros([nev, nvn, neta])
    for iev in range(nev):
        array_idx      = [True]*nev
        array_idx[iev] = False
        array_idx      = array(array_idx)
        vn_den         = mean((absolute(Qn_ref[array_idx, :, :]))**2., axis=0)
        vn_SP          = (mean(vn_SP_ev[array_idx, :, :], axis=0)
                          /((sqrt(vn_den) + 1e-16)
                             *mean(dN_array[array_idx, :, :], axis=0)))
        vn_SP_array[iev, :, :] = vn_SP
    vn_SP_mean = mean(vn_SP_array, axis=0)
    vn_SP_err  = sqrt((nev - 1.)/nev*sum((vn_SP_array - vn_SP_mean)**2., axis=0))
    return([vn_SP_mean, vn_SP_err])


def calculate_rn_eta(eta_array, Qn_rap_array, outputFileName):
    """
        This function computes the longitudinal factorization breaking ratios
        for all n passed from vn_array
            eta, rn(eta)
    """
    nev, nQn, neta = Qn_rap_array[:, 1:, :].shape
    Qn_array = Qn_rap_array[:, 1:, :]

    # calculate the reference flow vector for every event
    eta_b_min    = 2.5
    eta_b_max    = 4.0
    eta_ref1_tmp = linspace(eta_b_min, eta_b_max, 16)
    eta_ref2_tmp = linspace(-eta_b_max, -eta_b_min, 16)
    Qn_ref1      = []
    Qn_ref2      = []
    for iev in range(nev):
        Qn_ref1_vec = []
        Qn_ref2_vec = []
        for iorder in range(nQn):
            Qn1_interp = interp(eta_ref1_tmp, eta_array,
                                Qn_array[iev, iorder, :])
            Qn2_interp = interp(eta_ref2_tmp, eta_array,
                                Qn_array[iev, iorder, :])
            Qn_ref1_vec.append(sum(Qn1_interp))
            Qn_ref2_vec.append(sum(Qn2_interp))
        Qn_ref1.append(Qn_ref1_vec)
        Qn_ref2.append(Qn_ref2_vec)
    Qn_ref1 = array(Qn_ref1).reshape((nev, nQn, 1))
    Qn_ref2 = array(Qn_ref2).reshape((nev, nQn, 1))

    rn_num    = real(Qn_array[:, :, ::-1]*conj(Qn_ref1))
    rn_den    = real(Qn_array*conj(Qn_ref1))
    rnn_num   = real((Qn_ref2*conj(Qn_array))
                     *(Qn_array[:, :, ::-1]*conj(Qn_ref1)))
    rnn_den   = real((Qn_ref2*conj(Qn_array[:, :, ::-1]))
                     *(Qn_array*conj(Qn_ref1)))

    # compute the error using jack-knife
    rn_array  = zeros([nev, nQn, neta])
    rnn_array = zeros([nev, nQn, neta])
    for iev in range(nev):
        array_idx      = [True]*nev
        array_idx[iev] = False
        array_idx      = array(array_idx)
        rn_ev          = (mean(rn_num[array_idx], axis=0)
                          /(mean(rn_den[array_idx], axis=0) + 1e-15))
        rnn_ev         = (mean(rnn_num[array_idx], axis=0)
                          /(mean(rnn_den[array_idx], axis=0) + 1e-15))
        rn_array[iev, :, :]  = rn_ev
        rnn_array[iev, :, :] = rnn_ev
    rn_mean  = mean(rn_array, axis=0)
    rn_err   = sqrt((nev - 1.)/nev*sum((rn_array - rn_mean)**2., axis=0))
    rnn_mean = mean(rnn_array, axis=0)
    rnn_err  = sqrt((nev - 1.)/nev*sum((rnn_array - rnn_mean)**2., axis=0))

    f = open(outputFileName, 'w')
    f.write("#eta  rn(eta)  rn_err(eta)  rnn(eta)  rnn_err(eta)\n")
    for ieta in range(len(eta_array)-1):
        f.write("%.10e  " % eta_array[ieta])
        for iorder in range(0, n_order-1):
            f.write("%.10e  %.10e  %.10e  %.10e  "
                    % (rn_mean[iorder, ieta], rn_err[iorder, ieta],
                       rnn_mean[iorder, ieta], rnn_err[iorder, ieta]))
        f.write("\n")
    f.close()
    return


def calculate_meanpT_fluc(dN_array, pT_array, pT_min=0.0, pT_max=3.0):
    """
        This function computes the mean pT fluctuations
            returns sigma_pT/<pT>, sigma_pT/<pT>_err
        here sigma_pT is the standard deviation of the event-by-event mean pT
        This function accepts pT_cut through [pT_min, pT_max]
        dN_array is dN/(2\pi dy pT dpT)
    """
    npT_interp = 50
    pT_inte_array = linspace(pT_min, pT_max, npT_interp)

    nev, npT = dN_array.shape
    mean_pT_array = zeros(nev)
    for iev in range(nev):
        dN_interp = exp(interp(pT_inte_array, pT_array[iev, :],
                               log(dN_array[iev, :] + 1e-30)))
        mean_pT_array[iev] = (sum(pT_inte_array**2.*dN_interp)
                              /sum(pT_inte_array*dN_interp))

    # compute the error using jack-knife
    rn_array  = zeros(nev)
    for iev in range(nev):
        array_idx      = [True]*nev
        array_idx[iev] = False
        array_idx      = array(array_idx)
        rn_ev          = (std(mean_pT_array[array_idx])
                          /(mean(mean_pT_array[array_idx]) + 1e-15))
        rn_array[iev]  = rn_ev
    rn_mean  = mean(rn_array, axis=0)
    rn_err   = sqrt((nev - 1.)/nev*sum((rn_array - rn_mean)**2.))
    return([rn_mean, rn_err])


hf = h5py.File(path.join(data_path, data_name), "r")
event_list = list(hf.keys())
print("total number of events: {}".format(len(event_list)))

dNdyDict = {}
dNdyList = []
if CentralityFlag > 0:
    for ifolder, event_name in enumerate(event_list):
        file_name = "particle_9999_vndata_eta_-0.5_0.5.dat"
        if RapidityTrigger == 1:      # PHENIX BBC Trigger
            file_name = "particle_9999_vndata_eta_-3.9_-3.1.dat"
        event_group = hf.get(event_name)
        #eventStatus = check_an_event_is_good(event_group)
        eventStatus = True
        if eventStatus:
            temp_data   = event_group.get(file_name)
            temp_data   = nan_to_num(temp_data)
            dNdyDict[event_name] = temp_data[0, 1]
    dNdyList = -sort(-array(list(dNdyDict.values())))
print("Number of good events: {}".format(len(dNdyList)))


for icen in range(len(centralityCutList) - 1):
    if centralityCutList[icen+1] < centralityCutList[icen]: continue
    avg_folder = path.join(
        avg_folder_header, "{0:02.0f}-{1:02.0f}".format(
            centralityCutList[icen], centralityCutList[icen+1])
    )

    if path.isdir(avg_folder):
        print("{} already exists, skipped ...".format(avg_folder))
        continue
    else:
        mkdir(avg_folder)

    selected_events_list = []
    if CentralityFlag == 0:
        cen_label = ("C{0:d}-{1:d}_".format(
            int(centralityCutList[icen]),
            int(centralityCutList[icen+1]))
        )
        for ifolder, event_name in enumerate(event_list):
            if cen_label in event_name:
                selected_events_list.append(event_name)
    elif CentralityFlag > 0:
        dN_dy_cut_high = dNdyList[
            int(len(dNdyList)*centralityCutList[icen]/100.)
        ]
        dN_dy_cut_low  = dNdyList[
            min(len(dNdyList)-1,
                int(len(dNdyList)*centralityCutList[icen+1]/100.))
        ]
        if len(dNcutList) == len(centralityCutList):
            dN_dy_cut_high = dNcutList[icen]
            dN_dy_cut_low = dNcutList[icen+1]

        for event_name in dNdyDict.keys():
            if (dNdyDict[event_name] > dN_dy_cut_low
                and dNdyDict[event_name] <= dN_dy_cut_high):
                selected_events_list.append(event_name)

    nev = len(selected_events_list)
    print("analysis {}%-{}% nev = {}...".format(
            centralityCutList[icen], centralityCutList[icen+1], nev))
    print("dNdy: {0:.2f} - {1:.2f}".format(dN_dy_cut_low, dN_dy_cut_high))
    if nev == 0:
        print("Skip ...")
        continue

    for ipart, particle_id in enumerate(particle_list):
        print("processing %s ..." % particle_name_list[ipart])

        # first particle yield dN/dy
        if particle_id == '9999':
            n_order = 10
            file_name = 'particle_9999_vndata_eta_-0.5_0.5.dat'
        else:
            file_name = 'particle_%s_vndata_y_-0.5_0.5.dat' % particle_id
            n_order = 6

        dN_dy = zeros(len(selected_events_list))
        for ifolder, event_name in enumerate(selected_events_list):
            event_group = hf.get(event_name)
            temp_data = event_group.get(file_name)
            temp_data = nan_to_num(temp_data)
            dN_dy[ifolder] = temp_data[0, 1]
        dN_dy_avg = mean(dN_dy)
        dN_dy_avg_err = std(dN_dy)/sqrt(nev)

        # then <pT>, vn, dN/(2pi dy pT dpT), vn{SP}(pT)
        if particle_id == '9999':
            file_name = 'particle_9999_vndata_diff_eta_-0.5_0.5.dat'
        else:
            file_name = f'particle_{particle_id}_vndata_diff_y_-0.5_0.5.dat'
        file_name_ref1 = 'particle_9999_vndata_diff_eta_0.5_1.dat'
        file_name_ref2 = 'particle_9999_vndata_diff_eta_-1_-0.5.dat'

        pT_array = []
        dN_array = []
        vn_array_list = []
        vn_array_sub1_list = []; vn_array_sub2_list = []
        QnpT_diff_list = []; Qnref_list = []
        for iExp in kinematicCutsDict.keys():
            vn_array_list.append([])
            vn_array_sub1_list.append([])
            vn_array_sub2_list.append([])
            QnpT_diff_list.append([])
            Qnref_list.append([])

        for ifolder, event_name in enumerate(selected_events_list):
            event_group = hf.get(event_name)
            temp_data = event_group.get(file_name)
            temp_data = nan_to_num(temp_data)
            temp_data_ref1 = event_group.get(file_name_ref1)
            temp_data_ref1 = nan_to_num(temp_data_ref1)
            temp_data_ref2 = event_group.get(file_name_ref2)
            temp_data_ref2 = nan_to_num(temp_data_ref2)

            dN_event = temp_data[:, 1]  # dN/(2pi dy pT dpT)
            pT_event = temp_data[:, 0]

            # record particle spectra
            pT_array.append(pT_event)
            dN_array.append(dN_event)

            for iExp, expKey in enumerate(kinematicCutsDict.keys()):
                pTmin = kinematicCutsDict[expKey]['pTmin']
                pTmax = kinematicCutsDict[expKey]['pTmax']
                temp_vn_array = calcualte_inte_vn(pTmin, pTmax, temp_data)
                vn_array_list[iExp].append(temp_vn_array)
                temp_vn_array = calcualte_inte_vn(pTmin, pTmax, temp_data_ref1)
                vn_array_sub1_list[iExp].append(temp_vn_array)
                temp_vn_array = calcualte_inte_vn(pTmin, pTmax, temp_data_ref2)
                vn_array_sub2_list[iExp].append(temp_vn_array)
                temp_arr = calculate_diff_vn_single_event(
                                    pTmin, pTmax, temp_data, temp_data_ref1)
                QnpT_diff_list[iExp].append(temp_arr[0])
                Qnref_list[iExp].append(temp_arr[1])

        # now we perform event average
        dN_array = array(dN_array)
        pT_array = array(pT_array)

        n_pT = len(pT_array[0, :])
        pT_spectra = zeros([n_pT])
        for ipT in range(len(pT_array[0, :])):
            dN_temp = sum(dN_array[:, ipT]*pT_array[:, ipT])
            if(dN_temp > 0):
                pT_spectra[ipT] = (
                        sum(pT_array[:, ipT]**2.*dN_array[:, ipT])/dN_temp)
            else:
                pT_spectra[ipT] = mean(pT_array[:, ipT])
        # dN/(2pi dy pT dpT)
        dN_spectra = mean(pT_array*dN_array, 0)/pT_spectra
        dN_spectra_err = std(pT_array*dN_array, 0)/pT_spectra/sqrt(nev)

        # calculate mean pT from event-averaged particle spectrum
        pT_interp = linspace(0.05, 2.95, 30)
        dN_interp = exp(interp(pT_interp, pT_spectra, log(dN_spectra + 1e-30)))
        dN_interp_err = interp(pT_interp, pT_spectra, dN_spectra_err)
        mean_pT = sum(pT_interp**2.*dN_interp)/sum(pT_interp*dN_interp)
        mean_pT_upper = (sum(pT_interp**2.*(dN_interp+dN_interp_err))
                         /sum(pT_interp*(dN_interp+dN_interp_err)))
        mean_pT_lower = (sum(pT_interp**2.*(dN_interp-dN_interp_err))
                         /sum(pT_interp*(dN_interp-dN_interp_err)))
        mean_pT_err = max(abs(mean_pT_upper - mean_pT), 
                          abs(mean_pT - mean_pT_lower))
        pT_interp = linspace(0.15, 2.95, 30)
        dN_interp = exp(interp(pT_interp, pT_spectra, log(dN_spectra + 1e-30)))
        dN_interp_err = interp(pT_interp, pT_spectra, dN_spectra_err)
        mean_pT_1 = sum(pT_interp**2.*dN_interp)/sum(pT_interp*dN_interp)
        mean_pT_1_upper = (sum(pT_interp**2.*(dN_interp+dN_interp_err))
                         /sum(pT_interp*(dN_interp+dN_interp_err)))
        mean_pT_1_lower = (sum(pT_interp**2.*(dN_interp-dN_interp_err))
                         /sum(pT_interp*(dN_interp-dN_interp_err)))
        mean_pT_1_err = max(abs(mean_pT_1_upper - mean_pT_1), 
                          abs(mean_pT_1 - mean_pT_1_lower))

        # loop over all kinematic cuts
        vn_2_list = []
        vn_2_gap_list = []
        for iExp, expKey in enumerate(kinematicCutsDict.keys()):
            # calcualte vn{2}
            vn_2, vn_2_err = calcualte_vn_2(vn_array_list[iExp])
            vn_2_list.append(list(zip(vn_2, vn_2_err)))

            # calculate vn{2} with |\Delta \eta| > 1
            vn_2_gap, vn_2_gap_err = calcualte_vn_2_with_gap(
                    vn_array_sub1_list[iExp], vn_array_sub2_list[iExp])
            vn_2_gap_list.append(list(zip(vn_2_gap, vn_2_gap_err)))

            # calcualte and output vn{SP}(pT)
            vndiff_SP = calculate_vn_diff_SP(QnpT_diff_list[iExp],
                                             Qnref_list[iExp])
            output_filename = "{}_differential_observables_{}.dat".format(
                                            particle_name_list[ipart], expKey)
            f = open(path.join(avg_folder, output_filename), 'w')
            f.write("#pT  dN/(2pi dy pT dpT)  dN/(2pi dy pT dpT)_err  "
                    "vn{SP}  vn{SP}_err\n")
            for ipT in range(len(pT_spectra)):
                f.write("%.10e  %.10e  %.10e  "
                        % (pT_spectra[ipT], dN_spectra[ipT],
                           dN_spectra_err[ipT]))
                for iorder in range(1, n_order):
                    f.write("%.10e  %.10e  " % (vndiff_SP[2*iorder-2][ipT],
                                                vndiff_SP[2*iorder-1][ipT]))
                f.write("\n")
            f.close()

            if particle_id == '9999':
                vn_array2 = array(vn_array_list[iExp])
                # calculate non-linear response coefficents
                output_filename = path.join(
                    avg_folder,
                    "non_linear_response_coefficients_{}.dat".format(expKey)
                )
                calculate_nonlinear_reponse(vn_array2, output_filename)

                # calculate symmetric cumulant coefficents with ALICE pT cut
                output_filename = path.join(
                    avg_folder, "symmetric_cumulant_{}.dat".format(expKey))
                calculate_symmetric_cumulant(vn_array2, output_filename)

                # calculate vn{4}, vn{4}/vn{2}, and vn{6}/vn{4}
                output_filename1 = path.join(
                    avg_folder, "charged_hadron_vn4_{}.dat".format(expKey))
                output_filename2 = path.join(
                    avg_folder,
                    "charged_hadron_vn4_over_vn2_{}.dat".format(expKey))
                output_filename3 = path.join(
                    avg_folder,
                    "charged_hadron_vn6_over_vn4_{}.dat".format(expKey))
                calculate_vn4_vn6(vn_array2, output_filename1,
                                  output_filename2, output_filename3)

                if not FastFlag:
                    # calculate vn distribution for charged hadrons
                    output_filename = path.join(
                        avg_folder,
                        "charged_hadron_vn_distribution_{}.dat".format(expKey)
                    )
                    calculate_vn_distribution(vn_array2, output_filename)

                    # calculate rn ratios
                    rn_cms = calculate_rn_ratios(vn_cms_arrays_for_rn,
                                                 avg_folder)

                    # calculate flow event-plane correlation
                    output_filename = path.join(
                        avg_folder,
                        "charged_hadron_event_plane_correlation_{}.dat".format(
                                                                        expKey)
                    )
                    calcualte_event_plane_correlations(vn_array2,
                                                       output_filename)

                    # calculate three-particle correlations
                    output_filename = path.join(avg_folder,
                        "{}_three_particle_correlation_{}.dat".format(
                                            particle_name_list[ipart], expKey)
                    )
                    calcualte_event_plane_correlations_3sub(
                        vn_array_list[iExp], vn_array_sub1_list[iExp],
                        vn_array_sub2_list[iExp], output_filename)

                    vn4_pTdiff = calculate_vn4_diff(QnpT_diff_list[iExp],
                                                    Qnref_list[iExp])
                    # output vn{4}(pT)
                    output_filename = (
                        "{}_differential_observables_4particle_{}.dat".format(
                                            particle_name_list[ipart], expKey)
                    )
                    f = open(path.join(avg_folder, output_filename), 'w')
                    f.write("#pT  vn{4}  vn{4}_err\n")
                    for ipT in range(len(pT_spectra)):
                        f.write("%.10e  " % (pT_spectra[ipT]))
                        for iorder in range(1, 4):
                            f.write("%.10e  %.10e  " % (
                                        vn4_pTdiff[2*iorder-2][ipT],
                                        vn4_pTdiff[2*iorder-1][ipT]))
                        f.write("\n")
                    f.close()

        # then particle rapidity distribution
        if particle_id == '9999':
            file_name = 'particle_%s_dNdeta_pT_0.2_3.dat' % particle_id
        else:
            file_name = 'particle_%s_dNdy_pT_0.2_3.dat' % particle_id

        eta_array = []
        dN_array = []
        Qn_rap_array = []
        for ifolder, event_name in enumerate(selected_events_list):
            event_group = hf.get(event_name)
            temp_data = event_group.get(file_name)
            temp_data = nan_to_num(temp_data)

            eta_array.append(temp_data[:, 0])
            dN_array.append(temp_data[:, 1])
            temp_Qn_array = []
            for iorder in range(0, n_order):
                Qn_real = temp_data[:, 2*iorder+1]*temp_data[:, -1]
                Qn_imag = temp_data[:, 2*iorder+2]*temp_data[:, -1]
                if iorder == 0:
                    Qn = temp_data[:, -1]
                else:
                    Qn = Qn_real + 1j*Qn_imag
                temp_Qn_array.append(Qn)
            Qn_rap_array.append(temp_Qn_array)

        eta_array = array(eta_array)
        dN_array = array(dN_array)
        Qn_rap_array = array(Qn_rap_array)

        eta_point = mean(eta_array, 0)
        dNdeta = mean(dN_array, 0)
        dNdeta_err = std(dN_array, 0)/sqrt(nev)
        if RapidityTrigger == 0:
            vn_SP_eta, vn_SP_eta_err = calculate_vn_eta(
                                        eta_point, Qn_rap_array, -1.0, 1.0)
        elif RapidityTrigger == 1:
            vn_SP_eta, vn_SP_eta_err = calculate_vn_eta(
                                        eta_point, Qn_rap_array, -3.9, -3.1)
        vn_SP_eta_mid, vn_SP_eta_mid_err = calculate_vn_eta(
                                        eta_point, Qn_rap_array, -0.5, 0.5)

        output_filename = path.join(avg_folder, "{}_rn_eta.dat".format(
                                                    particle_name_list[ipart]))
        calculate_rn_eta(eta_point, Qn_rap_array, output_filename)

        ######################################################################
        # finally, output all the results
        ######################################################################
        output_filename = ("%s_integrated_observables.dat"
                           % particle_name_list[ipart])
        f = open(path.join(avg_folder, output_filename), 'w')
        f.write("dN/dy= %.5e +/- %.5e\n" % (dN_dy_avg, dN_dy_avg_err))
        f.write("<pT>= %.5e +/- %.5e\n" % (mean_pT, mean_pT_err))
        f.write("<pT(>0.15)>= %.5e +/- %.5e\n" % (mean_pT_1, mean_pT_1_err))
        for iExp, expKey in enumerate(kinematicCutsDict.keys()):
            for iorder in range(len(vn_2_list[iExp])):
                f.write("v_%d{2}(%s)= %.5e +/- %.5e\n"
                        % (iorder+1, expKey,
                           vn_2_list[iExp][iorder][0],
                           vn_2_list[iExp][iorder][1])
                )
        f.close()
        output_filename = ("%s_integrated_observables_with_rapidity_gap.dat"
                           % particle_name_list[ipart])
        f = open(path.join(avg_folder, output_filename), 'w')
        f.write("dN/dy= %.5e +/- %.5e\n" % (dN_dy_avg, dN_dy_avg_err))
        f.write("<pT>= %.5e +/- %.5e\n" % (mean_pT, mean_pT_err))
        f.write("<pT(>0.15)>= %.5e +/- %.5e\n" % (mean_pT_1, mean_pT_1_err))
        for iExp, expKey in enumerate(kinematicCutsDict.keys()):
            for iorder in range(len(vn_2_list[iExp])):
                f.write("v_%d{2}(%s)= %.5e +/- %.5e\n"
                        % (iorder+1, expKey,
                           vn_2_gap_list[iExp][iorder][0],
                           vn_2_gap_list[iExp][iorder][1])
                )
        f.close()

        output_filename = "{}_rapidity_distribution.dat".format(
                                                    particle_name_list[ipart])
        f = open(path.join(avg_folder, output_filename), 'w')
        if(particle_id == '9999'):
            f.write("#eta  dN/deta  dN/deta_err  vn{2}(eta)  vn{2}(eta)_err"
                    + "  Re{vn}(eta) Re{vn}(eta)_err\n")
        else:
            f.write("#y  dN/dy  dN/dy_err  vn{2}(y)  vn{2}(y)_err  "
                    + "Re{vn}(y)  Re{vn}(y)_err\n")
        for ieta in range(len(eta_point)):
            f.write("%.10e  %.10e  %.10e  "
                    % (eta_point[ieta], dNdeta[ieta], dNdeta_err[ieta]))
            for iorder in range(1, n_order):
                f.write("%.10e  %.10e  %.10e  %.10e  "
                        % (vn_SP_eta[iorder-1, ieta],
                           vn_SP_eta_err[iorder-1, ieta],
                           vn_SP_eta_mid[iorder-1, ieta],
                           vn_SP_eta_mid_err[iorder-1, ieta]))
            f.write("\n")
        f.close()

print("Analysis is done.")

