#!/usr/bin/env python3
"""
     This script performs event averaging for particle
     polarization observables calculated from event-by-event simulations

     All the errors are only statistic errors
"""

from sys import argv, exit
from os import path, mkdir
from glob import glob
from numpy import *
import h5py
import shutil
from random import randrange

# define colors
purple = "\033[95m"
green = "\033[92m"
blue = "\033[94m"
yellow = "\033[93m"
red = "\033[91m"
normal = "\033[0m"

Reg_centrality_cut_list = [
    0., 5., 10., 20., 30., 40., 50., 60., 70., 80., 90., 100.
]
PHOBOS_cen_list = [0., 6., 15., 25., 35., 45., 55.]  # PHOBOS AuAu 200
SPS_cen_list = [5., 12.5, 23.5, 33.5, 43.5]  # SPS PbPb
PHENIX_cen_list = [20., 40., 60., 88.]  # PHENIX dAu
STAR_cen_list = [0., 10., 40., 80]  # STAR v1

centrality_cut_list = Reg_centrality_cut_list + [20., 50., 20., 60.]
# predefined dNcut List
#dNcutList = [1000., 544.20, 448.98, 309.78, 204.20, 131.71, 83.02, 45.85,
#             24.89, 11.73, 4.58, 0.16, 309.78, 83.02, 309.78, 45.82]
dNcutList = []

n_order = 7
vorticityType = "Thermal"
rapType = "pseudorapidity"

try:
    data_path = path.abspath(argv[1])
    data_name = data_path.split("/")[-1]
    results_folder_name = data_name.split(".h5")[0]
    avg_folder_header = path.join(path.abspath(argv[2]), results_folder_name)
    print("output folder: %s" % avg_folder_header)
    if (path.isdir(avg_folder_header)):
        print("folder %s already exists!" % avg_folder_header)
        var = input("do you want to delete it? [y/N]")
        if 'y' in var:
            shutil.rmtree(avg_folder_header)
        else:
            print("please choose another folder path~")
            exit(0)
    mkdir(avg_folder_header)
except IndexError:
    print("Usage: {} database.h5 results_folder".format(argv[0]))
    exit(1)


def check_an_event_is_good(h5_event):
    """This function checks the given event contains all required files"""
    required_files_list = [
        'particle_9999_vndata_eta_-0.5_0.5.dat',
        'particle_9999_vndata_eta_-1_-0.5.dat',
        'particle_9999_vndata_eta_0.5_1.dat',
        "Smu_pT_Thermal_{}_3122.dat".format(rapType),
        "Smu_phi_Thermal_{}_3122.dat".format(rapType),
        "Smu_y_Thermal_{}_3122.dat".format(rapType),
    ]
    event_file_list = list(h5_event.keys())
    for ifile in required_files_list:
        if ifile not in event_file_list:
            print("event {} is bad, missing {} ...".format(
                h5_event.name, ifile))
            return False
    return True


def calculate_meanpT_inte_vn(pT_low: float, pT_high: float, data):
    """
        this function calculates the dN/dy, <pT>, pT-integrated vn in a
        given pT range (pT_low, pT_high) for every event in the data
    """
    npT = 50
    pT_inte_array = linspace(pT_low, pT_high, npT)
    dpT = pT_inte_array[1] - pT_inte_array[0]
    pT_event = data[:, 0]
    dN_event = data[:, 2]
    # dN/(2pi*pT*dpT*dy)
    dN_interp = exp(interp(pT_inte_array, pT_event, log(dN_event + 1e-30)))
    N = sum(dN_interp*pT_inte_array)*dpT*2.*pi
    meanpT = sum(dN_interp*pT_inte_array**2.)/sum(dN_interp*pT_inte_array)
    temp_vn_array = [
        N,
        meanpT,
    ]
    for iorder in range(1, n_order):
        vn_real_event = data[:, 2*iorder]
        vn_imag_event = data[:, 2*iorder + 1]
        vn_real_interp = interp(pT_inte_array, pT_event, vn_real_event)
        vn_imag_interp = interp(pT_inte_array, pT_event, vn_imag_event)
        vn_real_inte = (sum(vn_real_interp*dN_interp*pT_inte_array)
                        /sum(dN_interp*pT_inte_array))
        vn_imag_inte = (sum(vn_imag_interp*dN_interp*pT_inte_array)
                        /sum(dN_interp*pT_inte_array))
        vn_inte = vn_real_inte + 1j*vn_imag_inte
        temp_vn_array.append(vn_inte)
    return (temp_vn_array)


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
    Qn_array1 = dN1*vn_data_array_sub1[:, 2:]
    Qn_array2 = dN2*vn_data_array_sub2[:, 2:]

    num = sqrt(mean(real(Qn_array1*conj(Qn_array2)), axis=0))
    num_err = std(real(Qn_array1*conj(Qn_array2)), axis=0)/sqrt(nev)/(2.*num)
    denorm = sqrt(mean(dN1*dN2))
    denorm_err = std(dN1*dN2)/sqrt(nev)/(2.*denorm)
    vn_2 = num/denorm
    vn_2_err = sqrt((num_err/denorm)**2. + (num*denorm_err/denorm**2.)**2.)
    return (nan_to_num(vn_2), nan_to_num(vn_2_err))


def calculate_diff_vn_single_event(pT_ref_low, pT_ref_high, data, data_ref):
    """
        This function computes pT differential vn{4} for a single event
        It returns [Qn_pT_arr, Qn_ref_arr]
    """
    npT = 50
    pT_inte_array = linspace(pT_ref_low, pT_ref_high, npT)
    dpT = pT_inte_array[1] - pT_inte_array[0]
    dN_event = data[:, 2]
    dN_ref_event = data_ref[:, 2]
    pT_ref_event = data_ref[:, 0]
    dN_ref_interp = exp(
        interp(pT_inte_array, pT_ref_event, log(dN_ref_event + 1e-30)))
    dN_ref = sum(dN_ref_interp*pT_inte_array)*dpT*2.*pi
    temp_Qn_pT_array = [
        dN_event,
    ]
    temp_Qn_ref_array = [
        dN_ref,
    ]
    for iorder in range(1, n_order):
        vn_real_event = data[:, 2*iorder]
        vn_imag_event = data[:, 2*iorder + 1]
        vn_ref_real_event = data_ref[:, 2*iorder]
        vn_ref_imag_event = data_ref[:, 2*iorder + 1]
        vn_ref_real_interp = interp(pT_inte_array, pT_ref_event,
                                    vn_ref_real_event)
        vn_ref_imag_interp = interp(pT_inte_array, pT_ref_event,
                                    vn_ref_imag_event)
        vn_ref_real_inte = (sum(vn_ref_real_interp*dN_ref_interp)
                            /sum(dN_ref_interp))
        vn_ref_imag_inte = (sum(vn_ref_imag_interp*dN_ref_interp)
                            /sum(dN_ref_interp))
        Qn_ref = dN_ref*(vn_ref_real_inte + 1j*vn_ref_imag_inte)
        Qn_pt = dN_event*(vn_real_event + 1j*vn_imag_event)
        temp_Qn_pT_array.append(Qn_pt)
        temp_Qn_ref_array.append(Qn_ref)
    return (temp_Qn_pT_array, temp_Qn_ref_array)


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
                                  /mean(N2POIPairs[array_idx], 0)
                                  /sqrt(Cn2ref_arr))
        vnSPpT_mean = mean(vnSPpT_arr, 0)
        vnSPpT_err = sqrt((nev - 1.)/nev*sum((vnSPpT_arr - vnSPpT_mean)**2., 0))
        vn_diff_SP.append(vnSPpT_mean)
        vn_diff_SP.append(vnSPpT_err)
    return vn_diff_SP


def calculate_vn_eta(eta_array, dN_array, vn_array, eta_min, eta_max):
    """
        This function computes vn(eta).
        eta_min and eta_max specify the rapidity range of reference flow vector
    """
    nev, neta = dN_array.shape
    dN_array = dN_array.reshape((nev, 1, neta))
    idx = (eta_array > eta_min) & (eta_array < eta_max)
    vn_ref = (sum(dN_array[:, :, idx]*vn_array[:, :, idx], axis=2)/
              (sum(dN_array[:, :, idx], axis=2) + 1e-15))
    vnshape = vn_ref.shape
    nvn = vnshape[1]
    vn_ref = vn_ref.reshape((vnshape[0], vnshape[1], 1))
    vn_SP_ev = real(vn_array*conj(vn_ref))
    vn_SP_array = zeros([nev, nvn, neta])
    for iev in range(nev):
        array_idx = [True]*nev
        array_idx[iev] = False
        array_idx = array(array_idx)
        vn_den = mean((absolute(vn_ref[array_idx, :, :]))**2., axis=0)
        vn_SP = mean(vn_SP_ev[array_idx, :, :], axis=0)/sqrt(vn_den)
        vn_SP_array[iev, :, :] = vn_SP
    vn_SP_mean = mean(vn_SP_array, axis=0)
    vn_SP_err = sqrt((nev - 1.)/nev*sum((vn_SP_array - vn_SP_mean)**2., axis=0))
    return ([vn_SP_mean, vn_SP_err])


def analyze_Smu(hf_, eventList_, pTMin_, pTMax_, outputFolder_, icen_,
                vn_array_):
    """
        This function compute the event-averaged S^mu
        pT integated from pT_min, pT_max
        y integrated from -1 to 1
    """
    filelist = [
        "Smu_pT_{}_{}_3122.dat", "Smu_pT_{}_{}_3122_wSIP_BBPP.dat",
        "Smu_pT_{}_{}_3122_wSIP_LY.dat", "Smu_pT_{}_{}_3122_wMuIP_wSIP_LY.dat"
    ]

    pT_arr = []
    dN_list = []
    Sx_list = []
    Sy_list = []
    Sz_list = []
    for ifile in filelist:
        dN_list.append([])
        Sx_list.append([])
        Sy_list.append([])
        Sz_list.append([])

    nev = len(eventList_)
    for eventName in eventList_:
        event_group = hf_.get(eventName)
        for ifile, filename in enumerate(filelist):
            spin_data = nan_to_num(
                event_group.get(filename.format(vorticityType, rapType)))

            if ifile == 0:
                pT_arr = spin_data[:, 0]
            idx = (pT_arr > pTMin_) & (pT_arr < pTMax_)

            dN_list[ifile].append(sum(spin_data[idx, 1]))
            Sx_list[ifile].append(sum(spin_data[idx, 1]*spin_data[idx, 7]))
            Sy_list[ifile].append(sum(spin_data[idx, 1]*spin_data[idx, 8]))
            Sz_list[ifile].append(sum(spin_data[idx, 1]*spin_data[idx, 9]))

    Sx_avg = []
    Sx_err = []
    Sy_avg = []
    Sy_err = []
    Sz_avg = []
    Sz_err = []
    for ifile in range(len(filelist)):
        Sx_avg.append(
            mean(array(Sx_list[ifile]), axis=0)
            /mean(array(dN_list[ifile]), axis=0))
        Sx_err.append(
            std(array(Sx_list[ifile]), axis=0)
            /mean(array(dN_list[ifile]), axis=0)/sqrt(nev))
        Sy_avg.append(
            mean(array(Sy_list[ifile]), axis=0)
            /mean(array(dN_list[ifile]), axis=0))
        Sy_err.append(
            std(array(Sy_list[ifile]), axis=0)
            /mean(array(dN_list[ifile]), axis=0)/sqrt(nev))
        Sz_avg.append(
            mean(array(Sz_list[ifile]), axis=0)
            /mean(array(dN_list[ifile]), axis=0))
        Sz_err.append(
            std(array(Sz_list[ifile]), axis=0)
            /mean(array(dN_list[ifile]), axis=0)/sqrt(nev))

    # compute correlation between S^y and v_2^2 and S^y vs. mean pT
    dN_array = real(vn_array_[:, 0])
    dN_avg = mean(dN_array)
    meanpT_array = abs(vn_array_[:, 1])
    v2_array = abs(vn_array_[:, 3])**2.

    delta_N = dN_array - dN_avg
    varN = std(dN_array)**2.
    delta_meanpT = meanpT_array - mean(meanpT_array)
    hat_delta_meanpT = delta_meanpT - mean(delta_meanpT*delta_N)/varN*delta_N
    delta_v2 = v2_array - mean(v2_array)
    hat_delta_v2 = delta_v2 - mean(delta_v2*delta_N)/varN*delta_N

    rho_pTSy_avg = []
    rho_pTSy_err = []
    rho_v2Sy_avg = []
    rho_v2Sy_err = []
    for ifile in range(len(filelist)):
        Sy_array = array(Sy_list[ifile])/array(dN_list[ifile])

        delta_Sy = Sy_array - mean(Sy_array)
        hat_delta_Sy = delta_Sy - mean(delta_Sy*delta_N)/varN*delta_N

        # compute the error using the jackknife method
        rho_v2Sy = zeros(nev)
        rho_pTSy = zeros(nev)
        for iev in range(nev):
            array_idx = [True]*nev
            array_idx[iev] = False
            array_idx = array(array_idx)

            rho_v2Sy[iev] = (
                mean(hat_delta_v2[array_idx]*hat_delta_Sy[array_idx])/sqrt(
                    mean(hat_delta_v2[array_idx]**2.)
                    *mean(hat_delta_Sy[array_idx]**2.)))
            rho_pTSy[iev] = (
                mean(hat_delta_meanpT[array_idx]*hat_delta_Sy[array_idx])/sqrt(
                    mean(hat_delta_meanpT[array_idx]**2.)
                    *mean(hat_delta_Sy[array_idx]**2.)))
        rho_v2Sy_mean_tmp = mean(rho_v2Sy)
        rho_pTSy_mean_tmp = mean(rho_pTSy)
        rho_v2Sy_avg.append(rho_v2Sy_mean_tmp)
        rho_v2Sy_err.append(
            sqrt((nev - 1.)/nev*sum((rho_v2Sy - rho_v2Sy_mean_tmp)**2.)))
        rho_pTSy_avg.append(rho_pTSy_mean_tmp)
        rho_pTSy_err.append(
            sqrt((nev - 1.)/nev*sum((rho_pTSy - rho_pTSy_mean_tmp)**2.)))

    # output results to files
    f = open(
        path.join(
            outputFolder_,
            "averaged_Smu_{}_pT_{}_{}.txt".format(vorticityType, pTMin_,
                                                  pTMax_)), "a")
    if icen_ == 0:
        f.write("# cen  Nch  S^x  S^x_err  S^y  S^y_err  S^z  S^z_err  "
                + "({0}  {0}+SIP(BBPP) {0}+SIP(LY)  {0}+SIP+MuBIP(LY))\n".
                format(vorticityType))
    f.write("{0}  {1:.5e}  ".format(
        (centrality_cut_list[icen_] + centrality_cut_list[icen_ + 1])/2.,
        dN_avg))
    for icol in range(len(filelist)):
        f.write("%.5e  %.5e  " % (Sx_avg[icol], Sx_err[icol]))
        f.write("%.5e  %.5e  " % (Sy_avg[icol], Sy_err[icol]))
        f.write("%.5e  %.5e  " % (Sz_avg[icol], Sz_err[icol]))
    f.write("\n")
    f.close()

    f = open(
        path.join(
            outputFolder_,
            "Rho_v2Sy_{}_pT_{}_{}.txt".format(vorticityType, pTMin_, pTMax_)),
        "a")
    if icen_ == 0:
        f.write("# cen  Nch  rho(v2, Sy)  rho(v2, Sy)_err  "
                + "({0}  {0}+SIP(BBPP) {0}+SIP(LY)  {0}+SIP+MuBIP(LY))\n".
                format(vorticityType))
    f.write("{0}  {1:.5e}  ".format(
        (centrality_cut_list[icen_] + centrality_cut_list[icen_ + 1])/2.,
        dN_avg))
    for icol in range(len(filelist)):
        f.write("%.5e  %.5e  " % (rho_v2Sy_avg[icol], rho_v2Sy_err[icol]))
    f.write("\n")
    f.close()

    f = open(
        path.join(
            outputFolder_,
            "Rho_pTSy_{}_pT_{}_{}.txt".format(vorticityType, pTMin_, pTMax_)),
        "a")
    if icen_ == 0:
        f.write("# cen  Nch  rho([pT], Sy)  rho([pT], Sy)_err  "
                + "({0}  {0}+SIP(BBPP) {0}+SIP(LY)  {0}+SIP+MuBIP(LY))\n".
                format(vorticityType))
    f.write("{0}  {1:.5e}  ".format(
        (centrality_cut_list[icen_] + centrality_cut_list[icen_ + 1])/2.,
        dN_avg))
    for icol in range(len(filelist)):
        f.write("%.5e  %.5e  " % (rho_pTSy_avg[icol], rho_pTSy_err[icol]))
    f.write("\n")
    f.close()


def analyze_Smu_pT(hf_, eventList_, outputFolder_):
    """
        This function compute the event-averaged S^mu(pT)
    """
    filelist = [
        "Smu_pT_{}_{}_3122.dat", "Smu_pT_{}_{}_3122_wSIP_BBPP.dat",
        "Smu_pT_{}_{}_3122_wSIP_LY.dat", "Smu_pT_{}_{}_3122_wMuIP_wSIP_LY.dat"
    ]

    pT_arr = []
    dN_list = []
    Sx_list = []
    Sy_list = []
    Sz_list = []
    for ifile in filelist:
        dN_list.append([])
        Sx_list.append([])
        Sy_list.append([])
        Sz_list.append([])

    nev = len(eventList_)
    for eventName in eventList_:
        event_group = hf_.get(eventName)
        for ifile, filename in enumerate(filelist):
            spin_data = nan_to_num(
                event_group.get(filename.format(vorticityType, rapType)))

            if ifile == 0:
                pT_arr = spin_data[:, 0]

            dN_list[ifile].append(spin_data[:, 1])
            Sx_list[ifile].append(spin_data[:, 7])
            Sy_list[ifile].append(spin_data[:, 8])
            Sz_list[ifile].append(spin_data[:, 9])

    Sx_avg = []
    Sx_err = []
    Sy_avg = []
    Sy_err = []
    Sz_avg = []
    Sz_err = []
    for ifile in range(len(filelist)):
        Sx_avg.append(
            mean(array(dN_list[ifile])*array(Sx_list[ifile]), axis=0)
            /mean(array(dN_list[ifile]), axis=0))
        Sx_err.append(
            std(array(dN_list[ifile])*array(Sx_list[ifile]), axis=0)
            /mean(array(dN_list[ifile]), axis=0)/sqrt(nev))
        Sy_avg.append(
            mean(array(dN_list[ifile])*array(Sy_list[ifile]), axis=0)
            /mean(array(dN_list[ifile]), axis=0))
        Sy_err.append(
            std(array(dN_list[ifile])*array(Sy_list[ifile]), axis=0)
            /mean(array(dN_list[ifile]), axis=0)/sqrt(nev))
        Sz_avg.append(
            mean(array(dN_list[ifile])*array(Sz_list[ifile]), axis=0)
            /mean(array(dN_list[ifile]), axis=0))
        Sz_err.append(
            std(array(dN_list[ifile])*array(Sz_list[ifile]), axis=0)
            /mean(array(dN_list[ifile]), axis=0)/sqrt(nev))

    f = open(
        path.join(outputFolder_,
                  "averaged_Smu_pT_{}.txt".format(vorticityType)), "w")
    f.write("# pT  S^x  S^x_err  S^y  S^y_err  S^z  S^z_err  "
            + "({0}  {0}+SIP(BBPP) {0}+SIP(LY)  {0}+SIP+MuBIP(LY))\n".format(
                vorticityType))
    for ipT in range(len(pT_arr)):
        f.write("%.5e  " % pT_arr[ipT])
        for icol in range(len(filelist)):
            f.write("%.5e  %.5e  " % (Sx_avg[icol][ipT], Sx_err[icol][ipT]))
            f.write("%.5e  %.5e  " % (Sy_avg[icol][ipT], Sy_err[icol][ipT]))
            f.write("%.5e  %.5e  " % (Sz_avg[icol][ipT], Sz_err[icol][ipT]))
        f.write("\n")
    f.close()


def analyze_Smu_phi(hf_, eventList_, outputFolder_, vnArr_, vnRef1_, vnRef2_,
                    iorder_, globalOutputFolder_, icen_):
    """
        This function computes the event-averaged S^mu(phi) with respect to
        the anisotropic flow angle vn.

        iorder_ = 0: phi_n = 0 corresponds to the reaction plane
    """
    filelist = [
        "Smu_phi_{}_{}_3122.dat", "Smu_phi_{}_{}_3122_wSIP_BBPP.dat",
        "Smu_phi_{}_{}_3122_wSIP_LY.dat", "Smu_phi_{}_{}_3122_wMuIP_wSIP_LY.dat"
    ]

    vnArr_ = array(vnArr_)
    vnRef1_ = array(vnRef1_)
    vnRef2_ = array(vnRef2_)

    phi_arr = []
    dN_list = []
    Sx_list = []
    Sy_list = []
    Sz_list = []
    fnSx_list = []
    fnSy_list = []
    fnSz_list = []
    for ifile in filelist:
        dN_list.append([])
        Sx_list.append([])
        Sy_list.append([])
        Sz_list.append([])
        fnSx_list.append([])
        fnSy_list.append([])
        fnSz_list.append([])

    nev = len(eventList_)
    for ievent, eventName in enumerate(eventList_):
        Qn = 0.
        psi_n = 0.
        if iorder_ > 0:
            Qn = vnArr_[ievent, iorder_ + 1]
            psi_n = (arctan2(imag(Qn), real(Qn))
                     + randrange(iorder_)*2.*pi)/iorder_

        event_group = hf_.get(eventName)
        for ifile, filename in enumerate(filelist):
            spin_data = nan_to_num(
                event_group.get(filename.format(vorticityType, rapType)))

            if ifile == 0:
                phi_arr = spin_data[:, 0]
            dphi = phi_arr[1] - phi_arr[0]

            fnSx_list[ifile].append(
                (sum(spin_data[:, 7]*cos(iorder_*phi_arr))*dphi/
                 (2*pi) + 1j*sum(spin_data[:, 7]*sin(iorder_*phi_arr))*dphi/
                 (2*pi)))
            fnSy_list[ifile].append(
                (sum(spin_data[:, 8]*cos(iorder_*phi_arr))*dphi/
                 (2*pi) + 1j*sum(spin_data[:, 8]*sin(iorder_*phi_arr))*dphi/
                 (2*pi)))
            fnSz_list[ifile].append(
                (sum(spin_data[:, 9]*cos(iorder_*phi_arr))*dphi/
                 (2*pi) + 1j*sum(spin_data[:, 9]*sin(iorder_*phi_arr))*dphi/
                 (2*pi)))

            phi_local = phi_arr + psi_n  # the "+" sign is correct
            # phi_arr = phi - Psi_n
            for ii in range(len(phi_local)):
                while phi_local[ii] < 0:
                    phi_local[ii] += 2.*pi
                while phi_local[ii] > 2.*pi:
                    phi_local[ii] -= 2.*pi

            dN = interp(phi_local, phi_arr, spin_data[:, 1])
            Sx = interp(phi_local, phi_arr, spin_data[:, 7])
            Sy = interp(phi_local, phi_arr, spin_data[:, 8])
            Sz = interp(phi_local, phi_arr, spin_data[:, 9])

            dN_list[ifile].append(dN)
            Sx_list[ifile].append(Sx)
            Sy_list[ifile].append(Sy)
            Sz_list[ifile].append(Sz)

    dN_array = real(vnArr_[:, 0])
    dN_avg = mean(dN_array)
    delta_N = dN_array - dN_avg
    varN = std(dN_array)**2.
    Qn_magsq = abs(vnArr_[:, iorder_ + 1])**2.
    delta_Qn = Qn_magsq - mean(Qn_magsq)
    hat_delta_Qn = delta_Qn - mean(delta_Qn*delta_N)/varN*delta_N

    QnRef1 = vnRef1_[:, iorder_ + 1]
    QnRef2 = vnRef2_[:, iorder_ + 1]
    Sx_avg = []
    Sx_err = []
    Sy_avg = []
    Sy_err = []
    Sz_avg = []
    Sz_err = []
    fnSx_avg = []
    fnSx_err = []
    fnSy_avg = []
    fnSy_err = []
    fnSz_avg = []
    fnSz_err = []
    rho_vnfnSz_avg = []
    rho_vnfnSz_err = []
    for ifile in range(len(filelist)):
        Sx_avg.append(
            mean(array(dN_list[ifile])*array(Sx_list[ifile]), axis=0)
            /mean(array(dN_list[ifile]), axis=0))
        Sx_err.append(
            std(array(dN_list[ifile])*array(Sx_list[ifile]), axis=0)
            /mean(array(dN_list[ifile]), axis=0)/sqrt(nev))
        Sy_avg.append(
            mean(array(dN_list[ifile])*array(Sy_list[ifile]), axis=0)
            /mean(array(dN_list[ifile]), axis=0))
        Sy_err.append(
            std(array(dN_list[ifile])*array(Sy_list[ifile]), axis=0)
            /mean(array(dN_list[ifile]), axis=0)/sqrt(nev))
        Sz_avg.append(
            mean(array(dN_list[ifile])*array(Sz_list[ifile]), axis=0)
            /mean(array(dN_list[ifile]), axis=0))
        Sz_err.append(
            std(array(dN_list[ifile])*array(Sz_list[ifile]), axis=0)
            /mean(array(dN_list[ifile]), axis=0)/sqrt(nev))

        # compute the error using the jackknife method
        fnSxArr = array(fnSx_list[ifile])
        fnSyArr = array(fnSy_list[ifile])
        fnSzArr = array(fnSz_list[ifile])
        fnSxReal = zeros(nev)
        fnSxImag = zeros(nev)
        fnSyReal = zeros(nev)
        fnSyImag = zeros(nev)
        fnSzReal = zeros(nev)
        fnSzImag = zeros(nev)
        rho_vnfnSz = zeros(nev)
        fnSzMagsq = abs(fnSzArr)**2.
        delta_fnSz = fnSzMagsq - mean(fnSzMagsq)
        hat_delta_fnSz = delta_fnSz - mean(delta_fnSz*delta_N)/varN*delta_N
        for iev in range(nev):
            array_idx = [True]*nev
            array_idx[iev] = False
            array_idx = array(array_idx)

            fnSxReal[iev] = (
                real(
                    mean(
                        fnSxArr[array_idx]*
                        (conj(QnRef1[array_idx]) + conj(QnRef2[array_idx]))/2.))
                /sqrt(real(mean(QnRef1[array_idx]*conj(QnRef2[array_idx])))))
            fnSxImag[iev] = (
                imag(
                    mean(
                        fnSxArr[array_idx]*
                        (conj(QnRef1[array_idx]) + conj(QnRef2[array_idx]))/2.))
                /sqrt(real(mean(QnRef1[array_idx]*conj(QnRef2[array_idx])))))
            fnSyReal[iev] = (
                real(
                    mean(
                        fnSyArr[array_idx]*
                        (conj(QnRef1[array_idx]) + conj(QnRef2[array_idx]))/2.))
                /sqrt(real(mean(QnRef1[array_idx]*conj(QnRef2[array_idx])))))
            fnSyImag[iev] = (
                imag(
                    mean(
                        fnSyArr[array_idx]*
                        (conj(QnRef1[array_idx]) + conj(QnRef2[array_idx]))/2.))
                /sqrt(real(mean(QnRef1[array_idx]*conj(QnRef2[array_idx])))))
            fnSzReal[iev] = (
                real(
                    mean(
                        fnSzArr[array_idx]*
                        (conj(QnRef1[array_idx]) + conj(QnRef2[array_idx]))/2.))
                /sqrt(real(mean(QnRef1[array_idx]*conj(QnRef2[array_idx])))))
            fnSzImag[iev] = (
                imag(
                    mean(
                        fnSzArr[array_idx]*
                        (conj(QnRef1[array_idx]) + conj(QnRef2[array_idx]))/2.))
                /sqrt(real(mean(QnRef1[array_idx]*conj(QnRef2[array_idx])))))

            rho_vnfnSz[iev] = (
                mean(hat_delta_fnSz[array_idx]*hat_delta_Qn[array_idx])/sqrt(
                    mean(hat_delta_fnSz[array_idx]**2.)
                    *mean(hat_delta_Qn[array_idx]**2.)))
        fnSx_avg.append(mean(fnSxReal) + 1j*mean(fnSxImag))
        fnSx_err.append(
            sqrt((nev - 1.)/nev*sum((fnSxReal - mean(fnSxReal))**2.))
            + 1j*sqrt((nev - 1.)/nev*sum((fnSxImag - mean(fnSxImag))**2.)))
        fnSy_avg.append(mean(fnSyReal) + 1j*mean(fnSyImag))
        fnSy_err.append(
            sqrt((nev - 1.)/nev*sum((fnSyReal - mean(fnSyReal))**2.))
            + 1j*sqrt((nev - 1.)/nev*sum((fnSyImag - mean(fnSyImag))**2.)))
        fnSz_avg.append(mean(fnSzReal) + 1j*mean(fnSzImag))
        fnSz_err.append(
            sqrt((nev - 1.)/nev*sum((fnSzReal - mean(fnSzReal))**2.))
            + 1j*sqrt((nev - 1.)/nev*sum((fnSzImag - mean(fnSzImag))**2.)))
        rho_vnfnSz_avg.append(mean(rho_vnfnSz))
        rho_vnfnSz_err.append(
            sqrt((nev - 1.)/nev*sum((rho_vnfnSz - mean(rho_vnfnSz))**2.)))

    fileLabel = "Psi{}".format(iorder_)
    if iorder_ == 0:
        fileLabel = "RP"
    f = open(
        path.join(outputFolder_,
                  "averaged_Smu_phi_{}_{}.txt".format(fileLabel,
                                                      vorticityType)), "w")
    f.write("# phi  S^x  S^x_err  S^y  S^y_err  S^z  S^z_err  "
            + "({0}  {0}+SIP(BBP)  {0}+SIP(LY)  {0}+SIP+MuBIP(LY))\n".format(
                vorticityType))
    for iphi in range(len(phi_arr)):
        f.write("%.5e  " % phi_arr[iphi])
        for icol in range(len(filelist)):
            f.write("%.5e  %.5e  " % (Sx_avg[icol][iphi], Sx_err[icol][iphi]))
            f.write("%.5e  %.5e  " % (Sy_avg[icol][iphi], Sy_err[icol][iphi]))
            f.write("%.5e  %.5e  " % (Sz_avg[icol][iphi], Sz_err[icol][iphi]))
        f.write("\n")
    f.close()

    if iorder_ > 0:
        dN_avg = mean(real(vnArr_[:, 0]))
        f = open(
            path.join(globalOutputFolder_,
                      "f{}_{}.txt".format(iorder_, vorticityType)), "a")
        if icen_ == 0:
            f.write("# cen  Nch  Re{fn(S^x)}  Re{fn(S^x)}_err  "
                    + "Im{fn(S^x)}  Im{fn(S^x)}_err  "
                    + "Re{fn(S^y)}  Re{fn(S^y)}_err  "
                    + "Im{fn(S^y)}  Im{fn(S^y)}_err  "
                    + "Re{fn(S^z)}  Re{fn(S^z)}_err  "
                    + "Im{fn(S^z)}  Im{fn(S^z)}_err  "
                    + "({0}  {0}+SIP(BBP) {0}+SIP(LY)  {0}+SIP+MuBIP(LY))\n".
                    format(vorticityType))
        f.write("{0}  {1:.5e}  ".format(
            (centrality_cut_list[icen_] + centrality_cut_list[icen_ + 1])/2.,
            dN_avg))
        for icol in range(len(filelist)):
            f.write("{0:.5e}  {1:.5e}  {2:.5e}  {3:.5e}  ".format(
                nan_to_num(real(fnSx_avg[icol])),
                nan_to_num(real(fnSx_err[icol])),
                nan_to_num(imag(fnSx_avg[icol])),
                nan_to_num(imag(fnSx_err[icol]))))
            f.write("{0:.5e}  {1:.5e}  {2:.5e}  {3:.5e}  ".format(
                nan_to_num(real(fnSy_avg[icol])),
                nan_to_num(real(fnSy_err[icol])),
                nan_to_num(imag(fnSy_avg[icol])),
                nan_to_num(imag(fnSy_err[icol]))))
            f.write("{0:.5e}  {1:.5e}  {2:.5e}  {3:.5e}  ".format(
                nan_to_num(real(fnSz_avg[icol])),
                nan_to_num(real(fnSz_err[icol])),
                nan_to_num(imag(fnSz_avg[icol])),
                nan_to_num(imag(fnSz_err[icol]))))
        f.write("\n")
        f.close()

        f = open(
            path.join(globalOutputFolder_,
                      "Rho_v{0}f{0}Pz_{1}.txt".format(iorder_, vorticityType)),
            "a")
        if icen_ == 0:
            f.write("# cen  Nch  rho(vn, fnPz)  rho(vn, fnPz)_err  "
                    + "({0}  {0}+SIP(BBPP) {0}+SIP(LY)  {0}+SIP+MuBIP(LY))\n".
                    format(vorticityType))
        f.write("{0}  {1:.5e}  ".format(
            (centrality_cut_list[icen_] + centrality_cut_list[icen_ + 1])/2.,
            dN_avg))
        for icol in range(len(filelist)):
            f.write("%.5e  %.5e  " %
                    (rho_vnfnSz_avg[icol], rho_vnfnSz_err[icol]))
        f.write("\n")
        f.close()


def analyze_Smu_y(hf_, eventList_, outputFolder_):
    """
        This function computes the event-averaged S^mu(y)
    """
    filelist = [
        "Smu_y_{}_{}_3122.dat", "Smu_y_{}_{}_3122_wSIP_BBPP.dat",
        "Smu_y_{}_{}_3122_wSIP_LY.dat", "Smu_y_{}_{}_3122_wMuIP_wSIP_LY.dat"
    ]

    y_arr = []
    dN_list = []
    Sx_list = []
    Sy_list = []
    Sz_list = []
    for ifile in filelist:
        dN_list.append([])
        Sx_list.append([])
        Sy_list.append([])
        Sz_list.append([])

    nev = len(eventList_)
    for eventName in eventList_:
        event_group = hf_.get(eventName)
        for ifile, filename in enumerate(filelist):
            spin_data = nan_to_num(
                event_group.get(filename.format(vorticityType, rapType)))
            if ifile == 0:
                y_arr = spin_data[:, 0]

            dN_list[ifile].append(spin_data[:, 1])
            Sx_list[ifile].append(spin_data[:, 7])
            Sy_list[ifile].append(spin_data[:, 8])
            Sz_list[ifile].append(spin_data[:, 9])

    Sx_avg = []
    Sx_err = []
    Sy_avg = []
    Sy_err = []
    Sz_avg = []
    Sz_err = []
    for ifile in range(len(filelist)):
        Sx_avg.append(
            mean(array(dN_list[ifile])*array(Sx_list[ifile]), axis=0)
            /mean(array(dN_list[ifile]), axis=0))
        Sx_err.append(
            std(array(dN_list[ifile])*array(Sx_list[ifile]), axis=0)
            /mean(array(dN_list[ifile]), axis=0)/sqrt(nev))
        Sy_avg.append(
            mean(array(dN_list[ifile])*array(Sy_list[ifile]), axis=0)
            /mean(array(dN_list[ifile]), axis=0))
        Sy_err.append(
            std(array(dN_list[ifile])*array(Sy_list[ifile]), axis=0)
            /mean(array(dN_list[ifile]), axis=0)/sqrt(nev))
        Sz_avg.append(
            mean(array(dN_list[ifile])*array(Sz_list[ifile]), axis=0)
            /mean(array(dN_list[ifile]), axis=0))
        Sz_err.append(
            std(array(dN_list[ifile])*array(Sz_list[ifile]), axis=0)
            /mean(array(dN_list[ifile]), axis=0)/sqrt(nev))

    f = open(
        path.join(outputFolder_, "averaged_Smu_y_{}.txt".format(vorticityType)),
        "w")
    f.write("# eta  S^x  S^x_err  S^y  S^y_err  S^z  S^z_err  "
            + "({0}  {0}+SIP(BBPP)  {0}+SIP(LY)  {0}+SIP+MuBIP(LY))\n".format(
                vorticityType))
    for iy in range(len(y_arr)):
        f.write("%.5e  " % y_arr[iy])
        for icol in range(len(filelist)):
            f.write("%.5e  %.5e  " % (Sx_avg[icol][iy], Sx_err[icol][iy]))
            f.write("%.5e  %.5e  " % (Sy_avg[icol][iy], Sy_err[icol][iy]))
            f.write("%.5e  %.5e  " % (Sz_avg[icol][iy], Sz_err[icol][iy]))
        f.write("\n")
    f.close()


def analyze_spin_vn_pTdiff(hf_, eventList_, outputFolder_, vnArr_, vnRef1_,
                           vnRef2_, iorder_, globalOutputFolder_, icen_):
    """
        This function computes the event-averaged Fourier coefficents of
        S^mu_n(pT) with respect to the anisotropic flow angle vn.

        iorder_ = 0: phi_n = 0 corresponds to the reaction plane
    """
    filelist = [
        "Smu_dpTdphi_{}_{}_3122.dat", "Smu_dpTdphi_{}_{}_3122_wSIP_BBPP.dat",
        "Smu_dpTdphi_{}_{}_3122_wSIP_LY.dat",
        "Smu_dpTdphi_{}_{}_3122_wMuIP_wSIP_LY.dat"
    ]

    vnArr_ = array(vnArr_)
    vnRef1_ = array(vnRef1_)
    vnRef2_ = array(vnRef2_)

    NPHI = 48
    NPT = 30
    phi_arr = []
    pT_arr = []
    dN_list = []
    Sx_list = []
    Sy_list = []
    Sz_list = []
    fnSx_list = []
    fnSy_list = []
    fnSz_list = []
    for ifile in filelist:
        dN_list.append([])
        Sx_list.append([])
        Sy_list.append([])
        Sz_list.append([])
        fnSx_list.append([])
        fnSy_list.append([])
        fnSz_list.append([])

    nev = len(eventList_)
    for ievent, eventName in enumerate(eventList_):
        event_group = hf_.get(eventName)
        for ifile, filename in enumerate(filelist):
            spin_data = nan_to_num(
                event_group.get(filename.format(vorticityType, rapType)))

            if ifile == 0:
                pT_arr = (spin_data[:, 0].reshape(NPT, NPHI))[:, 0]
                phi_arr = spin_data[0:NPHI, 1]
                pT_arr = pT_arr.reshape(NPT, 1)
                phi_arr = phi_arr.reshape(1, NPHI)

            dpT = pT_arr[1, 0] - pT_arr[0, 0]
            dphi = phi_arr[0, 1] - phi_arr[0, 0]

            dN_list[ifile].append(
                sum(spin_data[:, 2].reshape(NPT, NPHI), axis=1)*dphi)

            Sx_mat = spin_data[:, 8].reshape(NPT, NPHI)
            fnSx_list[ifile].append(
                (sum(Sx_mat*cos(iorder_*phi_arr), axis=1)*dphi/
                 (2*pi) + 1j*sum(Sx_mat*sin(iorder_*phi_arr), axis=1)*dphi/
                 (2*pi)))
            Sy_mat = spin_data[:, 9].reshape(NPT, NPHI)
            fnSy_list[ifile].append(
                (sum(Sy_mat*cos(iorder_*phi_arr), axis=1)*dphi/
                 (2*pi) + 1j*sum(Sy_mat*sin(iorder_*phi_arr), axis=1)*dphi/
                 (2*pi)))
            Sz_mat = spin_data[:, 10].reshape(NPT, NPHI)
            fnSz_list[ifile].append(
                (sum(Sz_mat*cos(iorder_*phi_arr), axis=1)*dphi/
                 (2*pi) + 1j*sum(Sz_mat*sin(iorder_*phi_arr), axis=1)*dphi/
                 (2*pi)))

    QnRef1 = (vnRef1_[:, iorder_ + 1]).reshape(nev, 1)
    QnRef2 = (vnRef2_[:, iorder_ + 1]).reshape(nev, 1)

    fnSx_avg = []
    fnSx_err = []
    fnSy_avg = []
    fnSy_err = []
    fnSz_avg = []
    fnSz_err = []
    for ifile in range(len(filelist)):
        # compute the error using the jackknife method
        fnSxArr = array(fnSx_list[ifile])  # dimension: [nev, NPT]
        fnSyArr = array(fnSy_list[ifile])
        fnSzArr = array(fnSz_list[ifile])

        fnSxReal = zeros([nev, NPT])
        fnSxImag = zeros([nev, NPT])
        fnSyReal = zeros([nev, NPT])
        fnSyImag = zeros([nev, NPT])
        fnSzReal = zeros([nev, NPT])
        fnSzImag = zeros([nev, NPT])

        for iev in range(nev):
            array_idx = [True]*nev
            array_idx[iev] = False
            array_idx = array(array_idx)

            fnSxReal[iev, :] = (
                real(
                    mean(fnSxArr[array_idx, :]*
                         (conj(QnRef1[array_idx]) + conj(QnRef2[array_idx]))/2.,
                         axis=0))
                /sqrt(real(mean(QnRef1[array_idx]*conj(QnRef2[array_idx])))))
            fnSxImag[iev, :] = (
                imag(
                    mean(fnSxArr[array_idx, :]*
                         (conj(QnRef1[array_idx]) + conj(QnRef2[array_idx]))/2.,
                         axis=0))
                /sqrt(real(mean(QnRef1[array_idx]*conj(QnRef2[array_idx])))))
            fnSyReal[iev, :] = (
                real(
                    mean(fnSyArr[array_idx, :]*
                         (conj(QnRef1[array_idx]) + conj(QnRef2[array_idx]))/2.,
                         axis=0))
                /sqrt(real(mean(QnRef1[array_idx]*conj(QnRef2[array_idx])))))
            fnSyImag[iev, :] = (
                imag(
                    mean(fnSyArr[array_idx, :]*
                         (conj(QnRef1[array_idx]) + conj(QnRef2[array_idx]))/2.,
                         axis=0))
                /sqrt(real(mean(QnRef1[array_idx]*conj(QnRef2[array_idx])))))
            fnSzReal[iev, :] = (
                real(
                    mean(fnSzArr[array_idx, :]*
                         (conj(QnRef1[array_idx]) + conj(QnRef2[array_idx]))/2.,
                         axis=0))
                /sqrt(real(mean(QnRef1[array_idx]*conj(QnRef2[array_idx])))))
            fnSzImag[iev, :] = (
                imag(
                    mean(fnSzArr[array_idx, :]*
                         (conj(QnRef1[array_idx]) + conj(QnRef2[array_idx]))/2.,
                         axis=0))
                /sqrt(real(mean(QnRef1[array_idx]*conj(QnRef2[array_idx])))))

        fnSx_avg.append(mean(fnSxReal, axis=0) + 1j*mean(fnSxImag, axis=0))
        fnSx_err.append(
            sqrt((nev - 1.)/nev
                 *sum((fnSxReal - mean(fnSxReal, axis=0))**2., axis=0))
            + 1j*sqrt((nev - 1.)/nev*sum(
                (fnSxImag - mean(fnSxImag, axis=0))**2., axis=0)))
        fnSy_avg.append(mean(fnSyReal, axis=0) + 1j*mean(fnSyImag, axis=0))
        fnSy_err.append(
            sqrt((nev - 1.)/nev
                 *sum((fnSyReal - mean(fnSyReal, axis=0))**2., axis=0))
            + 1j*sqrt((nev - 1.)/nev*sum(
                (fnSyImag - mean(fnSyImag, axis=0))**2., axis=0)))
        fnSz_avg.append(mean(fnSzReal, axis=0) + 1j*mean(fnSzImag, axis=0))
        fnSz_err.append(
            sqrt((nev - 1.)/nev
                 *sum((fnSzReal - mean(fnSzReal, axis=0))**2., axis=0))
            + 1j*sqrt((nev - 1.)/nev*sum(
                (fnSzImag - mean(fnSzImag, axis=0))**2., axis=0)))

    f = open(
        path.join(
            globalOutputFolder_, "f{}_pTdiff_C{}-{}_{}.txt".format(
                iorder_, int(centrality_cut_list[icen_]),
                int(centrality_cut_list[icen_ + 1]), vorticityType)), "w")

    f.write("# pT[GeV]  Nch  Re{fn(S^x)}  Re{fn(S^x)}_err  "
            + "Im{fn(S^x)}  Im{fn(S^x)}_err  "
            + "Re{fn(S^y)}  Re{fn(S^y)}_err  "
            + "Im{fn(S^y)}  Im{fn(S^y)}_err  "
            + "Re{fn(S^z)}  Re{fn(S^z)}_err  "
            + "Im{fn(S^z)}  Im{fn(S^z)}_err  "
            + "({0}  {0}+SIP(BBP) {0}+SIP(LY)  {0}+SIP+MuBIP(LY))\n".format(
                vorticityType))
    dN_avg = mean(dN_list[0], axis=0)
    for ipT in range(NPT):
        f.write("{0}  {1:.5e}  ".format(pT_arr[ipT, 0], dN_avg[ipT]))
        for icol in range(len(filelist)):
            f.write("{0:.5e}  {1:.5e}  {2:.5e}  {3:.5e}  ".format(
                nan_to_num(real(fnSx_avg[icol][ipT])),
                nan_to_num(real(fnSx_err[icol][ipT])),
                nan_to_num(imag(fnSx_avg[icol][ipT])),
                nan_to_num(imag(fnSx_err[icol][ipT]))))
            f.write("{0:.5e}  {1:.5e}  {2:.5e}  {3:.5e}  ".format(
                nan_to_num(real(fnSy_avg[icol][ipT])),
                nan_to_num(real(fnSy_err[icol][ipT])),
                nan_to_num(imag(fnSy_avg[icol][ipT])),
                nan_to_num(imag(fnSy_err[icol][ipT]))))
            f.write("{0:.5e}  {1:.5e}  {2:.5e}  {3:.5e}  ".format(
                nan_to_num(real(fnSz_avg[icol][ipT])),
                nan_to_num(real(fnSz_err[icol][ipT])),
                nan_to_num(imag(fnSz_avg[icol][ipT])),
                nan_to_num(imag(fnSz_err[icol][ipT]))))
        f.write("\n")
    f.close()


def analyze_Rspin_y(hf_, eventList_, outputFolder_: str, pTmin_: float,
                    pTmax_: float) -> None:
    """
        This function computes the event averaged R_spin(pT, y)
    """
    filelist = [
        "Rspin_pTy_{}_{}_3122.dat", "Rspin_pTy_{}_{}_3122_wSIP_BBPP.dat",
        "Rspin_pTy_{}_{}_3122_wSIP_LY.dat",
        "Rspin_pTy_{}_{}_3122_wMuIP_wSIP_LY.dat"
    ]

    dNList = []
    RspinList = []
    for ifile in filelist:
        RspinList.append([])
        dNList.append([])

    nev = len(eventList_)
    for eventName in eventList_:
        event_group = hf_.get(eventName)
        for ifile, filename in enumerate(filelist):
            Rspin_data = nan_to_num(
                event_group.get(filename.format(vorticityType, rapType)))

            if ifile == 0:
                yarr = Rspin_data[:, 0].reshape(-1, 30)[:, 0]
                pTarr = Rspin_data[:, 1].reshape(-1, 30)[0, :]
            dNdpTdy = Rspin_data[:, 2].reshape(-1, 30)
            Rspin = Rspin_data[:, 3].reshape(-1, 30)

            idx = (pTarr > pTmin_) & (pTarr < pTmax_)
            RspinList[ifile].append(sum(dNdpTdy[:, idx]*Rspin[:, idx], axis=1))
            dNList[ifile].append(sum(dNdpTdy[:, idx], axis=1))

    R_spin_avg = []
    R_spin_err = []
    for ifile in range(len(filelist)):
        R_spin_avg.append(
            mean(array(RspinList[ifile]), axis=0)
            /mean(array(dNList[ifile]), axis=0))
        R_spin_err.append(
            std(array(RspinList[ifile]), axis=0)
            /mean(array(dNList[ifile]), axis=0)/sqrt(nev))
    fileName = (
        f"averaged_Rspin_{rapType}_pT_{pTmin_}_{pTmax_}_{vorticityType}.txt")
    f = open(path.join(outputFolder_, fileName), "w")
    f.write("# eta  Rspin  Rspin_err "
            + "({0}  {0}+SIP(BBPP)  {0}+SIP(LY)  {0}+SIP+MuBIP(LY))\n".format(
                vorticityType))
    for iy in range(len(yarr)):
        f.write("%.5e  " % yarr[iy])
        for icol in range(len(filelist)):
            f.write("%.5e  %.5e  " %
                    (R_spin_avg[icol][iy], R_spin_err[icol][iy]))
        f.write("\n")
    f.close()


def analyze_Rspin_pT(hf_, eventList_, outputFolder_: str, rapMin_: float,
                     rapMax_: float) -> None:
    """
        This function computes the event averaged R_spin(pT, y)
    """
    filelist = [
        "Rspin_pTy_{}_{}_3122.dat", "Rspin_pTy_{}_{}_3122_wSIP_BBPP.dat",
        "Rspin_pTy_{}_{}_3122_wSIP_LY.dat",
        "Rspin_pTy_{}_{}_3122_wMuIP_wSIP_LY.dat"
    ]

    dNList = []
    RspinList = []
    for ifile in filelist:
        RspinList.append([])
        dNList.append([])

    nev = len(eventList_)
    for eventName in eventList_:
        event_group = hf_.get(eventName)
        for ifile, filename in enumerate(filelist):
            Rspin_data = nan_to_num(
                event_group.get(filename.format(vorticityType, rapType)))

            if ifile == 0:
                yarr = Rspin_data[:, 0].reshape(-1, 30)[:, 0]
                pTarr = Rspin_data[:, 1].reshape(-1, 30)[0, :]
            dNdpTdy = Rspin_data[:, 2].reshape(-1, 30)
            Rspin = Rspin_data[:, 3].reshape(-1, 30)

            idx = (yarr > rapMin_) & (yarr < rapMax_)
            RspinList[ifile].append(sum(dNdpTdy[idx, :]*Rspin[idx, :], axis=0))
            dNList[ifile].append(sum(dNdpTdy[idx, :], axis=0))

    R_spin_avg = []
    R_spin_err = []
    for ifile in range(len(filelist)):
        R_spin_avg.append(
            mean(array(RspinList[ifile]), axis=0)
            /mean(array(dNList[ifile]), axis=0))
        R_spin_err.append(
            std(array(RspinList[ifile]), axis=0)
            /mean(array(dNList[ifile]), axis=0)/sqrt(nev))

    fileName = (
        f"averaged_Rspin_pT_{rapType}_{rapMin_}_{rapMax_}_{vorticityType}.txt")
    f = open(path.join(outputFolder_, fileName), "w")
    f.write("# pT (GeV)  Rspin  Rspin_err "
            + "({0}  {0}+SIP(BBPP)  {0}+SIP(LY)  {0}+SIP+MuBIP(LY))\n".format(
                vorticityType))
    for ipT in range(len(pTarr)):
        f.write("%.5e  " % pTarr[ipT])
        for icol in range(len(filelist)):
            f.write("%.5e  %.5e  " %
                    (R_spin_avg[icol][ipT], R_spin_err[icol][ipT]))
        f.write("\n")
    f.close()


###############################
# Analysis Starts from here ...
###############################
hf = h5py.File(data_path, "r")
event_list = list(hf.keys())
print("total number of events: {}".format(len(event_list)))

dNdyDict = {}
for ifolder, event_name in enumerate(event_list):
    file_name = "particle_9999_vndata_eta_-0.5_0.5.dat"
    event_group = hf.get(event_name)
    eventStatus = check_an_event_is_good(event_group)
    if eventStatus:
        temp_data = event_group.get(file_name)
        temp_data = nan_to_num(temp_data)
        dNdyDict[event_name] = temp_data[0, 1]
dNdyList = -sort(-array(list(dNdyDict.values())))
print("Number of good events: {}".format(len(dNdyList)))

for icen in range(len(centrality_cut_list) - 1):
    if centrality_cut_list[icen + 1] < centrality_cut_list[icen]:
        continue
    avg_folder = path.join(
        avg_folder_header,
        "{0:02.0f}-{1:02.0f}".format(centrality_cut_list[icen],
                                     centrality_cut_list[icen + 1]))
    mkdir(avg_folder)

    selected_events_list = []
    dN_dy_cut_high = (dNdyList[int(
        len(dNdyList)*centrality_cut_list[icen]/100.)])
    dN_dy_cut_low = dNdyList[min(
        len(dNdyList) - 1,
        int(len(dNdyList)*centrality_cut_list[icen + 1]/100.))]
    if len(dNcutList) == len(centrality_cut_list):
        dN_dy_cut_high = dNcutList[icen]
        dN_dy_cut_low = dNcutList[icen + 1]

    for event_name in dNdyDict.keys():
        if (dNdyDict[event_name] > dN_dy_cut_low
                and dNdyDict[event_name] <= dN_dy_cut_high):
            selected_events_list.append(event_name)

    nev = len(selected_events_list)
    print("analysis {}%-{}% nev = {}...".format(centrality_cut_list[icen],
                                                centrality_cut_list[icen + 1],
                                                nev))
    print("dNdy: {0:.2f} - {1:.2f}".format(dN_dy_cut_low, dN_dy_cut_high))
    if nev == 0:
        print("Skip ...")
        continue

    vnFileName = 'particle_9999_vndata_diff_eta_-0.5_0.5.dat'
    vnRefFileName1 = 'particle_9999_vndata_diff_eta_0.5_1.dat'
    vnRefFileName2 = 'particle_9999_vndata_diff_eta_-1_-0.5.dat'

    pT_array = []
    dN_array = []
    vn_star_array = []
    vn_star_array_ref1 = []
    vn_star_array_ref2 = []
    vn_alice_array = []
    vn_alice_array_ref1 = []
    vn_alice_array_ref2 = []
    QnpT_diff_star = []
    Qnref_star = []
    QnpT_diff_alice = []
    Qnref_alice = []
    for ifolder, event_name in enumerate(selected_events_list):
        event_group = hf.get(event_name)
        temp_data = nan_to_num(event_group.get(vnFileName))
        temp_data_ref1 = nan_to_num(event_group.get(vnRefFileName1))
        temp_data_ref2 = nan_to_num(event_group.get(vnRefFileName2))

        dN_event = temp_data[:, 1]  # dN/(2pi dy pT dpT)
        pT_event = temp_data[:, 0]

        # record particle spectra
        pT_array.append(pT_event)
        dN_array.append(dN_event)

        # pT-integrated vn
        # vn with STAR pT cut
        temp_vn_array = calculate_meanpT_inte_vn(0.2, 2.0, temp_data)
        vn_star_array.append(temp_vn_array)
        temp_vn_array = calculate_meanpT_inte_vn(0.2, 2.0, temp_data_ref1)
        vn_star_array_ref1.append(temp_vn_array)
        temp_vn_array = calculate_meanpT_inte_vn(0.2, 2.0, temp_data_ref2)
        vn_star_array_ref2.append(temp_vn_array)

        # vn with ALICE pT cut
        temp_vn_array = calculate_meanpT_inte_vn(0.2, 3.0, temp_data)
        vn_alice_array.append(temp_vn_array)
        temp_vn_array = calculate_meanpT_inte_vn(0.2, 3.0, temp_data_ref1)
        vn_alice_array_ref1.append(temp_vn_array)
        temp_vn_array = calculate_meanpT_inte_vn(0.2, 3.0, temp_data_ref2)
        vn_alice_array_ref2.append(temp_vn_array)

        # pT-differential vn using scalar-product method
        # vn{SP}(pT) with STAR pT cut
        temp_arr = calculate_diff_vn_single_event(0.2, 2.0, temp_data,
                                                  temp_data_ref2)
        QnpT_diff_star.append(temp_arr[0])
        Qnref_star.append(temp_arr[1])

        # vn{SP}(pT) with ALICE pT cut
        temp_arr = calculate_diff_vn_single_event(0.2, 3.0, temp_data,
                                                  temp_data_ref2)
        QnpT_diff_alice.append(temp_arr[0])
        Qnref_alice.append(temp_arr[1])

    # now we perform event average
    dN_array = array(dN_array)
    pT_array = array(pT_array)

    n_pT = len(pT_array[0, :])
    pT_spectra = zeros([n_pT])
    for ipT in range(len(pT_array[0, :])):
        dN_temp = sum(dN_array[:, ipT]*pT_array[:, ipT])
        if dN_temp > 0:
            pT_spectra[ipT] = (sum(pT_array[:, ipT]**2.*dN_array[:, ipT])
                               /dN_temp)
        else:
            pT_spectra[ipT] = mean(pT_array[:, ipT])
    # dN/(2pi dy pT dpT)
    dN_spectra = mean(pT_array*dN_array, 0)/pT_spectra
    dN_spectra_err = std(pT_array*dN_array, 0)/pT_spectra/sqrt(nev)

    # calcualte dN/dy and <pT>
    vn_star_array = array(vn_star_array)
    vn_alice_array = array(vn_alice_array)
    dNdy_avg_star = real(mean(vn_star_array[:, 0]))
    dNdy_avg_star_err = real(std(vn_star_array[:, 0]))/sqrt(nev)
    meanpT_star = real(mean(vn_star_array[:, 1]))
    meanpT_star_err = real(std(vn_star_array[:, 1]))/sqrt(nev)
    dNdy_avg_alice = real(mean(vn_alice_array[:, 0]))
    dNdy_avg_alice_err = real(std(vn_alice_array[:, 0]))/sqrt(nev)
    meanpT_alice = real(mean(vn_alice_array[:, 1]))
    meanpT_alice_err = real(std(vn_alice_array[:, 1]))/sqrt(nev)

    # calcualte vn{2}
    vn_star_2_gap, vn_star_2_gap_err = calcualte_vn_2_with_gap(
        vn_star_array_ref1, vn_star_array_ref2)
    vn_alice_2_gap, vn_alice_2_gap_err = calcualte_vn_2_with_gap(
        vn_alice_array_ref1, vn_alice_array_ref2)

    # calcualte vn{SP}(pT)
    vn_diff_SP_star = calculate_vn_diff_SP(QnpT_diff_star, Qnref_star)
    vn_diff_SP_alice = calculate_vn_diff_SP(QnpT_diff_alice, Qnref_alice)

    # then particle rapidity distribution
    rapFileName = 'particle_9999_dNdeta_pT_0.2_3.dat'

    eta_array = []
    dN_array = []
    vn_array = []
    for ifolder, event_name in enumerate(selected_events_list):
        event_group = hf.get(event_name)
        temp_data = event_group.get(rapFileName)
        temp_data = nan_to_num(temp_data)

        eta_array.append(temp_data[:, 0])
        dN_array.append(temp_data[:, 1])
        temp_vn_array = []
        for iorder in range(1, n_order):
            vn_real = temp_data[:, 2*iorder + 1]
            vn_imag = temp_data[:, 2*iorder + 2]
            vn = vn_real + 1j*vn_imag
            temp_vn_array.append(vn)
        vn_array.append(temp_vn_array)

    eta_array = array(eta_array)
    dN_array = array(dN_array)
    vn_array = array(vn_array)

    eta_point = mean(eta_array, 0)
    dNdeta = mean(dN_array, 0)
    dNdeta_err = std(dN_array, 0)/sqrt(nev)
    vn_SP_eta, vn_SP_eta_err = calculate_vn_eta(eta_point, dN_array, vn_array,
                                                -5.1, -2.8)
    vn_SP_eta_mid, vn_SP_eta_mid_err = calculate_vn_eta(eta_point, dN_array,
                                                        vn_array, -0.5, 0.5)

    # analysis spin observables
    analyze_Smu(hf, selected_events_list, 0.5, 3.0, avg_folder_header, icen,
                vn_alice_array)
    analyze_Smu_pT(hf, selected_events_list, avg_folder)
    analyze_Smu_y(hf, selected_events_list, avg_folder)
    for iorder in range(n_order):
        analyze_Smu_phi(hf, selected_events_list, avg_folder, vn_alice_array,
                        vn_alice_array_ref1, vn_alice_array_ref2, iorder,
                        avg_folder_header, icen)
        analyze_spin_vn_pTdiff(hf, selected_events_list, avg_folder,
                               vn_alice_array, vn_alice_array_ref1,
                               vn_alice_array_ref2, iorder, avg_folder_header,
                               icen)
    analyze_Rspin_y(hf, selected_events_list, avg_folder, 0.5, 3.0)
    analyze_Rspin_pT(hf, selected_events_list, avg_folder, -0.5, 0.5)
    analyze_Rspin_pT(hf, selected_events_list, avg_folder, -2.5, -1.5)
    analyze_Rspin_pT(hf, selected_events_list, avg_folder, 1.5, 2.5)

    ######################################################################
    # finally, output all the results
    ######################################################################

    output_filename = "charged_hadrons_integrated_observables.dat"
    f = open(path.join(avg_folder, output_filename), 'w')
    f.write("dN/dy(STAR)= %.5e +/- %.5e\n" % (dNdy_avg_star, dNdy_avg_star_err))
    f.write("dN/dy(ALICE)= %.5e +/- %.5e\n" %
            (dNdy_avg_alice, dNdy_avg_alice_err))
    f.write("<pT>(STAR)= %.5e +/- %.5e\n" % (meanpT_star, meanpT_star_err))
    f.write("<pT>(ALICE)= %.5e +/- %.5e\n" % (meanpT_alice, meanpT_alice_err))
    for iorder in range(1, n_order):
        f.write(
            "v_%d{2}(STAR)= %.5e +/- %.5e\n" %
            (iorder, vn_star_2_gap[iorder - 1], vn_star_2_gap_err[iorder - 1]))
        f.write("v_%d{2}(ALICE)= %.5e +/- %.5e\n" %
                (iorder, vn_alice_2_gap[iorder - 1],
                 vn_alice_2_gap_err[iorder - 1]))
    f.close()

    output_filename = "charged_hadrons_differential_observables_STAR.dat"
    f = open(path.join(avg_folder, output_filename), 'w')
    f.write("#pT  dN/(2pi dy pT dpT)  dN/(2pi dy pT dpT)_err  "
            + "vn{SP}  vn{SP}_err\n")
    for ipT in range(len(pT_spectra)):
        f.write("%.10e  %.10e  %.10e  " %
                (pT_spectra[ipT], dN_spectra[ipT], dN_spectra_err[ipT]))
        for iorder in range(1, n_order):
            f.write("%.10e  %.10e  " % (vn_diff_SP_star[2*iorder - 2][ipT],
                                        vn_diff_SP_star[2*iorder - 1][ipT]))
        f.write("\n")
    f.close()

    output_filename = "charged_hadrons_differential_observables_ALICE.dat"
    f = open(path.join(avg_folder, output_filename), 'w')
    f.write("#pT  dN/(2pi dy pT dpT)  dN/(2pi dy pT dpT)_err  "
            + "vn{SP}  vn{SP}_err\n")
    for ipT in range(len(pT_spectra)):
        f.write("%.10e  %.10e  %.10e  " %
                (pT_spectra[ipT], dN_spectra[ipT], dN_spectra_err[ipT]))
        for iorder in range(1, n_order):
            f.write("%.10e  %.10e  " % (vn_diff_SP_alice[2*iorder - 2][ipT],
                                        vn_diff_SP_alice[2*iorder - 1][ipT]))
        f.write("\n")
    f.close()

    output_filename = "charged_hadrons_rapidity_distribution.dat"
    f = open(path.join(avg_folder, output_filename), 'w')
    f.write("#eta  dN/deta  dN/deta_err  vn{2}(eta)  vn{2}(eta)_err"
            + "  Re{vn}(eta) Re{vn}(eta)_err\n")
    for ieta in range(len(eta_point)):
        f.write("%.10e  %.10e  %.10e  " %
                (eta_point[ieta], dNdeta[ieta], dNdeta_err[ieta]))
        for iorder in range(1, n_order):
            f.write(
                "%.10e  %.10e  %.10e  %.10e  " %
                (vn_SP_eta[iorder - 1, ieta], vn_SP_eta_err[iorder - 1, ieta],
                 vn_SP_eta_mid[iorder - 1, ieta], vn_SP_eta_mid_err[iorder - 1,
                                                                    ieta]))
        f.write("\n")
    f.close()

print("Analysis is done.")
