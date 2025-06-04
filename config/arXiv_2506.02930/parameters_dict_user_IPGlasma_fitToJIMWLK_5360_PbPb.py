#!/usr/bin/env python3
"""
    This script contains all the user modified parameters in
    the iEBE-MUSIC package.
"""

# control parameters
control_dict = {
    'initial_state_type': "IPGlasma",  # 3DMCGlauber, IPGlasma
    'walltime': "120:00:00",           # walltime to run
    'save_ipglasma_results': False,
    'save_hydro_surfaces': False,   # flag to save hydro surfaces
    'save_UrQMD_files': False,      # flag to save UrQMD files
}


# IPGlasma
ipglasma_dict = {
    'type': "self",
    # all parameters below are for (type == self)
    'nucleonPositionsFromFile': 0,
    'bmin': 0.,
    'bmax': 25.,
    'size': 800,              # number of grid points of IP-Glasma computation
    'Projectile': "Pb",
    'Target': "Pb",
    'roots': 5360.,
    'SigmaNN': 68.3,
    'setWSDeformParams': 1,
    'R_WS': 6.647,     # 2108.09578
    'a_WS': 0.537,
    'beta2': 0.006,
    'beta3': 0.,
    'beta4': 0.,
    'gamma': 0.,
    'SubNucleonParamType': 4,    # 0: do not use posterior parameter sets
                                 # 1: use subnucleon parameters from variant Nq posterior distribution
                                 # 2: use subnucleon parameters from fixed Nq = 3 posterior distribution
                                 # 3: use subnucleon parameters from variant Nq posterior distribution (fit to LHC x)
                                 # 4: use subnucleon parameters from fixed Nq = 3 posterior distribution (fit to LHC x)
    'SubNucleonParamSet': 0,     # -1: choose a random set from the posterior distribution
                                 # 0: choose the MAP parameter set
                                 # positive intergers: choose a fixed set of parameter for sub-nucleonic structure
    'useConstituentQuarkProton': 3,
    'maxtime': 0.4,
    'LOutput': 34,
    'sizeOutput': 340,
}


# MUSIC
music_dict = {
    'beastMode': 2,
    'Initial_profile': 9,   # type of initial condition 
                            # 9: IPGlasma (full Tmunu),
                            #   -- 91: e and u^\mu,
                            #   -- 92: e only,
                            #   -- 93: e, u^\mu, and pi^\munu
    's_factor': 0.085,      # normalization factor read in initial data file
    'Initial_time_tau_0': 0.4,  # starting time of the hydrodynamic evolution (fm/c)
    'Delta_Tau': 0.005,         # time step to use in the evolution [fm/c]
    'boost_invariant':  1,      # whether the simulation is boost-invariant
    'EOS_to_use': 9,            # type of the equation of state
                                # 9: hotQCD EOS with UrQMD
    # transport coefficients
    'quest_revert_strength': 1.0,          # the strength of the viscous regulation
    'Viscosity_Flag_Yes_1_No_0': 1,        # turn on viscosity in the evolution
    'Include_Shear_Visc_Yes_1_No_0': 1,    # include shear viscous effect
    'Shear_to_S_ratio': 0.120,             # value of \eta/s
    'T_dependent_Shear_to_S_ratio': 3,     # flag to use temperature dep. \eta/s(T)
    'shear_viscosity_3_eta_over_s_T_kink_in_GeV': 0.18,
    'shear_viscosity_3_eta_over_s_low_T_slope_in_GeV': -4,
    'shear_viscosity_3_eta_over_s_high_T_slope_in_GeV': 0,
    'shear_viscosity_3_eta_over_s_at_kink': 0.12,
    'Include_Bulk_Visc_Yes_1_No_0': 1,     # include bulk viscous effect
    'T_dependent_zeta_over_s': 10,         # parameterization of \zeta/s(T)
    'bulk_viscosity_10_max': 0.12,         # the peak value of \zeta/s(T)
    'bulk_viscosity_10_T_peak': 0.180,     # the peak temperature for \zeta/s(T)
    'bulk_viscosity_10_width_high': 0.12,
    'bulk_viscosity_10_width_low': 0.025,
    'Include_second_order_terms': 1,       # include second order non-linear coupling terms
    'Include_vorticity_terms': 0,          # include vorticity coupling terms

    # parameters for freeze out and Cooper-Frye
    'N_freeze_out': 1,
    'eps_switch': 0.18,
}


# iSS
iss_dict = {
    'hydro_mode': 1,    # mode for reading in freeze out information 
    'include_deltaf_shear': 1,      # include delta f contribution from shear
    'include_deltaf_bulk': 1,       # include delta f contribution from bulk
    'bulk_deltaf_kind': 20,         # 21: relaxation time approximation (both shear and bulk)
    'sample_upto_desired_particle_number': 1,  # 1: flag to run sampling until desired
                                               # particle numbers is reached
    'number_of_particles_needed': 100000,      # number of hadrons to sample
    'local_charge_conservation': 0,  # flag to impose local charge conservation
    'global_momentum_conservation': 0,  # flag to impose GMC
}


# hadronic afterburner toolkit
hadronic_afterburner_toolkit_dict = {
    'analyze_flow': 4,                 # 0/1: flag to perform flow analysis
    'order_max': 16,
    'event_buffer_size': 100000,       # the number of events read in at once
    'compute_correlation': 0,       # flag to compute correlation function
    'flag_charge_dependence': 0,    # flag to compute charge dependence correlation
    'compute_corr_rap_dep': 0,      # flag to compute the rapidity dependent multi-particle correlation
    'resonance_weak_feed_down_flag': 0,  # include weak feed down contribution
    'npT': 20,          # number of pT points for pT-differential spectra and vn
    'pT_min': 0.0,      # the minimum value of transverse momentum (GeV)
    'pT_max': 3.8,      # the maximum value of transverse momentum (GeV)
    'n_rap': 71,                  # numpber of points in rapidity distr.
    'rapidity_dis_min': -7.0,     # minimum value of particle rapidity distribution
    'rapidity_dis_max': 7.0,      # maximum value of particle rapidity distribution
}
