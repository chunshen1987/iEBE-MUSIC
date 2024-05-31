#!/usr/bin/env python3
"""
    This script contains all the user modified parameters in
    the iEBE-MUSIC package.
"""


# control parameters
control_dict = {
    'initial_state_type': "3DMCGlauber_dynamical",
    'walltime': "120:00:00",        # walltime to run
    'save_hydro_surfaces': False,   # flag to save hydro surfaces
    'save_UrQMD_files': False,      # flag to save UrQMD files
}


# 3DMCGlauber model
mcglauber_dict = {
    'database_name': "self",     # self: generate initial condition on the fly
    'Projectile':  "Ru",         # projectile nucleus name
    'Target'    :  "Ru",         # target nucleus name
    'resetProjWS':  1,
    'resetTargWS':  1,
    'ProjWS_R':  5.09,
    'ProjWS_a':  0.46,
    'ProjWS_beta2':  0.16,
    'ProjWS_beta3':  0.0,
    'ProjWS_beta4':  0.,
    'ProjWS_gamma':  0.0,
    'ProjWS_da':  0.01,
    'ProjWS_dR':  0.015,
    'd_min': 0.9,
    'useQuarks': 1,
    'Q2': 1.,
    'roots' : 200,               # collision energy (GeV)
    'seed' : -1,                 # random seed (-1: system)
    'baryon_junctions': 1,       # 0: baryon number assumed to be at string end
                                 # 1: baryon number transported assuming baryon
                                 # junctions (at smaller x)
                                 # see arXiv:nucl-th/9602027
    'electric_junctions':  1,
    'integer_electric_charge': 1,
    'electricChargeinStringProb': 1.0,
    'lambdaB': 0.2,              # parameter the controls the strength of
    'lambdaQ': 0.0,              # parameter the controls the strength of
                                 # the baryon junction stopping
    'lambdaBs': 1.0,             # fraction of single-to-double string junction stopping
    'lambdaQs': 1.0,             # fraction of single-to-double string junction stopping
    'baryonInStringProb': 1.0,

    'BG': 16.,
    'shadowing_factor':  0.60,   # a shadowning factor for producing strings from multiple scatterings
    'rapidity_loss_method': 4,
    'remnant_energy_loss_fraction': 0.5,         # nucleon remnants energy loss fraction (fraction of string's y_loss) [0, 1]
    'ylossParam4At2': 1.70,
    'ylossParam4At4': 2.00,
    'ylossParam4At6': 2.20,
    'ylossParam4var': 0.5,
    'evolve_QCD_string_mode': 4,        # string evolution mode
                                        # 1: deceleration with fixed rapidity loss (m/sigma = 1 fm, dtau = 0.5 fm)
                                        # 2: deceleration with LEXUS sampled rapidit loss (both dtau and sigma fluctuate)
                                        # 3: deceleration with LEXUS sampled rapidit loss (m/sigma = 1 fm, dtau fluctuates)
                                        # 4: deceleration with LEXUS sampled rapidit loss (dtau = 0.5 fm, m/sigma fluctuates)
}


# MUSIC
music_dict = {
    'Initial_profile': 131,   # type of initial condition 
                              # 13: dynamical initialization (3dMCGlauber_dynamical)
                              #   -- 131: 3dMCGlauber with zero nucleus thickness
    'string_source_sigma_eta': 0.5,  # the smearning size of the hotspot in eta
    'string_source_sigma_x': 0.4,    # the smearning size of the hotspot in eta
    'stringTransverseShiftFrac': 1.0,  # control the shift of transverse coord as a function of eta for string
    'stringPreEqFlowFactor': 0.50,     # pre-Eq. flow factor

    's_factor': 1.000,      # normalization factor read in initial data file
    'Delta_Tau': 0.020,     # time step to use in the evolution [fm/c]
    'boost_invariant':  0,  # whether the simulation is boost-invariant
    'Eta_grid_size': 16.0,  # spatial rapidity range 
                            # [-Eta_grid_size/2, Eta_grid_size/2 - delta_eta]
    'Grid_size_in_eta': 64,     # number of the grid points in spatial rapidity direction
    'X_grid_size_in_fm': 30,    # spatial range along x direction in the transverse plane 
                                # [-X_grid_size_in_fm/2, X_grid_size_in_fm/2]
    'Y_grid_size_in_fm': 30,    # spatial range along x direction in the transverse plane 
                                # [-X_grid_size_in_fm/2, X_grid_size_in_fm/2]
    'Grid_size_in_x': 96,      # number of the grid points in x direction
    'Grid_size_in_y': 96,      # number of the grid points in y direction
    'EOS_to_use': 20,           # type of the equation of state
                                # 14: neos_BQS lattice EoS at finite mu_B
                                # 17: BEST lattice EoS at finite mu_B
    # transport coefficients
    'quest_revert_strength': 1.0,
    'Viscosity_Flag_Yes_1_No_0': 1,        # turn on viscosity in the evolution
    'Include_Shear_Visc_Yes_1_No_0': 1,    # include shear viscous effect
    'T_dependent_Shear_to_S_ratio': 3,     # flag to use temperature dep. \eta/s(T)
    'shear_viscosity_3_eta_over_s_T_kink_in_GeV': 0.17,
    'shear_viscosity_3_eta_over_s_low_T_slope_in_GeV': -0.0,
    'shear_viscosity_3_eta_over_s_high_T_slope_in_GeV': 0.,
    'shear_viscosity_3_eta_over_s_at_kink': 0.08,
    'muB_dependent_Shear_to_S_ratio': 10,
    'shear_muBDep_alpha': 0.7,
    'shear_muBDep_slope': 2.0,
    'shear_muBDep_scale': 0.6,             # GeV
    'Include_Bulk_Visc_Yes_1_No_0': 1,     # include bulk viscous effect
    'T_dependent_zeta_over_s': 10,         # parameterization of \zeta/s(T)
    'bulk_viscosity_10_max': 0.10,         # the peak value of \zeta/s(T)
    'bulk_viscosity_10_T_peak': 0.17,      # the peak temperature for \zeta/s(T)
    'bulk_viscosity_10_T_peak_muBcurv': -0.15,  # Tpeak = Tpeak_0 + curv*muB^2
    'bulk_viscosity_10_width_high': 0.08,
    'bulk_viscosity_10_width_low': 0.010,
    'Include_second_order_terms': 1,       # include second order non-linear coupling terms
    'Include_vorticity_terms': 0,          # include vorticity coupling terms
    'Include_Rhob_Yes_1_No_0': 1,
    'Include_QS_Yes_1_No_0': 1,            # use nq and ns in hydro evolution.
    'use_BQ_ratios': 0,                    # use nq = 0.4 nb in hydro source strings.
    'turn_on_baryon_diffusion': 0,
    'kappa_coefficient': 0.4,

    # parameters for freeze out and Cooper-Frye
    'freeze_out_tau_start_max': 2,      # the maximum freeze-out starting time [fm/c]
    'N_freeze_out': 1,
    'eps_switch': 0.35,                 # GeV/fm^3
    'beastMode': 2,                     # Release the Kraken
}


# iSS
iss_dict = {
    'hydro_mode': 2,                # mode for reading in freeze out information 
    'MC_sampling': 4,
    'include_deltaf_shear': 1,      # include delta f contribution from shear
    'include_deltaf_bulk': 1,       # include delta f contribution from bulk
    'bulk_deltaf_kind': 20,         # 20: 22-momentum approximation, 21: CE relaxation time approximation
    'include_deltaf_diffusion': 0,  # include delta f contribution from diffusion
    'sample_upto_desired_particle_number': 1,  # 1: flag to run sampling until desired
                                               # particle numbers is reached
    'number_of_particles_needed': 100000,      # number of hadrons to sample
    'local_charge_conservation': 0,     # flag to impose local charge conservation
    'global_momentum_conservation': 0,  # flag to impose GMC
}


# hadronic afterburner toolkit
hadronic_afterburner_toolkit_dict = {
    'analyze_HBT': 0,                    # 0/1: flag to perform HBT analysis
    'event_buffer_size': 500000,         # the number of events read in at once
    'rapidity_shift':  0.0,              # 0.5*log(Z_P*A_T/(A_P*Z_T))
    'compute_correlation': 0,            # flag to compute correlation function
    'flag_charge_dependence': 0,         # flag to compute charge dependence correlation
    'compute_corr_rap_dep': 0,           # flag to compute the rapidity dependent multi-particle correlation
    'resonance_weak_feed_down_flag': 1,  # include weak feed down contribution
    'collect_neutral_particles': 1,      # Collect neutral particles for net-baryon number checks.
}
