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
    'Projectile':  "dipole",         # projectile nucleus name
    'Target'    :  "Pb",         # target nucleus name
    'resetProjWS':  0,                                                          
    'resetTargWS':  1,                                                          
    'TargWS_R':  6.68,                                                          
    'TargWS_a':  0.528,                                                          
    'TargWS_beta2':  0.0,                                                      
    'TargWS_beta3':  0.0,                                                      
    'TargWS_beta4':  0.,                                                        
    'TargWS_gamma':  0.0,
    'TargWS_da':  0.16,
    'TargWS_dR':  0.01,
    'd_min': 0.9,
    'useQuarks': 1,
    'Q2': 1.,
    'roots'     :   5020,         # collision energy (GeV)
    'seed'      :   -1,          # random seed (-1: system)
    'baryon_junctions': 1,       # 0: baryon number assumed to be at string end
                                 # 1: baryon number transported assuming baryon
                                 # junctions (at smaller x)
                                 # see arXiv:nucl-th/9602027
    'electric_junctions':  1,
    'integer_electric_charge': 1,
    'electricChargeInStringProb': 0.0, 
    'lambdaB': 1.0,              # parameter the controls the strength of
    'lambdaQ': 1.0,              # parameter the controls the strength of
                                 # the baryon junction stopping
    'lambdaBs': 1.0,             # fraction of single-to-double string junction stopping
    'lambdaQs': 1.0,             # fraction of single-to-double string junction stopping
    'baryonInStringProb': 0.1,
    'batch_density_output': 0,
    'BG': 16.,
    'BG_proj': 5.,
    'BG_targ': 17.,
    'use_roots_cut': 0,
    'use_roots_distribution': 1,
    'roots_low_cut': 2.,
    'roots_up_cut': 700.,
    'use_E_dependent_LB': 0,
    'CB': 1.58,
    'shadowing_factor':  0.20,   # a shadowning factor for producing strings from multiple scatterings
    'rapidity_loss_method': 4,
    'remnant_energy_loss_fraction': 0.5,         # nucleon remnants energy loss fraction (fraction of string's y_loss) [0, 1]
    'ylossParam4At2': 1.70,
    'ylossParam4At4': 2.00,
    'ylossParam4At6': 2.20,
    'ylossParam4At10': 2.6,
    'ylossParam4var': 0.675,
    'evolve_QCD_string_mode': 4,        # string evolution mode
                                        # 1: deceleration with fixed rapidity loss (m/sigma = 1 fm, dtau = 0.5 fm)
                                        # 2: deceleration with LEXUS sampled rapidit loss (both dtau and sigma fluctuate)
                                        # 3: deceleration with LEXUS sampled rapidit loss (m/sigma = 1 fm, dtau fluctuates)
                                        # 4: deceleration with LEXUS sampled rapidit loss (dtau = 0.5 fm, m/sigma fluctuates)
}


# MUSIC
music_dict = {
    'beastMode': 2,
    'Initial_profile': 13,   # type of initial condition 
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
    'Include_QS_Yes_1_No_0': 1, # use nq and ns in hydro evolution.
    'use_BQ_ratios': 0, # use nq = 0.4 nb in hydro source strings.
    'turn_on_baryon_diffusion': 0,
    'kappa_coefficient': 0.4,

    # parameters for freeze out and Cooper-Frye
    'freeze_out_tau_start_max': 2,      # the maximum freeze-out starting time [fm/c]
    'N_freeze_out': 1,
    'eps_switch': 0.2,                 # GeV/fm^3
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
    'regulateEOS': 0,          # flag to regulate T and mu with pure HRG EOS
    'maximum_sampling_events': 10000,

}

hadronic_afterburner_toolkit_dict = {
    'rap_shift': 0.465,                  # The rapidity shift
    'rap_type_for_pid': 1,               # 0: for pseudo-rapidity; 1: for rapidity
    'rap_min_for_pid': -1.0,             # minimum value of rapidity integration
                                         # range for mid-rapidity observables
    'rap_max_for_pid': 0.0,              # maximum value of rapidity integration
    'output_y_pt_spectra': 1,
    'y0': -0.7,
    'y1': -0.5,
    'y2': -0.3,
    'y3': -0.1,
    'y4':  0.1,
    'y5':  0.3,
    'y6':  0.5,
    'y7':  0.7,
    'event_buffer_size': 500000,         # the number of events read in at once
    'compute_correlation': 1,            # flag to compute correlation function
    'flag_charge_dependence': 1,         # flag to compute charge dependence correlation
    'compute_corr_rap_dep': 0,           # flag to compute the rapidity dependent multi-particle correlation
    'resonance_weak_feed_down_flag': 0,  # include weak feed down contribution
}

