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
    'Projectile':  "He3",         # projectile nucleus name
    'Target'    :  "Au",         # target nucleus name
    'roots'     :   200.,        # collision energy (GeV)
    'seed'      :   -1,          # random seed (-1: system)
    'baryon_junctions': 1,       # 0: baryon number assumed to be at string end
                                 # 1: baryon number transported assuming baryon
                                 # junctions (at smaller x)
                                 # see arXiv:nucl-th/9602027
    'lambdaB': 0.2,              # parameter the controls the strength of
                                 # the baryon junction stopping
    'BG': 5.,
    'nucleon_configuration_from_file': 1,
    'shadowing_factor': 0.45,     # a shadowning factor for producing strings from multiple scatterings
    'rapidity_loss_method': 3,
    'remnant_energy_loss_fraction': 0.5,         # nucleon remnants energy loss fraction (fraction of string's y_loss) [0, 1]
    'yloss_param_slope': 1.32,         # the slope parameter for yloss parameterization [0., 1.]
    'yloss_param_alpha1': 1.8,         # the small y ~ y^alpha1 for yloss parameterization (>=1.)
    'yloss_param_alpha2': 0.34,        # the large y ~ y^alpha2 for yloss parameterization [0., 1.]
    'yloss_param_fluct_var_RHIC': 0.4,
    'yloss_param_fluct_var_LHC': 0.4,
    'evolve_QCD_string_mode': 4,        # string evolution mode
                                        # 1: deceleration with fixed rapidity loss (m/sigma = 1 fm, dtau = 0.5 fm)
                                        # 2: deceleration with LEXUS sampled rapidit loss (both dtau and sigma fluctuate)
                                        # 3: deceleration with LEXUS sampled rapidit loss (m/sigma = 1 fm, dtau fluctuates)
                                        # 4: deceleration with LEXUS sampled rapidit loss (dtau = 0.5 fm, m/sigma fluctuates)
}


# MUSIC
music_dict = {
    'Initial_profile': 13,   # type of initial condition 
                             # 13: dynamical initialization (3dMCGlauber_dynamical)
                             #   -- 131: 3dMCGlauber with zero nucleus thickness
    'string_source_sigma_eta': 0.6,  # the smearning size of the hotspot in eta
    'string_source_sigma_x': 0.2,    # the smearning size of the hotspot in eta
    's_factor': 1.000,      # normalization factor read in initial data file
    'Delta_Tau': 0.015,     # time step to use in the evolution [fm/c]
    'boost_invariant':  0,  # whether the simulation is boost-invariant
    'Eta_grid_size': 16.0,  # spatial rapidity range 
                            # [-Eta_grid_size/2, Eta_grid_size/2 - delta_eta]
    'Grid_size_in_eta': 64,     # number of the grid points in spatial rapidity direction
    'X_grid_size_in_fm': 30,    # spatial range along x direction in the transverse plane 
                                # [-X_grid_size_in_fm/2, X_grid_size_in_fm/2]
    'Y_grid_size_in_fm': 30,    # spatial range along x direction in the transverse plane 
                                # [-X_grid_size_in_fm/2, X_grid_size_in_fm/2]
    'Grid_size_in_x': 128,      # number of the grid points in x direction
    'Grid_size_in_y': 128,      # number of the grid points in y direction
    'EOS_to_use': 14,           # type of the equation of state
                                # 14: neos_BQS lattice EoS at finite mu_B
                                # 17: BEST lattice EoS at finite mu_B
    # transport coefficients
    'quest_revert_strength': 1.0,
    'Viscosity_Flag_Yes_1_No_0': 1,        # turn on viscosity in the evolution
    'Include_Shear_Visc_Yes_1_No_0': 1,    # include shear viscous effect
    'Shear_to_S_ratio': 0.13,              # value of \eta/s
    'T_dependent_Shear_to_S_ratio': 0,     # flag to use temperature dep. \eta/s(T)
    'muB_dependent_Shear_to_S_ratio': 10,
    'shear_muBDep_alpha': 1.5,
    'shear_muBDep_slope': 1.0,
    'shear_muBDep_scale': 0.6,             # GeV
    'Include_Bulk_Visc_Yes_1_No_0': 1,     # include bulk viscous effect
    'T_dependent_zeta_over_s': 10,         # parameterization of \zeta/s(T)
    'bulk_viscosity_10_max': 0.08,         # the peak value of \zeta/s(T)
    'bulk_viscosity_10_T_peak': 0.17,      # the peak temperature for \zeta/s(T)
    'Include_second_order_terms': 1,       # include second order non-linear coupling terms
    'Include_vorticity_terms': 0,          # include vorticity coupling terms
    'Include_Rhob_Yes_1_No_0': 1,
    'turn_on_baryon_diffusion': 0,
    'kappa_coefficient': 0.4,
    'stringPreEqFlowFactor': 0.15,

    # parameters for freeze out and Cooper-Frye
    'freeze_out_tau_start_max': 2,      # the maximum freeze-out starting time [fm/c]
    'N_freeze_out': 1,
    'eps_switch': 0.5, 
}


# iSS
iss_dict = {
    'hydro_mode': 2,    # mode for reading in freeze out information 
    'MC_sampling': 4,
    'include_deltaf_shear': 1,      # include delta f contribution from shear
    'include_deltaf_bulk': 1,       # include delta f contribution from bulk
    'bulk_deltaf_kind': 20,         # 20: 22-momentum approximation, 21: CE relaxation time approximation
    'include_deltaf_diffusion': 0,  # include delta f contribution from diffusion
    'sample_upto_desired_particle_number': 1,  # 1: flag to run sampling until desired
                                               # particle numbers is reached
    'number_of_particles_needed': 100000,      # number of hadrons to sample
    'local_charge_conservation': 0,  # flag to impose local charge conservation
    'global_momentum_conservation': 0,  # flag to impose GMC
    'maximum_sampling_events': 10000, 
}


# hadronic afterburner toolkit
hadronic_afterburner_toolkit_dict = {
    'event_buffer_size': 500000,         # the number of events read in at once
    'compute_correlation': 0,            # flag to compute correlation function
    'flag_charge_dependence': 0,         # flag to compute charge dependence correlation
    'compute_corr_rap_dep': 0,           # flag to compute the rapidity dependent multi-particle correlation
    'resonance_weak_feed_down_flag': 1,  # include weak feed down contribution
}

