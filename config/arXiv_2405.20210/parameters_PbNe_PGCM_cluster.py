#!/usr/bin/env python3
"""
    This script contains all the default parameters in the iEBE-MUSIC package.
"""

# control parameters
control_dict = {
    'initial_state_type': "3DMCGlauber_dynamical",  # options: IPGlasma, IPGlasma+KoMPoST,
                                                    #          3DMCGlauber_dynamical, 3DMCGlauber_consttau
    'walltime': "120:00:00",  # walltime to run
}



# 3DMCGlauber model
mcglauber_dict = {
    'database_name': "self",     # self: generate initial condition on the fly
    'Projectile':  "Ne20",       # projectile nucleus name
    'nucleon_configuration_from_file': 1,
    'lightNucleusOption': 2,
    'Target'    :  "Pb",         # target nucleus name
    'roots'     :   68.5,        # collision energy (GeV)
    'seed'      :   -1,          # random seed (-1: system)
    'baryon_junctions': 1,       # 0: baryon number assumed to be at string end
                                 # 1: baryon number transported assuming baryon
                                 # junctions (at smaller x)
                                 # see arXiv:nucl-th/9602027
    'lambdaB': 0.20,             # parameter the controls the strength of
                                 # the baryon junction stopping
    'BG': 5.,
    'shadowing_factor': 0.45,     # a shadowning factor for producing strings from multiple scatterings
    'rapidity_loss_method': 3,
    'remnant_energy_loss_fraction': 0.5,     # nucleon remnants energy loss fraction (fraction of string’s y_loss) [0, 1]
    'yloss_param_slope': 1.32,               # the slope parameter for yloss parameterization [0., 1.]
    'yloss_param_alpha1': 1.8,               # the small y ~ y^alpha1 for yloss parameterization (>=1.)
    'yloss_param_alpha2': 0.34,              # the large y ~ y^alpha2 for yloss parameterization [0., 1.]
    'yloss_param_fluct_var_RHIC': 0.6,
    'yloss_param_fluct_var_LHC': 0.6,
    'evolve_QCD_string_mode': 4,        # string evolution mode
                                        # 1: deceleration with fixed rapidity loss (m/sigma = 1 fm, dtau = 0.5 fm)
                                        # 2: deceleration with LEXUS sampled rapidit loss (both dtau and sigma fluctuate)
                                        # 3: deceleration with LEXUS sampled rapidit loss (m/sigma = 1 fm, dtau fluctuates)
                                        # 4: deceleration with LEXUS sampled rapidit loss (dtau = 0.5 fm, m/sigma fluctuates)
}


# MUSIC
music_dict = {
    'beastMode': 2,
    'Initial_profile': 13,  # type of initial condition 
                            # 13: dynamical initialization (3dMCGlauber_dynamical)
                            #   -- 131: 3dMCGlauber with zero nucleus thickness

    # parameters for Initial_profile == 13 or 131
    'string_source_sigma_x': 0.2,   # the transverse size of the hotspot [fm]
    'string_source_sigma_eta': 0.6, # the smearning size of the hotspot in eta
    'stringTransverseShiftFrac': 0.0,  # control the shift of transverse coord as a function of eta for string
    'stringPreEqFlowFactor': 0.15,

    # read in initial conditions from external file (Initial_profile == 9x)
    's_factor': 1.0,      # normalization factor read in initial data file

    'Initial_time_tau_0': 0.4,          # starting time of the hydrodynamic evolution (fm/c)
    'Delta_Tau': 0.015,                 # time step to use in the evolution [fm/c]
    'Total_evolution_time_tau': 30.,    # the maximum allowed running evolution time (fm/c)

    'boost_invariant':  0,      # whether the simulation is boost-invariant 
    'Eta_grid_size': 16.0,      # spatial rapidity range
                                # [-Eta_grid_size/2, Eta_grid_size/2 - delta_eta]
    'Grid_size_in_eta': 64,      # number of the grid points in spatial rapidity direction
    'X_grid_size_in_fm': 30.0,  # spatial range along x direction in the transverse plane
                                # [-X_grid_size_in_fm/2, X_grid_size_in_fm/2]
    'Y_grid_size_in_fm': 30.0,  # spatial range along x direction in the transverse plane
                                # [-X_grid_size_in_fm/2, X_grid_size_in_fm/2]
    'Grid_size_in_x': 128,      # number of the grid points in x direction
    'Grid_size_in_y': 128,      # number of the grid points in y direction

    'EOS_to_use': 14,            # type of the equation of state
                                # 0: ideal gas
                                # 1: EOS-Q from azhydro
                                # 7: lattice EOS s95p-v1.2 for UrQMD
                                # 9: hotQCD EOS with UrQMD
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

    # parameters for freeze out and Cooper-Frye 
    'Do_FreezeOut_Yes_1_No_0': 1,       # flag to find freeze-out surface
    'freeze_out_tau_start_max': 2,      # the maximum freeze-out starting time [fm/c]
    'eps_switch': 0.45,
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
    'maximum_sampling_events': 10000,
    'local_charge_conservation': 0,     # flag to impose local charge conservation
    'global_momentum_conservation': 0,  # flag to impose GMC
}


# hadronic afterburner toolkit
hadronic_afterburner_toolkit_dict = {
    'analyze_flow': 1,                   # 0/1: flag to perform flow analysis
    'rapidity_shift': -4.2945,           # The rapidity shift
    'rapidityPTDistributionFlag': 0,     # output Qn vectors in (eta, pT)
    'event_buffer_size': 500000,         # the number of events read in at once
    'compute_correlation': 0,            # flag to compute correlation function
    'flag_charge_dependence': 0,         # flag to compute charge dependence correlation
    'compute_corr_rap_dep': 0,           # flag to compute the rapidity dependent multi-particle correlation
    'resonance_weak_feed_down_flag': 0,  # include weak feed down contribution
}

