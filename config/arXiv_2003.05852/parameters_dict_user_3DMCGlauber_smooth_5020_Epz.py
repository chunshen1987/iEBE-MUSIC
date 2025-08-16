#!/usr/bin/env python3
"""
    This script contains all the user modified parameters in
    the iEBE-MUSIC package.
"""


# control parameters
control_dict = {
    'initial_state_type': "3DMCGlauber_consttau",  # 3DMCGlauber, IPGlasma
    'walltime': "100:00:00",  # walltime to run
    'save_hydro_surfaces': True,   # flag to save hydro surfaces
    'save_UrQMD_files': False,      # flag to save UrQMD files
}


# 3DMCGlauber model
mcglauber_dict = {
    'database_name': "initialcondition_database/MCGlbPbPb5020_sigmaNN_gauss_d0.9_withMultFluct",  # path for initial conditions
}


# MUSIC
music_dict = {
    'Initial_profile': 111,     # type of initial condition 
                                # 111: 3dMCGlauber smooth initial condition (arXiv:2003.05852)
                                # 112: generalization of arXiv:2003.05852 (arXiv:2106.08125)
                                # 113: event-by-event TA TB version of 112 (arXiv:2203.15718)
    'Initial_TA_Distribution_Filename': 'initial/initial_TA.dat',
    'Initial_TB_Distribution_Filename': 'initial/initial_TB.dat',
    # parameters for the eta profiles in entropy density and net baryon density
    'ecm': 5020.,                   # collision energy
    'Eta_plateau_size': 2.8,        # [-Eta_plateau_size/2, Eta_plateau_size/2] for entropy density
    'Eta_fall_off': 2.15,            # Gaussian width fall off for entropy density
    'eta_rhob_0': 6.0,              # peak position of the net baryon density
    'eta_rhob_width_1': 0.1,        # Gaussian width for |eta| > |eta_0|
    'eta_rhob_width_2': 2.0,        # Gaussian width for |eta| < |eta_0|
    'yL_frac': 0.0, 
    'Initial_time_tau_0': 1.0,  # starting time of the hydrodynamic evolution (fm/c)
    'Delta_Tau': 0.010,     # time step to use in the evolution [fm/c]
    'boost_invariant':  0,  # whether the simulation is boost-invariant
    'Eta_grid_size': 20.0,  # spatial rapidity range 
                            # [-Eta_grid_size/2, Eta_grid_size/2 - delta_eta]
    'Grid_size_in_eta': 64,  # number of the grid points in spatial rapidity direction
    'X_grid_size_in_fm': 26.0,  # spatial range along x direction in the transverse plane 
                                # [-X_grid_size_in_fm/2, X_grid_size_in_fm/2]
    'Y_grid_size_in_fm': 26.0,  # spatial range along x direction in the transverse plane 
                                # [-X_grid_size_in_fm/2, X_grid_size_in_fm/2]
    'Grid_size_in_x': 261,      # number of the grid points in x direction
    'Grid_size_in_y': 261,      # number of the grid points in y direction
    'EOS_to_use': 14,           # type of the equation of state
                                # 14: neos_BQS lattice EoS at finite mu_B
                                # 17: BEST lattice EoS at finite mu_B
    # transport coefficients
    'quest_revert_strength': 1.0,          # the strength of the viscous regulation
    'Viscosity_Flag_Yes_1_No_0': 1,        # turn on viscosity in the evolution
    'Include_Shear_Visc_Yes_1_No_0': 1,    # include shear viscous effect
    'Shear_to_S_ratio': 0.10,              # value of \eta/s
    'T_dependent_Shear_to_S_ratio': 1,     # flag to use temperature dep. \eta/s(T)
    'muB_dependent_Shear_to_S_ratio': 10,  # flag to use temperature dep. \eta/s(T)
    'Include_Bulk_Visc_Yes_1_No_0': 0,     # include bulk viscous effect
    'Include_second_order_terms': 1,       # include second order non-linear coupling terms
    'Include_vorticity_terms': 1,          # include vorticity coupling terms
    'Include_Rhob_Yes_1_No_0': 1,
    'turn_on_baryon_diffusion': 0,
    'kappa_coefficient': 0.4,
    
    # switches to output evolution information
    'output_evolution_data': 2,     # flag to output evolution history to file
    'output_movie_flag': 0,
    'output_evolution_T_cut': 0.145,
    'outputBinaryEvolution': 1,     # output evolution file in binary format
    'output_evolution_every_N_eta': 1,  # output evolution file every Neta steps
    'output_evolution_every_N_x':  2,   # output evolution file every Nx steps
    'output_evolution_every_N_y': 2,    # output evolution file every Ny steps
    'output_evolution_every_N_timesteps': 10,  # output evolution every Ntime steps

    # parameters for freeze out and Cooper-Frye
    'N_freeze_out': 1,
    'eps_freeze_max': 0.2,
    'eps_freeze_min': 0.2,
}


# iSS
iss_dict = {
    'hydro_mode': 2,    # mode for reading in freeze out information 
    'include_deltaf_shear': 1,      # include delta f contribution from shear
    'include_deltaf_bulk': 0,       # include delta f contribution from bulk
    'include_deltaf_diffusion': 0,  # include delta f contribution from diffusion
    'sample_upto_desired_particle_number': 1,  # 1: flag to run sampling until desired
                                               # particle numbers is reached
    'number_of_particles_needed': 100000,      # number of hadrons to sample
    'local_charge_conservation': 0,  # flag to impose local charge conservation
    'global_momentum_conservation': 0,  # flag to impose GMC
}


# hadronic afterburner toolkit
hadronic_afterburner_toolkit_dict = {
    'event_buffer_size': 100000,       # the number of events read in at once
    'compute_correlation': 0,       # flag to compute correlation function
    'flag_charge_dependence': 0,    # flag to compute charge dependence correlation
    'resonance_weak_feed_down_flag': 0,  # include weak feed down contribution
}
