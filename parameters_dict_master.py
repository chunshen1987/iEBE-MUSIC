#!/usr/bin/env python3
"""
    This script contains all the default parameters in the iEBE-MUSIC package.
"""

# control parameters
control_dict = {
    'initial_state_type': "3DMCGlauber_dynamical",  # options: IPGlasma, IPGlasma+KoMPoST,
                                                    #          3DMCGlauber_dynamical, 3DMCGlauber_consttau
    'walltime': "120:00:00",  # walltime to run
    'save_ipglasma_results': False,   # flag to save IPGlasma results
    'save_kompost_results': False,    # flag to save kompost results
    'save_hydro_surfaces': False,     # flag to save hydro surfaces
    'save_UrQMD_files': False,        # flag to save UrQMD files
    'compute_photon_emission': False,   # flag to compute EM radiation from hydrodynamic medium
    'compute_polarization': False,       # flag to save spin polarization results
}



# 3DMCGlauber model
mcglauber_dict = {
    'database_name': "self",     # self: generate initial condition on the fly
    'Projectile':  "dipole",          # projectile nucleus name
    'Target'    :  "Au",         # target nucleus name
    'roots'     :   54.4,        # collision energy (GeV)
    'lightNucleusOption':   2,   # O-16: 4. profile from PGCM
    'seed'      :   -1,          # random seed (-1: system)
    'baryon_junctions': 1,       # 0: baryon number assumed to be at string end
                                 # 1: baryon number transported assuming baryon
                                 # junctions (at smaller x)
                                 # see arXiv:nucl-th/9602027
    'lambdaB': 0.40,              # parameter the controls the strength of
                                 # the baryon junction stopping
    'use_roots_distribution': 1, 
    'BG_proj': 16.,
    'BG_targ': 16.,
    'BG': 16.,
    'use_roots_cut': 0,
    'roots_low_cut': 600.,
    'roots_up_cut': 700.,
    'use_E_dependent_LB': 1,
    'CB': 1.58,
    "remnant_x_is_ori": 1,
    'shadowing_factor': 0.25,     # a shadowning factor for producing strings from multiple scatterings
    'rapidity_loss_method': 3,
    'remnant_energy_loss_fraction': 0.5,     # nucleon remnants energy loss fraction (fraction of string’s y_loss) [0, 1]
    'yloss_param_slope': 1.32,               # the slope parameter for yloss parameterization [0., 1.]
    'yloss_param_alpha1': 1.8,               # the small y ~ y^alpha1 for yloss parameterization (>=1.)
    'yloss_param_alpha2': 0.35,              # the large y ~ y^alpha2 for yloss parameterization [0., 1.]
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
    'echo_level':  1,       # control the mount of message output to screen
    'beastMode': 1,
    'mode': 2,              # MUSIC running mode 2: Evolution only.
    'Initial_profile': 13,   # type of initial condition 
                            # 9: IPGlasma (full Tmunu),
                            #   -- 91: e and u^\mu,
                            #   -- 92: e only,
                            #   -- 93: e, u^\mu, and pi^\munu
                            # 11: 3dMCGlauber initial condition at a constant tau surface
                            #     based on the nuclear thickness funciton TA and TB
                            #   -- 111: second parameterization of eta profile
                            # 13: dynamical initialization (3dMCGlauber_dynamical)
                            #   -- 131: 3dMCGlauber with zero nucleus thickness

    # parameters for the eta profiles in entropy density and net baryon density
    # Initial_profile == 11, 111, 112, 113
    'Initial_participantList_Filename': 'initial/participants_event.dat',
    'ecm': 200.,                    # collision energy
    'Eta_plateau_size': 5.4,        # [-Eta_plateau_size/2, Eta_plateau_size/2] for entropy density
    'Eta_fall_off': 0.3,            # Gaussian width fall off for entropy density
    'eta_rhob_0': 1.5,              # peak position of the net baryon density
    'eta_rhob_width_1': 0.2,        # Gaussian width for |eta| > |eta_0|
    'eta_rhob_width_2': 1.0,        # Gaussian width for |eta| < |eta_0|

    # parameters for Initial_profile == 13 or 131
    'string_source_sigma_x': 0.2,   # the transverse size of the hotspot [fm]
    'string_source_sigma_eta': 0.6, # the smearning size of the hotspot in eta
    'stringTransverseShiftFrac': 0.0,  # control the shift of transverse coord as a function of eta for string

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
    'bulk_viscosity_10_max': 0.30,         # the peak value of \zeta/s(T)
    'bulk_viscosity_10_T_peak': 0.17,      # the peak temperature for \zeta/s(T)
    'Include_second_order_terms': 1,       # include second order non-linear coupling terms
    'Include_vorticity_terms': 0,          # include vorticity coupling terms
    'Include_Rhob_Yes_1_No_0': 1,
    'turn_on_baryon_diffusion': 0,
    'kappa_coefficient': 0.4,
    'stringPreEqFlowFactor': 0.15,
    # switches to output evolution information
    'output_hydro_debug_info': 0,   # flag to output debug information
    'output_evolution_data': 0,     # flag to output evolution history to file
    'output_movie_flag': 0,
    'output_evolution_T_cut': 0.145,
    'outputBinaryEvolution': 1,     # output evolution file in binary format
    'output_evolution_every_N_eta': 1,  # output evolution file every Neta steps
    'output_evolution_every_N_x':  1,   # output evolution file every Nx steps
    'output_evolution_every_N_y': 1,    # output evolution file every Ny steps
    'output_evolution_every_N_timesteps': 10,  # output evolution every Ntime steps

    # parameters for freeze out and Cooper-Frye 
    'Do_FreezeOut_Yes_1_No_0': 1,       # flag to find freeze-out surface
    'Do_FreezeOut_lowtemp': 1,          # flag to include cold corona
    'freeze_out_tau_start_max': 2,      # the maximum freeze-out starting time [fm/c]
    'freeze_surface_in_binary': 1,      # switch to output surface file in binary format
    'average_surface_over_this_many_time_steps': 10,   # the step skipped in the tau
    'freeze_Ncell_x_step': 1,
    'freeze_Ncell_eta_step': 1,
    'N_freeze_out': 1,
    'eps_switch': 0.5,
    'eps_freeze_max': 0.5,
    'eps_freeze_min': 0.5,
    'use_eps_for_freeze_out': 1,  # find freeze-out surface 
                                  # 0: use temperature, 1: use energy density
}


# photon_emission
photon_dict = {
    'hydro_flag': 2,          # read in mode for hydro medium
                              # 0: read in hdf5 file
                              # 1: read in binary file output from MUSIC
                              # 2: read in binary file output from new MUSIC (no grid)
                              # 3: read in binary file output from new MUSIC (on grid)
    'hydro_nskip_tau': 1,     # read in hydro slice every hydro_nskip_tau
                              # steps from the medium file
                              # (only works for hydro_flag = 1)
    'Xmin': -15.0,            # minimum points along x direction
    'dx': 0.1,                # lattice spacing along x direction
    'Ymin': -15.0,            # minimum points along y direction
    'dy': 0.1,                # lattice spacing along y direction
    'tau_start': 0.4,         # emission start time (fm)
    'tau_end': 30.0,          # emission end time (fm)
    'dTau': 0.1,              # lattice spacing along tau direction

    'neta': 10,               # number of points in eta direction
    'eta_i': 0.0,             # beginning value of eta slice
    'eta_f': 3.0,             # end value of eta slice

    'np': 20,                 # number of points for photon momentum
    'nphi': 40,               # number of points for angles of photons momenta
    'nrapidity': 1,           # number of points for photon rapidity

    'photon_q_i': 0.2,        # the smallest photon momentum to be calculated
    'photon_q_f': 4.0,        # the largest photon momentum to be calculated
    'photon_phi_q_i': 0.0,    # the smallest angle of photon momentum
    'photon_phi_q_f': 6.2831853,    # the largest angle of photon momentum
    'photon_y_i': 0.0,        # the smallest photon rapidity
    'photon_y_f': 0.0,        # the largest photon rapidity

    'norder': 10,             # calculate photon vn to norder
    'turn_on_muB': 1,         # flag to include muB dependence in photon rates

    'T_dec': 0.105,           # freeze out temperature (GeV)
    'T_sw_high': 0.180,       # high end of the switching temperature
    'T_sw_low': 0.1795,       # low end of the switching temperature
    'T_cuthigh': 0.80,        # maximum allowed emission T (GeV)
    'T_cutlow': 0.10,         # minimum allowed emission T (GeV)

    'calHGIdFlag': 0,         # Flag to decide whether to calculate individual HG channels

    'PhotonemRatetableInfo_Emin': 0.05,   # minimum photon energy in the photon rate tables
    'PhotonemRatetableInfo_Tmin': 0.10,   # minimum temperature in the photon rate tables
    'PhotonemRatetableInfo_dE': 0.05,     # lattice space of energy in the photon rate tables
    'PhotonemRatetableInfo_dT': 0.002,    # lattice space of temperature in the photon rate tables

    'HydroinfoVisflag': 1,         # determine whether to read in the viscous evolution information
    'HydroinfoBuffersize': 500,    # set the buffer size for hydro evolution profile

    'turn_off_transverse_flow': 0,      # flag to turn off transverse flow in the photon calculation
    'enable_polyakov_suppression': 0,   # apply the polyakov suppression to QGP photon rates

    'differential_flag': 0,  # determine whether to output differential photon yield and vn
                             # 1: differential in T and tau
                             # 2: differential in x and tau
                             # 10: differeitial in all options above
    'nTaucut': 50,           # number of points in tau (range of tau is specified by tau_start and tau_end)
    'nTcut': 50,             # number of points in T (range of T is specified by T_cuthigh and T_cutlow)
    'n_xperp_cut': 101,      # number of points in x
    'xperp_cuthigh': 10.0,   # maximum value in x (fm)
    'xperp_cutlow': -10.0,   # minimum value in x (fm)
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
    'local_charge_conservation': 0,  # flag to impose local charge conservation
    'global_momentum_conservation': 0,  # flag to impose GMC
}


# hadronic afterburner toolkit
hadronic_afterburner_toolkit_dict = {
    'rap_shift': -1.71,                  # The rapidity shift
    'rap_type_for_pid': 0,               # 0: for pseudo-rapidity; 1: for rapidity

    'rap_min_for_pid': -1.0,             # minimum value of rapidity integration  
                                         # range for mid-rapidity observables 
    'rap_max_for_pid': 1.0,              # maximum value of rapidity integration
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
    'resonance_weak_feed_down_flag': 1,  # include weak feed down contribution
}
