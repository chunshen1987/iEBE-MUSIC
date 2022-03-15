#!/usr/bin/env python3
"""
    This script contains all the default parameters in the iEBE-MUSIC package.
"""

from os import path, makedirs
import sys
import shutil
import argparse

# control parameters
control_dict = {
    'initial_state_type': "3DMCGlauber_dynamical",  # options: IPGlasma, IPGlasma+KoMPoST,
                                                    #          3DMCGlauber_dynamical, 3DMCGlauber_consttau
    'walltime': "10:00:00",  # walltime to run
    'save_ipglasma_results': False,   # flag to save IPGlasma results
    'save_kompost_results': False,    # flag to save kompost results
    'save_hydro_surfaces': False,     # flag to save hydro surfaces
    'save_UrQMD_files': False,        # flag to save UrQMD files
    'compute_photon_emission': False,   # flag to compute EM radiation from hydrodynamic medium
    'compute_polarization': False,       # flag to save spin polarization results
}


# IPGlasma
ipglasma_dict = {
    'type': "self",  # minimumbias or fixed (pre-generated)
                     # self (generate on the fly)

    # path for the database file (for type == minimumbias or fixed)
    'database_name_pattern': "IPGlasma_database/AuAu_C{0:s}.h5",

    # all parameters below are for (type == self)
    'mode': 1,              # run mode
    'readMultFromFile': 0,
    'size': 720,            # number of grid points of IP-Glasma computation
    'L': 30.,               # grid size in the transverse plane
    'Nc': 3,                # number of color
    'm': 0.2,               # infrared cut-off mass (GeV)
    'rmax': 10.,
    'UVdamp': 0.,
    'Jacobianm': 0.35,
    'g': 1.,                # strong coupling constant
    'BG': 4.,
    'BGq': 0.3,
    'BGqVar': 0.0,
    'dqMin': 0.0,
    'useSmoothNucleus': 0,
    'useConstituentQuarkProton': 0,
    'NqFluc': 0.0,
    'shiftConstituentQuarkProtonOrigin': 1,
    'runningCoupling': 0,
    'muZero': 0.3,
    'minimumQs2ST': 0.,
    'beta2': 0.28,
    'c': 0.2,
    'g2mu': 0.1,
    'useFatTails': 0,
    'tDistNu': 3,
    'smearQs': 1,
    'smearingWidth': 0.6,
    'protonAnisotropy': 0,
    'roots': 200.,
    'usePseudoRapidity': 0,
    'Rapidity': 0.,
    'useFluctuatingx': 1,
    'xFromThisFactorTimesQs': 1,
    'useNucleus': 1,
    'useGaussian': 0,
    'nucleonPositionsFromFile': 0,
    'NucleusQsTableFileName': "qs2Adj_vs_Tp_vs_Y_200.in",
    'QsmuRatio': 0.8,
    'samplebFromLinearDistribution': 1,
    'runWith0Min1Avg2MaxQs': 2,
    'runWithThisFactorTimesQs': 0.5,
    'runWithLocalQs': 0,
    'runWithkt': 0,
    'Ny': 50,
    'useSeedList': 0,
    'seed': 3,
    'useTimeForSeed': 0,
    'Projectile': "Au",
    'Target': "Au",
    'bmin': 0.,
    'bmax': 20.,
    'lightNucleusOption': 1,
    'useFixedNpart': 0,
    'averageOverThisManyNuclei': 1,
    'SigmaNN': 42.,
    'gaussianWounding': 1,
    'inverseQsForMaxTime': 0,
    'maxtime': 0.6,
    'dtau': 0.1,
    'LOutput': 30,
    'sizeOutput': 512,
    'etaSizeOutput': 1,
    'detaOutput': 0,
    'writeOutputs': 5,
    'writeEvolution': 0,
    'writeInitialWilsonLines': 0,
    'writeOutputsToHDF5': 0
}


# 3DMCGlauber model
mcglauber_dict = {
    'database_name': "self",     # self: generate initial condition on the fly
    'Projectile':  "Pb",         # projectile nucleus name
    'Target'    :  "Pb",         # target nucleus name
    'nucleon_configuration_from_file': 0,
    'roots'     :   17.3,        # collision energy (GeV)
    'useQuarks' :   1,           # switch to use valence quarks
    'Q2'        :   1.,          # the scale when evaluating the pdf
    'b_min'     :   0.,          # minimum impact parameter (fm)
    'b_max'     :   20.,         # maximum impact parameter (fm)
    'seed'      :   -1,          # random seed (-1: system)
    'only_event_statistics': 0,  # flag to only output the event_summary file
    'cache_tables': 1,           # 1: use pre-generated tables for valence quark x
                                 # 0: re-generate tables for valence quark x
    'baryon_junctions': 0,       # 0: baryon number assumed to be at string end
                                 # 1: baryon number transported assuming baryon
                                 # junctions (at smaller x)
                                 # see arXiv:nucl-th/9602027
    'lambdaB': 0.2,              # parameter the controls the strength of
                                 # the baryon junction stopping
    'BG': 4.,                    # Gaussian width for sampling the valence quark positions
    'shadowing_factor': 1.0,     # a shadowning factor for producing strings from multiple scatterings
    'fluct_Nstrings_per_NN_collision': 1,        # fluctuate number of strings produced per NN collision
    'QCD_string_production_mode': 1,    # string production mode
                                        # 1: strings are produced randomly in the binary collision list
                                        # 2: strings are produced at the last binary collision
    'rapidity_loss_method': 2,          # 1: LEXUS
                                        # 2: parameterization
                                        # 3: parameterization with logit-normal fluctuation
    'remnant_energy_loss_fraction': 0.5,         # nucleon remnants energy loss fraction (fraction of string's y_loss) [0, 1]
    'yloss_param_slope': 1.50,          # the slope parameter for yloss parameterization [0., 1.]
    'yloss_param_alpha1': 2.50,         # the small y ~ y^alpha1 for yloss parameterization (>=1.)
    'yloss_param_alpha2': 0.25,         # the large y ~ y^alpha2 for yloss parameterization [0., 1.]
    'yloss_param_fluct_var_RHIC': 0.60, # the variance ofthe logit-normal parameterized y_loss fluctuation
    'yloss_param_fluct_var_LHC': 0.80,  # the variance ofthe logit-normal parameterized y_loss fluctuation
    'evolve_QCD_string_mode': 2,        # string evolution mode
                                        # 1: deceleration with fixed rapidity loss (m/sigma = 1 fm, dtau = 0.5 fm)
                                        # 2: deceleration with LEXUS sampled rapidit loss (both dtau and sigma fluctuate)
                                        # 3: deceleration with LEXUS sampled rapidit loss (m/sigma = 1 fm, dtau fluctuates)
                                        # 4: deceleration with LEXUS sampled rapidit loss (dtau = 0.5 fm, m/sigma fluctuates)
    'tau_form_mean': 0.5,
    'tau_form_fluct_gamma_beta': 1.0,
}

# KoMPoST
kompost_dict = {
    'KoMPoSTInputs': {
        'tIn': 0.1,
        'tOut': 0.8,
        'InputFile': "Tmunu.dat",
        'OutputFileTag': "ekt_tIn01_tOut08",
    },
    'KoMPoSTParameters': {
        'EtaOverS': 0.12,                   # specific shear viscosity
        'EtaOverSTemperatureScale': 0.1,
        'EVOLUTION_MODE': 1,                # 0 for free-streaming, 1: for "KoMPoST" EKT evolution
        'ENERGY_PERTURBATIONS': 1,
        'MOMENTUM_PERTURBATIONS': 1,
        'DECOMPOSITION_METHOD': 1,
    },
    'EventInput': {
        'normFactor': 0.287,    # Tmunu is normalized by this factor after being read in
        'afm': 0.0664062,       # lattice spacing in fm
        'Ns': 512,              # number of grid points on a square lattice
        'xSTART': 0,            # The first grid point to include in the x direction
        'xEND': 511,            # The last grid point to include in the x direction
        'ySTART': 0,            # The first grid point to include in the y direction
        'yEND': 511,            # The last grid point to include in the y direction
    },
}


# MUSIC
music_dict = {
    'echo_level':  1,       # control the mount of message output to screen
    'mode': 2,              # MUSIC running mode 2: Evolution only.
    'Initial_profile': 9,   # type of initial condition 
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
    'string_source_sigma_x': 0.5,   # the transverse size of the hotspot [fm]
    'string_source_sigma_eta': 0.5, # the smearning size of the hotspot in eta

    # read in initial conditions from external file (Initial_profile == 9x)
    'Initial_Distribution_input_filename': 'initial/epsilon-u-Hydro.dat',
    's_factor': 0.190,      # normalization factor read in initial data file

    'Initial_time_tau_0': 0.4,          # starting time of the hydrodynamic evolution (fm/c)
    'Delta_Tau': 0.005,                 # time step to use in the evolution [fm/c]
    'Total_evolution_time_tau': 30.,    # the maximum allowed running evolution time (fm/c)

    'boost_invariant':  1,      # whether the simulation is boost-invariant 
    'Eta_grid_size': 14.0,      # spatial rapidity range
                                # [-Eta_grid_size/2, Eta_grid_size/2 - delta_eta]
    'Grid_size_in_eta': 1,      # number of the grid points in spatial rapidity direction
    'X_grid_size_in_fm': 20.0,  # spatial range along x direction in the transverse plane
                                # [-X_grid_size_in_fm/2, X_grid_size_in_fm/2]
    'Y_grid_size_in_fm': 20.0,  # spatial range along x direction in the transverse plane
                                # [-X_grid_size_in_fm/2, X_grid_size_in_fm/2]
    'Grid_size_in_x': 200,      # number of the grid points in x direction
    'Grid_size_in_y': 200,      # number of the grid points in y direction

    'EOS_to_use': 9,            # type of the equation of state
                                # 0: ideal gas
                                # 1: EOS-Q from azhydro
                                # 7: lattice EOS s95p-v1.2 for UrQMD
                                # 9: hotQCD EOS with UrQMD
                                # 14: neos_BQS lattice EoS at finite mu_B
                                # 17: BEST lattice EoS at finite mu_B
    # transport coefficients
    'quest_revert_strength': 10.0,         # the strength of the viscous regulation
    'Viscosity_Flag_Yes_1_No_0': 1,        # turn on viscosity in the evolution
    'Include_Shear_Visc_Yes_1_No_0': 1,    # include shear viscous effect
    'Shear_to_S_ratio': 0.12,              # value of \eta/s
    'T_dependent_Shear_to_S_ratio': 0,     # flag to use temperature dep. \eta/s(T)
    'muB_dependent_Shear_to_S_ratio': 1,   # flag to use temperature dep. \eta/s(T, muB)
    'Include_Bulk_Visc_Yes_1_No_0': 1,     # include bulk viscous effect
    'T_dependent_zeta_over_s': 7,          # parameterization of \zeta/s(T)
    'Include_second_order_terms': 1,       # include second order non-linear coupling terms
    'Include_vorticity_terms': 0,          # include vorticity coupling terms
    'Include_Rhob_Yes_1_No_0': 0,
    'turn_on_baryon_diffusion': 0,
    'kappa_coefficient': 0.4,

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
    'freeze_out_method': 4,             # method for hyper-surface finder
                                        # 4: Cornelius
    'freeze_surface_in_binary': 1,      # switch to output surface file in binary format
    'average_surface_over_this_many_time_steps': 10,   # the step skipped in the tau
    'freeze_Ncell_x_step': 1,
    'freeze_Ncell_eta_step': 1,
    'freeze_eps_flag': 0,
    'N_freeze_out': 1,
    'eps_freeze_max': 0.18,
    'eps_freeze_min': 0.18,
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
    'afterburner_type': 1,     # 0: PDG_Decay, 1: UrQMD, 2: SMASH
    'turn_on_bulk': 0,  # read in bulk viscous pressure
    'turn_on_rhob': 0,  # read in net baryon chemical potential
    'turn_on_diff': 0,  # read in baryon diffusion current

    'include_deltaf_shear': 1,      # include delta f contribution from shear
    'include_deltaf_bulk': 1,       # include delta f contribution from bulk
    'include_deltaf_diffusion': 0,  # include delta f contribution from diffusion

    'bulk_deltaf_kind': 1,     # 0: 14-momentum approximation, 1: relaxation time approximation
    'restrict_deltaf': 0,      # flag to apply restriction on the size of delta f
    'deltaf_max_ratio': 1.0,   # the maximum allowed size of delta f w.r.t f0
    'quantum_statistics': 1,   # include quantum statistics (1: yes, 0: no)

    'calculate_polarization': 0,   # switch to compute Lambda's polarization
    'polarizationRapType': 1,      # 0: rapidity; 1: pseudorapidity; 2: both

    'randomSeed': -1,   # If <0, use system clock.
    'calculate_vn': 0,  # 1/0: whether to calculate the 
    'MC_sampling': 4,   # 0/1/2/3: whether to perform Monte-Carlo sampling
                        # (not required for spectra calculation). 
                        # 0: No sampling. 
                        # 1: use dN_dxtdetady to sample. 
                        # 2: use dN_dxtdy to sample. 
                        # 3: use dN_pTdpTdphidy to sample 
                        #    (overwrites calculate_vn to be 1). 
    'sample_upto_desired_particle_number': 0,  # 1: flag to run sampling until desired
                                               # particle numbers is reached
    'number_of_repeated_sampling': 10,         # number of repeaded sampling
    'number_of_particles_needed': 100000,      # number of hadrons to sample

    'sample_y_minus_eta_s_range': 4,    # y_minus_eta_s will be sampled
    'sample_pT_up_to': -1,  # Up to this value will pT be sampled; 
                            # if<0 then use the largest value in the pT table.
    'dN_dy_sampling_model': 30, # 30: Use Poisson distribution to sample the 
                                #      whole dN_dy.
    'dN_dy_sampling_para1': 0.16,  # Additional parameters for dN/dy sampling. 
                                   # -- For dN_dy_sampling_model==10 or 20, 

    'perform_decays': 0,             # flag to perform resonance decay
    'perform_checks': 0,             # flag to perform tests for the sampler
    'include_spectators': 0,         # flag to include spectators
    'local_charge_conservation': 0,  # flag to impose local charge conservation
    'global_momentum_conservation': 0,  # flag to impose GMC

    'y_LB': -5.0,       # lower bound for y-sampling; 
    'y_RB': 5.0,        # upper bound for y-sampling;

    'output_samples_into_files': 0,  # output particle samples into individual files
    'store_samples_in_memory': 1,    # flag to store particle samples in memory
    'use_OSCAR_format': 1,           # output results in OSCAR format
    'use_gzip_format': 0,  # output results in gzip format (only works with
                           # store_samples_in_memory = 1)
    'use_binary_format': 0,
    'calculate_vn_to_order': 9,     # v_n's are calculated up to this order
    'use_pos_dN_only': 0,           # 1: all negative emission functions will be skipped. 
    'grouping_particles': 1,  # 0/1: Particles will be re-order according to 
                              # their mass. This parameter combined with 
                              # grouping_tolerance parameter can make particles 
                              # with similar mass and chemical potentials to be 
                              # sampled together.
    'grouping_tolerance': 0.01,  # If two particles adjacent in the table have 
                                 # mass and chemical potentials close within this 
                                 # relative tolerance, they are considered to be 
                                 # identical and will be sampled successively 
                                 # without regenerating the dN / (dxt deta dy) 
                                 # matrix for efficiency.
    'minimum_emission_function_val': 1e-30,  # If dN/(dx_t deta dy) is evaluated to 
                                             # be smaller than this value, then it 
                                             # is replaced by this value.
    'use_historic_flow_output_format': 0,
    'eta_s_LB': -0.5,      # lower bound for eta_s sampling; used only when 
                           # sampling using total energy flux
    'eta_s_RB': 0.5,       # upper bound for eta_s sampling.
    'use_dynamic_maximum': 0,   # 0/1: Whether to automatically reduce the 
                                # guessed maximum after some calculations. 
                                # Work only when MC_sampling is set to 2.
    'adjust_maximum_after': 100000,  # Used only when use_dynamic_maximum=1. 
                                     # After the number of sampling given by 
                                     # this parameter the guessed maximum is 
                                     # adjusted.
    'adjust_maximum_to': 1.2,   # [1,inf]: When guessed maximum is adjusted, 
                                # it is adjusted to the "observed maximum" 
                                # multiplied by this value. Note that the 
                                # "observed maximum" is measured relative to 
                                # the guessed maximum. See code for details.
    'calculate_dN_dtau': 0,     # Output dN_dtau table. Only applicable 
                                # if MC_sampling parameter is set to 1.
    'bin_tau0': 0.6,            # used to generate bins for 
                                # calculate_dN_dtau_using_dN_dxtdeta function
    'bin_dtau': 0.2,            # used to generate bins for 
                                # calculate_dN_dtau_using_dN_dxtdeta function
    'bin_tau_max': 17.0,        # used to generate bins for 
                                # calculate_dN_dtau_using_dN_dxtdeta function
    'calculate_dN_dx': 0,       # Output dN_dx table. Only applicable 
                                # if MC_sampling parameter is set to 1.
    'bin_x_min': -10.0,         # used to generate bins for 
                                # calculate_dN_dx_using_dN_dxtdeta function
    'bin_dx': 0.5,              # used to generate bins 
                                # for calculate_dN_dx_using_dN_dxtdeta function
    'bin_x_max': 10.0,          # used to generate bins for 
                                # calculate_dN_dx_using_dN_dxtdeta function
    'calculate_dN_dphi': 0,     # Output dN_dphi table. Only applicable 
                                # if calculate_vn parameter is set to 1.
    'calculate_dN_deta':1,      # Output dN_deta table. Only applicable 
                                # if MC_sampling parameter is set to 1.
    'calculate_dN_dxt':1,       # Output dN_dxt table. Only applicable 
                                # if MC_sampling parameter is set to 1.
    'output_dN_dxtdy_4all': 0,  # Output dN_dxtdy table. Only applicable 
                                # if MC_sampling parameter is set to 2.
}


# hadronic afterburner toolkit
hadronic_afterburner_toolkit_dict = {
    'echo_level': 9,    # control the mount of print messages
    'read_in_mode': 2,  # mode for reading in particle information
                        # 0: reads outputs from OSCAR outputs
                        # 1: reads outputs from UrQMD outputs
                        # 2: reads outputs from zipped UrQMD outputs
                        # 3: reads outputs from Sangwook's UrQMD outputs 
                        #    (without header lines)
                        # 4: reads outputs from UrQMD 3.3p2 outputs
                        # 10: reads outputf from gzip outputs
    'analyze_flow': 1,                  # 0/1: flag to perform flow analysis
    'analyze_HBT': 0,                   # 0/1: flag to perform HBT analysis
    'analyze_balance_function': 0,      # 0/1: flag to analyze Balance function
    'analyze_ebe_yield': 0,             # 0/1: flag to analyze ebe dis. of particle yield
    'read_in_real_mixed_events': 0,     # 0/1: read in real mixed events
    'randomSeed': -1,
    'particle_monval': 211,     # particle Monte-Carlo number
    'distinguish_isospin': 1,   # flag whether to distinguish the isospin of particles
    'event_buffer_size': 100000,       # the number of events read in at once
    'resonance_weak_feed_down_flag': 0,  # include weak feed down contribution
    'resonance_feed_down_flag': 0,  # perform resonance feed down
                                    # (will read in all hadrons and filter particle
                                    #  after decays are performed)
    'select_resonances_flag': 0,    # perform resonance decays only for selected particle species
    'resonance_weak_feed_down_Sigma_to_Lambda_flag': 0,    # include weak feed down contribution
                                                           # turn on only for Lambda (monval=3122)
                                                           # for Sigma^0 -> Lambda + gamma
    'net_particle_flag': 0,         # flag to collect net particle yield distribution
    # Parameters for single particle spectra and vn
    'rapidity_shift': 0.,
    'order_max': 9,     # the maximum harmonic order of anisotropic flow
    'compute_correlation': 0,       # flag to compute correlation function
    'flag_charge_dependence': 0,    # flag to compute charge dependence correlation
    'compute_corr_rap_dep': 0,      # flag to compute the rapidity dependent multi-particle correlation
    'npT': 41,          # number of pT points for pT-differential spectra and vn
    'pT_min': 0.05,     # the minimum value of transverse momentum (GeV)
    'pT_max': 4.05,     # the maximum value of transverse momentum (GeV)
    'rap_min': -0.5,    # minimum value of rapidity integration range for mid-rapidity observables 
    'rap_max': 0.5,     # maximum value of rapidity integration range for mid-rapidity observables 

    'rap_type': 1,      # 0: for pseudo-rapidity; 1: for rapidity
    'rapidity_distribution': 1,   # 1: output particle rapidity distribution 
    'n_rap': 141,                 # numpber of points in rapidity distr.
    'rapidity_dis_min': -7.0,     # minimum value of particle rapidity distribution
    'rapidity_dis_max': 7.0,      # maximum value of particle rapidity distribution
    'vn_rapidity_dis_pT_min': 0.20,  # the minimum value of pT for vn rap. distr.
    'vn_rapidity_dis_pT_max': 3.0,   # the maximum value of pT for vn rap. distr.

    'check_spatial_dis': 0,         # flag to check dN/dtau distribution
    'intrinsic_detas': 0.1,         # deta_s in the output samples
    'intrinsic_dtau': 0.01,         # dtau in the output samples
    'intrinsic_dx': 0.1,            # dx in the output samples
    # Parameters for HBT correlation functions
    'long_comoving_boost': 0,              # whether qlong will be boost by the pair velocity
    'needed_number_of_pairs': 30000000,    # number of pairs for eack K point
    'number_of_oversample_events': 100,    # number of the combined events in the numerator
    'number_of_mixed_events': 50,          # number of the mixed events in the denorminator
    'invariant_radius_flag': 0,     # 0: compute 3D HBT correlation function
                                    # 1: compute 1D HBT correlation function for q_inv
    'azimuthal_flag': 0,            # 0: compute the azimuthal averaged HBT correlation function
                                    # 1: compute azimuthal dependent HBT correlation function
    'kT_differenitial_flag': 1,     # 0: integrate the pair momentum k_T over  a given kT range for correlation function
                                    # 1: compute the correlation function at each specifiec kT point
    'n_KT': 6,      # number of the pair momentum k_T to calculate
    'KT_min': 0.0,  # minimum value of the pair momentum k_T 
    'KT_max': 1.0,  # maximum value of the pair momentum k_T 
    'n_Kphi': 48,   # number of the azimuthal angles for the pair momentum k_T 
                    # (range is assumed to be from 0 to 2*pi)
    'Krap_min': -0.5,   # minimum accept pair momentum rapidity
    'Krap_max': 0.5,    # maximum accept pair momentum rapidity
    'buffer_rapidity': 5.0,     # collect particles with rapidity from [Krap_min - buffer_rapidity, Krap_max + buffer_rapidity]
    'qnpts': 31,    # number of points for momentum q (difference of the pair momentum) for correlaction function
    'q_min': -0.15,     # minimum value for momentum q (GeV)
    'q_max': 0.15,      # maximum value for momentum q (GeV)
    'reject_decay_flag': 0,         # reject particles from resonance decays
                                    # 0: no rejection
                                    # 1: reject particles from all decay channels
                                    # 2: reject particles only from long lived resonance decays (future)
    'tau_reject': 10.,              # reject decay particle whose tau_f > tau_reject
                                    # only effective when reject_decay_flag == 2
    # options for calculting Balance function
    'particle_alpha': 9998,     # monte carlo number for particle alpha
    'particle_beta': -9998,     # monte carlo number for particle beta
    'Bnpts': 21,                # number of bins for the balance function
    'Brap_max': 2.0,    # the maximum \Delta y rapidity for balance function
    'BpT_min': 0.2,     # the minimum pT cut for particles used in balance function
    'BpT_max': 3.0,     # the maximum pT cut for particles used in balance function
}


Parameters_list = [
    (ipglasma_dict, "input", 3),
    (kompost_dict, "setup.ini", 4),
    (mcglauber_dict, "input", 0),
    (music_dict, "music_input_mode_2", 2),
    (photon_dict, "parameters.dat", 1),
    (iss_dict, "iSS_parameters.dat", 1),
    (hadronic_afterburner_toolkit_dict, "parameters.dat", 1)
]

path_list = [
    'model_parameters/IPGlasma/',
    'model_parameters/KoMPoST/',
    'model_parameters/3dMCGlauber/',
    'model_parameters/MUSIC/',
    'model_parameters/photonEmission_hydroInterface/',
    'model_parameters/iSS/',
    'model_parameters/hadronic_afterburner_toolkit/'
]


def update_parameters_dict(par_dict_path, ran_seed):
    """This function update the parameters dictionaries with user's settings"""
    par_diretory = path.dirname(par_dict_path)
    sys.path.insert(0, par_diretory)
    print(par_diretory)
    parameters_dict = __import__(par_dict_path.split('.py')[0].split('/')[-1])
    initial_condition_type = (
                    parameters_dict.control_dict['initial_state_type'])
    if initial_condition_type in ("IPGlasma", "IPGlasma+KoMPoST"):
        ipglasma_dict.update(parameters_dict.ipglasma_dict)

        # set random seed
        if ran_seed == -1:
            ipglasma_dict['useTimeForSeed'] = 1
        else:
            ipglasma_dict['seed'] = ran_seed
            ipglasma_dict['useTimeForSeed'] = 0

        if 'Initial_profile' not in parameters_dict.music_dict:
            parameters_dict.music_dict['Initial_profile'] = 9
        if 'Initial_Distribution_input_filename' not in parameters_dict.music_dict:
            parameters_dict.music_dict[
                'Initial_Distribution_input_filename'] = (
                        'initial/epsilon-u-Hydro.dat')
        if 'boost_invariant' not in parameters_dict.music_dict:
            parameters_dict.music_dict['boost_invariant'] = 1

        if 'Include_Rhob_Yes_1_No_0' not in parameters_dict.music_dict:
            parameters_dict.music_dict['Include_Rhob_Yes_1_No_0'] = 0
        if initial_condition_type == "IPGlasma+KoMPoST":
            # update KoMPoST parameters
            for subdict in parameters_dict.kompost_dict:
                kompost_dict[subdict].update(
                            parameters_dict.kompost_dict[subdict])
            parameters_dict.music_dict['s_factor'] = 1.0
            parameters_dict.music_dict['Initial_time_tau_0'] = (
                    kompost_dict['KoMPoSTInputs']['tOut'])
    else:
        mcglauber_dict.update(parameters_dict.mcglauber_dict)

        mcglauber_dict['seed'] = ran_seed       # set random seed

        if 'Initial_profile' not in parameters_dict.music_dict:
            parameters_dict.music_dict['Initial_profile'] = 13
        if 'Initial_Distribution_input_filename' not in parameters_dict.music_dict:
            parameters_dict.music_dict[
                'Initial_Distribution_input_filename'] = 'initial/strings.dat'
        if 'boost_invariant' not in parameters_dict.music_dict:
            parameters_dict.music_dict['boost_invariant'] = 0

        if 'Include_Rhob_Yes_1_No_0' not in parameters_dict.music_dict:
            parameters_dict.music_dict['Include_Rhob_Yes_1_No_0'] = 1

    if parameters_dict.music_dict['boost_invariant'] == 1:
        parameters_dict.iss_dict['hydro_mode'] = 1
    else:
        parameters_dict.iss_dict['hydro_mode'] = 2

    if 'calculate_polarization' in parameters_dict.iss_dict.keys():
        if parameters_dict.iss_dict['calculate_polarization'] == 1:
            parameters_dict.music_dict['output_vorticity'] = 1

    music_dict.update(parameters_dict.music_dict)
    if 'compute_photon_emission' in parameters_dict.control_dict:
        if parameters_dict.control_dict['compute_photon_emission']:
            music_dict['output_evolution_data'] = 2

    if 'photon_dict' in dir(parameters_dict):
        photon_dict.update(parameters_dict.photon_dict)
    iss_dict.update(parameters_dict.iss_dict)
    iss_dict['randomSeed'] = ran_seed
    hadronic_afterburner_toolkit_dict.update(
        parameters_dict.hadronic_afterburner_toolkit_dict)
    hadronic_afterburner_toolkit_dict['randomSeed'] = ran_seed


def update_parameters_bayesian(bayes_file):
    parfile = open(bayes_file, "r")
    for line in parfile:
        key, val = line.split()
        if key in mcglauber_dict.keys():
            mcglauber_dict[key] = float(val)
        if key in music_dict.keys():
            music_dict[key] = float(val)


def output_parameters_to_files(workfolder="."):
    """This function outputs parameters in dictionaries to files"""
    workfolder = path.abspath(workfolder)
    print("\U0001F375  Output input parameter files to {}...".format(
                                                                workfolder))
    for idict, (parameters_dict, fname, itype) in enumerate(Parameters_list):
        output_folder = path.join(workfolder, path_list[idict])
        if not path.exists(output_folder):
            makedirs(output_folder)
        f = open(path.join(output_folder, fname), "w")
        for key_name in parameters_dict:
            if itype in (0, 2):
                f.write("{parameter_name}  {parameter_value}\n".format(
                    parameter_name=key_name,
                    parameter_value=parameters_dict[key_name]))
            elif itype == 1:
                f.write("{parameter_name} = {parameter_value}\n".format(
                    parameter_name=key_name,
                    parameter_value=parameters_dict[key_name]))
            elif itype == 3:
                if key_name in ("type", "database_name_pattern"): continue
                f.write("{parameter_name}  {parameter_value}\n".format(
                    parameter_name=key_name,
                    parameter_value=parameters_dict[key_name]))
            elif itype == 4:
                f.write("[{}]\n".format(key_name))
                for subkey_name in parameters_dict[key_name]:
                    f.write("{parameter_name} = {parameter_value}\n".format(
                        parameter_name=subkey_name,
                        parameter_value=parameters_dict[key_name][subkey_name])
                    )
        if itype == 2:
            f.write("EndOfData")
        elif itype == 3:
            f.write("EndOfFile")
        f.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
            description='\U0000269B Welcome to iEBE-MUSIC parameter master',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-path', '--path', metavar='',
                        type=str, default='.',
                        help='output folder path')
    parser.add_argument('-par', '--par_dict', metavar='',
                        type=str, default='parameters_dict_user',
                        help='user-defined parameter dictionary filename')
    parser.add_argument('-b', '--bayes_file', metavar='',
                        type=str, default='',
                        help='parameters from bayesian analysis')
    parser.add_argument('-seed', '--random_seed', metavar='',
                        type=int, default=-1,
                        help='input random seed')
    args = parser.parse_args()
    update_parameters_dict(path.abspath(args.par_dict), args.random_seed)
    if args.bayes_file != "":
        update_parameters_bayesian(args.bayes_file)
    output_parameters_to_files(args.path)
