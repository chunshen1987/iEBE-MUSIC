#!/usr/bin/env python
"""
    This script contains all the user modified parameters in
    the iEBE-MUSIC package.
"""

# initial condition
initial_dict = {
    'initial_state_type': "3DMCGlauber",  # 3DMCGlauber, IPGlasma
}


# IPGlasma
ipglasma = {
    'database_name': "IPGlasma_database/AuAu_C0-5.h5",  # path for the database file
}


# 3DMCGlauber model
mcglauber_dict = {
    'database_name': "self",     # self: generate initial condition on the fly
    'Projectile':  "Pb",         # projectile nucleus name
    'Target'    :  "Pb",         # target nucleus name
    'roots'     :   17.3,        # collision energy (GeV)
    'seed'      :   -1,          # random seed (-1: system)
    'baryon_junctions': 1,       # 0: baryon number assumed to be at string end
                                 # 1: baryon number transported assuming baryon
                                 # junctions (at smaller x)
                                 # see arXiv:nucl-th/9602027
    'lambdaB': 0.4,              # parameter the controls the strength of
                                 # the baryon junction stopping
    'evolve_QCD_string_mode': 2,        # string evolution mode
                                        # 1: deceleration with fixed rapidity loss (m/sigma = 1 fm, dtau = 0.5 fm)
                                        # 2: deceleration with LEXUS sampled rapidit loss (both dtau and sigma fluctuate)
                                        # 3: deceleration with LEXUS sampled rapidit loss (m/sigma = 1 fm, dtau fluctuates)
                                        # 4: deceleration with LEXUS sampled rapidit loss (dtau = 0.5 fm, m/sigma fluctuates)
}


# MUSIC
music_dict = {
    'Initial_profile': 13,  # type of initial condition 
                            # 9: IPGlasma (full Tmunu),
                            #   -- 91: e and u^\mu,
                            #   -- 92: e only,
                            #   -- 93: e, u^\mu, and pi^\munu
                            # 13: dynamical initialization (3dMCGlauber)
    # read in initial conditions from external file
    'Initial_Distribution_input_filename': 'initial/strings.dat',
    's_factor': 1.400,      # normalization factor read in initial data file
    'Initial_time_tau_0': 0.4,  # starting time of the hydrodynamic evolution (fm/c)
    'Delta_Tau': 0.010,         # time step to use in the evolution [fm/c]
    'boost_invariant':  0,  # whether the simulation is boost-invariant 
    'Eta_grid_size': 12.0,  # spatial rapidity range 
                            # [-Eta_grid_size/2, Eta_grid_size/2 - delta_eta]
    'Grid_size_in_eta': 64,  # number of the grid points in spatial rapidity direction
    'X_grid_size_in_fm': 24.0,  # spatial range along x direction in the transverse plane 
                                # [-X_grid_size_in_fm/2, X_grid_size_in_fm/2]
    'Y_grid_size_in_fm': 24.0,  # spatial range along x direction in the transverse plane 
                                # [-X_grid_size_in_fm/2, X_grid_size_in_fm/2]
    'Grid_size_in_x': 256,      # number of the grid points in x direction
    'Grid_size_in_y': 256,      # number of the grid points in y direction
    'EOS_to_use': 14,           # type of the equation of state
                                # 0: ideal gas
                                # 1: EOS-Q from azhydro
                                # 7: lattice EOS s95p-v1.2 for UrQMD
                                # 9: hotQCD EOS with UrQMD
                                # 14: neos_BQS lattice EoS at finite mu_B
                                # 17: BEST lattice EoS at finite mu_B
    # transport coefficients
    'quest_revert_strength': 1.0,          # the strength of the viscous regulation
    'Viscosity_Flag_Yes_1_No_0': 1,        # turn on viscosity in the evolution
    'Include_Shear_Visc_Yes_1_No_0': 1,    # include shear viscous effect
    'Shear_to_S_ratio': 0.08,              # value of \eta/s
    'T_dependent_Shear_to_S_ratio': 0,     # flag to use temperature dep. \eta/s(T)
    'Include_Bulk_Visc_Yes_1_No_0': 0,     # include bulk viscous effect
    'T_dependent_zeta_over_s': 7,          # parameterization of \zeta/s(T)
    'Include_second_order_terms': 1,       # include second order non-linear coupling terms
    'Include_vorticity_terms': 0,          # include vorticity coupling terms
    'Include_Rhob_Yes_1_No_0': 0,
    'turn_on_baryon_diffusion': 0,
    'kappa_coefficient': 0.4,

    # parameters for freeze out and Cooper-Frye 
    'N_freeze_out': 1,
    'eps_freeze_max': 0.26,
    'eps_freeze_min': 0.26,
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
}


# hadronic afterburner toolkit
hadronic_afterburner_toolkit_dict = {
    'event_buffer_size': 100,       # the number of events read in at once
    'compute_correlation': 0,       # flag to compute correlation function
    'flag_charge_dependence': 0,    # flag to compute charge dependence correlation
}
