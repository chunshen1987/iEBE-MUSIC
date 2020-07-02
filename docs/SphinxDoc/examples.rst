Examples
========

There are example parameter files under the folder :code:`config` for
different running modes of the iEBE-MUSIC framework.

3DMCGlauber_consttau
--------------------

In this mode, the iEBE-MUSIC will perform 3D hydrodynamics + hadronic
casade simulations with given initial condition at a fixed proper time,
:math:`\tau_0`.
For the current settings, users needs to provide the nuclear thickness
function :math:`T_A(x, y)` and :math:`T_B(x, y)`.
A 3D initial condition for energy
density and net bayron density profiles at a fixed proper time tau will
be constructed based on the given :math:`T_A` and :math:`T_B`.
Please see arXiv: 2003.05852 [nucl-th] for a detailed physics model.

To run this mode, users can start with the example parameter file
:code:`config/parameters_dict_user_3DMCGlauber_consttau.py`.
You need to set :code:`initial_state_type` to "3DMCGlauber_consttau"
in the :code:`control_dict`. Then the user needs to provide the absolute
directory path where the nuclear thickeness functions are stored. This path
need to be set to :code:`database_name` in the :code:`mcglauber_dict`. The
script to generate jobs will run all the files with names
"nuclear_thickness_TA_fromSd_order_2_C*.dat" and
"nuclear_thickness_TB_fromSd_order_2_C*.dat".
If users want to use a different naming convension, they can modify the
:code:`generate_jobs.py` accordingly.

All the hydrodynamic parameters are specified in the :code:`music_dict`
dictionary. The mode "3DMCGlauber_consttau" has it unique parameters

- :code:`Eta_plateau_size`
- :code:`Eta_fall_off`
- :code:`eta_rhob_0`
- :code:`eta_rhob_width_1`
- :code:`eta_rhob_width_2`

These paramters specify the longitudinal profile of energy density and net
baryon density. The proper values of those parameters as different collsion
energies can be found in arXiv: 2003.05852 [nucl-th].


3DMCGlauber_dynamical
---------------------


IPGlasma
--------


IPGlasma + KoMPoST
------------------
