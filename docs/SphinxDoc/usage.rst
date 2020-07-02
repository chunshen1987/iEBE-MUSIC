Usage
=====

Setting up jobs
---------------

After all the code packages are compiled, the user can go back to the root
directory of iEBE-MUSIC and setup simulations with the script
:code:`generate_jobs.py`.

One can type :code:`./generate_jobs.py -h` for help information

The iEBE-MUSIC framework uses openMP and MPI to releaze parallelization at
multiple levels to simulte event-by-event simulations. For example, if a
user wants to generate 1024 minimum bias events for Au+Au collisions
at 200 GeV, the :code:`-n` or :code:`--n_jobs` option will divide the total
number of events into n seperate jobs. If we take n = 64, then each job
will run 16 collision events. This number 16 needs to be set to the parameter
:code:`-n_hydro` or :code:`--n_hydro_per_job`. With these two parameters, the
users can set up the framework to simulation 1024 events in total. The n = 64
jobs will be launched on the Cluster as 64 independent jobs and run in a
parallel fashion. For the NERSC Clusters, a MPI job wrapper is provided such
that these 64 jobs are launched inside the MPI with only one big job request
in the job submission system. In each job, the :code:`-n_hydro` events will
be run in a sequential manner. If there are plenty of CPU available on the
cluster, for example on Open Science Grid, we recommend to set
:code:`-n_hydro` to 1 to maximally parallize in the collision event direction.

In each job, the user can use openMP to speed up the 3D hydrodynamic
simulations. The number of openMP threads can be set using the parameter
:code:`-n_th` or :code:`--n_threads`. Because hadronic transport code only
supports single thread, we parallize :code:`-n_urqmd` transport simulation
at the framework level. Users can set :code:`n_urqmd` = :code:`n_th` to use
all the available resource available after hydrodynamic simualtions.

After setting up jobs, one can use the script :code:`submit_all_jobs.sh` to
submit all the jobs to cluster. On NERSC, the job submission script will be
generated at the work_folder. One can go to that directory and submit the job
using slurm command.


Collecting results after simulations
------------------------------------

After the all simulations finish, the framework will save the
particle spectra, flow Q_n vectors, multiparticle correlations (optional)
in a hdf5 file for every event. These results are analyzed by the code
package, :code:`hadronic_afterburner_toolkit`. Please see the code readme
for details about the analysis. Then the users can use a script
:code:`collect_events.sh` to combine all the hdf5 file from inidividual
event into one database. With the API provided by h5py, the users can
easily access the results for every event.

To perform event averaging, one can use the provided python script
:code:`average_event_spvn_h5.py`. This script will output the final event
averaged results in ascii format for users to make plots. If one needs to
perform centrality selection and event averging, there is an integrated
python script :code:`average_event_spvn_h5_minimumbias.py` provided by
the :code:`hadronic_afterburner_toolkit` package. The user can find this
script under :code:`codes/hadronic_afterburner_toolkit_code/ebe_scripts/`

Additional information can be saved by setting options in the user's
parameter_dict file. In the :code:`control_dict` the options,

- :code:`save_ipglasma_results`

  This boolean decides to save the initial condition from IPGlasma

- :code:`save_kompost_results`

  This boolean decides whether to save the T^{\mu\nu} tensor at the
  starting time of hydrodynamics for the IPGlasma+KoMPoST running mode

- :code:`save_hydro_surfaces`

  This boolean decides to save the particlization hyper-surface

- :code:`save_UrQMD_files`

  This boolean decides to save the final state particle list after the UrQMD
  finishes


Data generation for Bayesian Analysis
-------------------------------------

In additional to the user defined parameter_dict file, he can provide an
additional parameter file to change a small set of parameters that will be
used for Bayesian analysis. This is acheived by pass the parameter file
with :code:`-b` or :code:`--bayes_file` option in the :code:`generate_jobs.py`.

