# iEBE-MUSIC
This is a repository that stores all the scripts for running event-by-event relativistic heavy-ion simulations.

## Ingradients:
All the code packages can be downloaded from online git repositories. Please use `codes/get_code_packages.sh` to download the code packages and `codes/compile_code_packages.sh` to compile all the packages before event-by-event simulations.

Initial conditions:

- IPGlasma
- 3DMCGlauber: a 3D Monte-Carlo Glauber model for heavy-ion collisions

Hydrodynamics:

- MUSIC

Particlization & hadronic transport:

- iSS + UrQMD

## Usage:

./generate\_jobs.py initial\_condition\_type initial\_condition\_filename working\_folder cluster\_name n\_jobs n\_hydro\_per\_job n\_UrQMD\_per\_hydro [n\_threads]

cluster_name:

- nersc
- wsugrid
- local
- guillimin
- McGill

initial\_condition\_type:

- IPGlasma
- 3DMCGlauber

1. If initial\_condition\_filename == "self", initial condition will be generated on the fly with the rest part of the simulations.

2. If initial\_condition\_filename == database_file, the code package will use the pre-generated initial coniditions stored in the database (HDF5).


## Settings on NERSC:

The clusters on NERSC does not use utf-8 as default. To run the script properly, one needs to add the following commands in the ~/.bashrc.ext file,

```
export PYTHONIOENCODING=utf-8
export LC_CTYPE=en_US.UTF8
```



