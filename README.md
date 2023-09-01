# iEBE-MUSIC
This is a repository is an overarching numerical framework for event-by-event simulations of relativistic heavy-ion collisions.

If you have any questions, please email to the iEBE-MUSIC google groups, iebe-music@googlegroups.com


## Setup & Ingradients:
All the code packages can be downloaded from online git repositories. Please use `codes/get_code_packages.sh` to download the code packages and `codes/compile_code_packages.sh` to compile all the packages before event-by-event simulations.

Initial conditions:

- [IPGlasma](https://github.com/schenke/ipglasma)
- [3DMCGlauber](https://github.com/chunshen1987/3dMCGlauber): a 3D Monte-Carlo Glauber model for heavy-ion collisions

Pre-equilibrium evolution:

- [KoMPoST](https://github.com/KMPST/KoMPoST)

Hydrodynamics:

- [MUSIC](https://github.com/MUSIC-fluid/MUSIC)

Particlization & hadronic transport:

- [iSS](https://github.com/chunshen1987/iSS) + [UrQMD](https://Chunshen1987@bitbucket.org/Chunshen1987/urqmd_afterburner.git)


## Usage:

type `./generate_jobs.py -h` for help information

## Parameters:
Users can pass model parameters through a python script `parameters_dict_user.py`. It contains multiple dictionaries, which are related to each code module inside the iEBE-MUSIC framework. This script will update the master parameters dictionaries in `config/parameters_dict_master.py`. One can read `config/parameters_dict_master.py` for all the available parameters options for each module. If a user want to modify any parameters, he can add it in the `parameters_dict_user.py`.

#### initial_dict
The `initial_dict` dictionary specify the type of the initial condition to use,

- IPGlasma
- 3DMCGlauber

For the IPGlasma initial condition, the user needs to specify the HDF5 database filename. 

For the 3DMCGlauber initial condition,

1. If `database_name : "self"`, initial condition will be generated on the fly with the rest part of the simulations.

2. If `database_name : database_file`, the code package will use the pre-generated initial coniditions stored in the database (HDF5).

#### music_dict
The `music_dict` has all the parameters to run MUSIC for (3+1)D hydrodynamic simulations. 

    'Initial_profile': 9,  for IPGlasma initial condition
        - 9: IPGlasma (full Tmunu),
        - 91: e and u^\mu,
        - 92: e only,
        - 93: e, u^\mu, and pi^\munu
    'Initial_profile': 13,  for 3D MCGlauber initial condition
        - 13: 3D MCGlauber initial condition with dynamical initialization,
        - 131: 3D MCGlauber initial condition with instantaneous initialization

## Settings on NERSC:

The clusters on NERSC does not use utf-8 as default. To run the script properly, one needs to add the following commands in the ~/.bashrc.ext file,

```
export PYTHONIOENCODING=utf-8
export LC_CTYPE=en_US.UTF8
```

## Docker Support

The iEBE-MUSIC has its official docker image on docker hub [iebe-music](https://hub.docker.com/r/chunshen1987/iebe-music).

## Coding Style

We use YAPF to impose coding format for the python scripts. Before every commit, please use

    yapf -i filename.py

to apply the uniformed format to the source code files
