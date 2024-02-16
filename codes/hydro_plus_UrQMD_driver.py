#!/usr/bin/env python3
"""This is a drive script to run hydro + hadronic cascade simulation"""

from multiprocessing import Pool
from subprocess import call
from os import path, mkdir, remove, makedirs, stat
import tarfile
from glob import glob
import sys
import time
import shutil
import re
import h5py
import numpy as np
from fetch_IPGlasma_event_from_hdf5_database import fecth_an_IPGlasma_event, fecth_an_IPGlasma_event_Tmunu
from fetch_3DMCGlauber_event_from_hdf5_database import fecth_an_3DMCGlauber_event


def print_usage():
    """This function prints out help messages"""
    print("\U0001F3B6  " + "Usage: {} ".format(sys.argv[0])
          + "initial_condition_type initial_condition_database "
          + "n_hydro_events hydro_event_id n_UrQMD n_threads "
          + "save_ipglasma_flag save_kompost_flag save_hydro_flag "
          + "save_urqmd_flag seed_add tau0 compute_polarization_flag "
          + "compute_photons_flag enableCheckPoint")


def fecth_an_3DMCGlauber_smooth_event(database_path, iev):
    """This function returns the filename of an initial condition in the
       database_path folder
    """
    filelist = glob(path.join(database_path, 'nuclear_thickness_TA_*.dat'))
    return (filelist[iev])


def mapEventIdToCentrality(event_id):
    """This function maps the event_id to centrality bin for 3D-Glauber model
       It makes sure the centrality are evenly distributed over different
       simulation jobs when nev > 1
    """
    event_id = event_id % 100
    b = event_id % 10
    a = int((event_id - b)/10)
    centrality = 0
    if (b % 2) == 0:
        centrality = b/2*10 + a
    else:
        centrality = (9 - (b - 1)/2)*10 + a
    return(centrality)


def get_initial_condition(database, initial_type, iev, event_id, seed_add,
                          final_results_folder, time_stamp_str="0.4"):
    """This funciton get initial conditions"""
    status = True
    if "IPGlasma" in initial_type:
        ipglasma_local_folder = "ipglasma/ipglasma_results"
        res_path = path.join(path.abspath(final_results_folder),
                             "ipglasma_results_{}".format(event_id))
        file_name = ("epsilon-u-Hydro-t{0:s}-{1}.dat".format(
                                                time_stamp_str, event_id))
        if "KoMPoST" in initial_type:
            file_name = ("Tmunu-t{0:s}-{1}.dat".format(time_stamp_str,
                                                       event_id))
        if database == "self":
            # check existing events ...
            if not path.exists(path.join(res_path, file_name)):
                run_ipglasma(event_id)
                collect_ipglasma_event(res_path)
                if not path.exists(path.join(res_path, file_name)):
                    # IPGlasma event failed
                    print("IPGlasma event failed ... ")
                    status = False
            else:
                print("IPGlasma event exists ...")
                print("No need to rerun ...")
        else:
            if "KoMPoST" in initial_type:
                file_temp = fecth_an_IPGlasma_event_Tmunu(
                                            database, time_stamp_str, event_id)
            else:
                file_temp = fecth_an_IPGlasma_event(database, time_stamp_str,
                                                    event_id)
            makedirs(ipglasma_local_folder, exist_ok=True)
            shutil.move(file_temp,
                        path.join(ipglasma_local_folder, file_name))
            collect_ipglasma_event(res_path)
        connect_ipglasma_event(res_path, initial_type, file_name)
        return status, file_name
    elif initial_type == "3DMCGlauber_dynamical":
        if database == "self" or "fixCentrality":
            file_name = "strings_event_{}.dat".format(event_id)
            specFilename = "spectators_event_{}.dat".format(event_id)
            ran = np.random.default_rng().integers(1e8)
            if not path.exists(file_name):
                if database == "self":
                    cenMin = mapEventIdToCentrality(event_id)
                    call("(cd 3dMCGlauber; ./3dMCGlb.e 1 input "
                         + "{} cenMin={} cenMax={};)".format(seed_add+iev*ran,
                                                             cenMin, cenMin+1),
                         shell=True)
                else:
                    call("(cd 3dMCGlauber; ./3dMCGlb.e 1 input "
                         + "{};)".format(seed_add+iev*ran), shell=True)
                call("mv 3dMCGlauber/strings_event_0.dat {}".format(file_name),
                     shell=True)
                call("mv 3dMCGlauber/spectators_event_0.dat {}".format(
                                                specFilename), shell=True)
            else:
                print("3D MC-Glauber event exists ...")
                print("No need to rerun ...")
            shutil.copy(file_name, "MUSIC/initial/strings.dat")
            shutil.copy(file_name, path.join(final_results_folder,
                                             "strings_{}.dat".format(event_id)))
            shutil.copy(specFilename,
                        path.join(final_results_folder,
                                  "spectators_{}.dat".format(event_id)))
            filePatterns = ["ed_etas", "nB_etas", "ecc_ed"]
            call("mv 3dMCGlauber/ed_etas_*.dat {}".format(
                                    final_results_folder), shell=True)
            call("mv 3dMCGlauber/nB_etas_*.dat {}".format(
                                    final_results_folder), shell=True)
            call("mv 3dMCGlauber/nQ_etas_*.dat {}".format(
                                    final_results_folder), shell=True)
            call("mv 3dMCGlauber/ecc_ed*.dat {}".format(
                                    final_results_folder), shell=True)
            return status, file_name
        else:
            file_name = fecth_an_3DMCGlauber_event(database, event_id)
            return status, file_name
    elif initial_type == "3DMCGlauber_participants":
        file_name = "participants_event_{}.dat".format(event_id)
        specFilename = "spectators_event_{}.dat".format(event_id)
        ran = np.random.default_rng().integers(1e8)
        if not path.exists(file_name):
            if database == "self":
                cenMin = mapEventIdToCentrality(event_id)
                call("(cd 3dMCGlauber; ./3dMCGlb.e 1 input "
                     + "{} cenMin={} cenMax={};)".format(seed_add + iev*ran,
                                                         cenMin, cenMin+1),
                     shell=True)
            else:
                call("(cd 3dMCGlauber; ./3dMCGlb.e 1 input "
                     + "{};)".format(seed_add + iev*ran), shell=True)
            call("mv 3dMCGlauber/participants_event_0.dat {}".format(file_name),
                 shell=True)
            call("mv 3dMCGlauber/spectators_event_0.dat {}".format(
                                            specFilename), shell=True)
        else:
            print("3D MC-Glauber event exists ...")
            print("No need to rerun ...")
        shutil.copy(file_name, "MUSIC/initial/participants_event.dat")
        shutil.copy(file_name, path.join(final_results_folder, file_name))
        shutil.copy(specFilename,
                    path.join(final_results_folder,
                              "spectators_{}.dat".format(event_id)))
        return status, file_name
    elif initial_type == "3DMCGlauber_consttau":
        file_name = fecth_an_3DMCGlauber_smooth_event(database, event_id)
        return status, file_name
    else:
        print("\U0001F6AB  "
              + "Do not recognize the initial condition type: {}".format(
                  initial_type))
        sys.exit(1)


def run_ipglasma(event_id):
    """This functions run IPGlasma"""
    print("\U0001F3B6  Run IPGlasma ... ")
    call("bash ./run_ipglasma.sh {}".format(event_id), shell=True)


def collect_ipglasma_event(final_results_folder):
    """This function collects the ipglasma results"""
    if path.exists(final_results_folder):
        shutil.rmtree(final_results_folder)
    shutil.move("ipglasma/ipglasma_results", final_results_folder)


def connect_ipglasma_event(res_path, initial_type, filename):
    if initial_type == "IPGlasma":
        hydro_initial_file = "MUSIC/initial/epsilon-u-Hydro.dat"
        if path.islink(hydro_initial_file):
            remove(hydro_initial_file)
        call("ln -s {0:s} {1:s}".format(path.join(res_path, filename),
                                        hydro_initial_file),
             shell=True)
    elif initial_type == "IPGlasma+KoMPoST":
        kompost_initial_file = "kompost/Tmunu.dat"
        if path.islink(kompost_initial_file):
            remove(kompost_initial_file)
        call("ln -s {0:s} {1:s}".format(path.join(res_path, filename),
                                        kompost_initial_file),
             shell=True)


def run_hydro_event(final_results_folder, event_id):
    """This functions run hydro"""
    logo = "\U0001F3B6"
    hydro_folder_name = "hydro_results_{}".format(event_id)
    results_folder = path.join(final_results_folder, hydro_folder_name)
    hydro_success = False

    if path.exists(results_folder):
        print("{}  Hydrodynaimc results {} exist ... ".format(
            logo, hydro_folder_name),
              flush=True)
        # check hydro finishes properly
        try:
            ftmp = open(path.join(results_folder, "run.log"),
                        'r',
                        encoding="utf-8")
            hydro_status = ftmp.readlines()[-1].split()[3]
            ftmp.close()
            if hydro_status == "Finished.":
                print("{} Hydrodynamic run finished properly ... ".format(logo),
                      flush=True)
                hydro_success = True
        except FileNotFoundError:
            hydro_success = False

        if not hydro_success:
            print("{} Hydrodynamic run failed, rerun ... ".format(logo),
                 flush=True)
            shutil.rmtree(results_folder)

    if not hydro_success:
        curr_time = time.asctime()
        print("{}  [{}] Playing MUSIC ... ".format(logo, curr_time), flush=True)
        call("bash ./run_hydro.sh", shell=True)

        # check hydro finishes properly
        ftmp = open("MUSIC/hydro_results/run.log", 'r', encoding="utf-8")
        hydro_status = ftmp.readlines()[-1].split()[3]
        if hydro_status == "Finished.":
            hydro_success = True

        # collect hydro results
        shutil.move("MUSIC/hydro_results", results_folder)

    return (hydro_success, hydro_folder_name)


def run_kompost(final_results_folder, event_id):
    """This functions run KoMPoST simulation"""
    logo = "\U0001F3B6"
    kompost_folder_name = "kompost_results_{}".format(event_id)
    results_folder = path.join(final_results_folder, kompost_folder_name)
    kompost_success = False

    if path.exists(results_folder):
        # check whether KoMPoST has already run or not
        print("{} KoMPoST results {} exist ...".format(logo,
                                                       kompost_folder_name),
              flush=True)
        kompost_success = True
        if kompost_success:
            print("{} no need to rerun KoMPoST".format(logo), flush=True)
        else:
            print("{} KoMPoST simulation failed, rerun ...".format(logo),
                  flush=True)
            shutil.rmtree(results_folder)

    if not kompost_success:
        curr_time = time.asctime()
        print("\U0001F3B6  [{}] Run KoMPoST ... ".format(curr_time), flush=True)
        call("bash ./run_kompost.sh", shell=True)

        kompost_success = True
        if kompost_success:
            # collect results
            shutil.move("kompost/kompost_results", results_folder)

    return (kompost_success, kompost_folder_name)


def prepare_evolution_files_for_photon(final_results_folder, hydro_folder_name):
    """This function prepares hydro evolution file for photon radiation"""
    evoFileName = "evolution_all_xyeta.dat"
    evo_file = path.join(final_results_folder, hydro_folder_name,
                         evoFileName)
    photonFolderPath = path.join('photonEmission_hydroInterface', 'results')
    if path.exists(photonFolderPath):
        shutil.rmtree(photonFolderPath)
    mkdir(photonFolderPath)
    call("ln -s {0:s} {1:s}".format(path.abspath(evo_file),
                                    path.join(photonFolderPath, evoFileName)),
         shell=True)
    shutil.copy(path.join(final_results_folder, hydro_folder_name,
                          "music_input"),
                photonFolderPath)


def run_photon(final_results_folder, event_id):
    """This functions run photon radiation"""
    logo = "\U0001F3B6"
    photon_folder_name = "photon_results_{}".format(event_id)
    results_folder = path.join(final_results_folder, photon_folder_name)
    photon_success = False

    if path.exists(results_folder):
        # check whether photon has already run or not
        print("{} photon results {} exist ...".format(logo,
                                                      photon_folder_name),
              flush=True)
        photon_success = True
        if photon_success:
            print("{} no need to rerun photon".format(logo), flush=True)
        else:
            print("{} photon radiation failed, rerun ...".format(logo),
                  flush=True)
            shutil.rmtree(results_folder)

    if not photon_success:
        curr_time = time.asctime()
        print("\U0001F3B6  [{}] Run photon ... ".format(curr_time), flush=True)
        call("bash ./run_photon.sh", shell=True)

        photon_success = True
        if photon_success:
            # collect results
            shutil.move("photonEmission_hydroInterface/results",
                        results_folder)

    return (photon_success, photon_folder_name)


def prepare_surface_files_for_urqmd(final_results_folder, hydro_folder_name,
                                    n_urqmd):
    """This function prepares hydro surface for hadronic casade"""
    surface_file = glob(
        path.join(final_results_folder, hydro_folder_name, "surface*.dat"))
    spectatorFileList = glob(path.join(final_results_folder, "spectator*.dat"))
    spectatorFile = ""
    if spectatorFileList != []:
        spectatorFile = spectatorFileList[0]
    if stat(surface_file[0]).st_size == 0:
        return False
    for iev in range(n_urqmd):
        hydro_surface_folder = "UrQMDev_{0:d}/hydro_event".format(iev)
        if path.exists(hydro_surface_folder):
            shutil.rmtree(hydro_surface_folder)
        mkdir(hydro_surface_folder)
        call("ln -s {0:s} {1:s}".format(
            path.abspath(surface_file[0]),
            path.join(hydro_surface_folder, "surface.dat")),
             shell=True)
        shutil.copy(
            path.join(final_results_folder, hydro_folder_name, "music_input"),
            hydro_surface_folder)
        if spectatorFile != "":
            shutil.copy(spectatorFile,
                        path.join(hydro_surface_folder, "spectators.dat"))
    return True


def run_urqmd_event(event_id):
    """This function runs hadornic afterburner"""
    call("bash ./run_afterburner.sh {0:d}".format(event_id), shell=True)


def run_urqmd_shell(n_urqmd, final_results_folder, event_id, para_dict,
                    startTime, checkPointFileName):
    """This function runs urqmd events in parallel"""
    logo = "\U0001F5FF"
    urqmdResults = "particle_list_{}.bin".format(event_id)
    results_folder = path.join(final_results_folder, urqmdResults)
    urqmd_success = False

    if path.exists(results_folder):
        print("{} UrQMD results {} exist ... ".format(logo, urqmdResults),
              flush=True)
        urqmd_success = True

    if not urqmd_success:
        curr_time = time.asctime()
        if para_dict["compute_polarization"]:
            spin_folder_name = "spin_results_{}".format(event_id)
            spin_folder = path.join(final_results_folder, spin_folder_name)
            if path.exists(spin_folder):
                print("{} spin results {} exist ... ".format(logo,
                                                             spin_folder),
                      flush=True)
            else:
                print("{}  [{}] Running spin calculations ... ".format(logo,
                                                                   curr_time),
                      flush=True)
                call("bash ./run_spinPol.sh {}".format(n_urqmd), shell=True)
                shutil.move("UrQMDev_{}/iSS/results".format(n_urqmd),
                            spin_folder)
                if para_dict["check_point_flag"]:
                    checkPoint(startTime, checkPointFileName,
                               final_results_folder)

        print("{}  [{}] Running UrQMD ... ".format(logo, curr_time),
              flush=True)
        with Pool(processes=n_urqmd) as pool1:
            pool1.map(run_urqmd_event, range(n_urqmd))

        urqmdResFile = "particle_list.bin"
        for iev in range(1, n_urqmd):
            call("cat UrQMDev_{}/UrQMD_results/{} ".format(iev, urqmdResFile)
                 + ">> UrQMDev_0/UrQMD_results/{}".format(urqmdResFile),
                 shell=True)
            remove("UrQMDev_{}/UrQMD_results/{}".format(iev, urqmdResFile))
        urqmd_success = True
        shutil.move("UrQMDev_0/UrQMD_results/{}".format(urqmdResFile),
                    results_folder)

    return (urqmd_success, results_folder)


def run_spvn_analysis(urqmd_file_path, n_threads, final_results_folder,
                      event_id):
    """This function runs analysis"""
    final_results_folder = path.join(final_results_folder,
                                     "spvn_results_{0:s}".format(event_id))
    if path.exists(final_results_folder):
        shutil.rmtree(final_results_folder)
    spvn_folder = "hadronic_afterburner_toolkit/results"
    if path.exists(spvn_folder):
        shutil.rmtree(spvn_folder)
    mkdir(spvn_folder)
    call("ln -s {0:s} {1:s}".format(path.abspath(urqmd_file_path),
                                    path.join(spvn_folder,
                                              "particle_list.bin")),
         shell=True)
    # finally collect results
    curr_time = time.asctime()
    print(f"\U0001F3CD  [{curr_time}] Running spvn analysis ... ", flush=True)

    call("bash ./run_analysis_spvn.sh", shell=True)

    curr_time = time.asctime()
    print(f"\U0001F3CD  [{curr_time}] Finished spvn analysis ... ", flush=True)

    call(f"rm {spvn_folder}/particle_list.bin", shell=True)
    shutil.move(spvn_folder, final_results_folder)


def check_an_event_is_good(event_folder):
    """This function checks the given event contains all required files"""
    required_files_list = [
        'particle_9999_vndata_eta_-0.5_0.5.dat',
        'particle_211_vndata_diff_y_-0.5_0.5.dat',
        'particle_321_vndata_diff_y_-0.5_0.5.dat',
        'particle_2212_vndata_diff_y_-0.5_0.5.dat',
    ]
    event_file_list = glob(path.join(event_folder, "*"))
    for ifile in required_files_list:
        filename = path.join(event_folder, ifile)
        if filename not in event_file_list:
            print("event {} is bad, missing {} ...".format(
                event_folder, filename),
                  flush=True)
            return False
    return True


def zip_results_into_hdf5(final_results_folder, event_id, para_dict):
    """This function combines all the results into hdf5"""
    results_name = "spvn_results_{}".format(event_id)
    time_stamp = para_dict['time_stamp_str']
    initial_state_filelist1 = [
        'NcollList{}.dat'.format(event_id),
        'NpartList{}.dat'.format(event_id),
        'NpartdNdy-t0.6-{}.dat'.format(event_id),
        'NgluonEstimators{}.dat'.format(event_id),
        'usedParameters{}.dat'.format(event_id),
    ]
    initial_state_filelist2 = [
        'epsilon-u-Hydro-t0.1-{}.dat'.format(event_id),
        'epsilon-u-Hydro-t{0}-{1}.dat'.format(time_stamp, event_id),
    ]
    glauber_filelist = ["strings_{}.dat".format(event_id),
                        "spectators_{}.dat".format(event_id),
                        "participants_event_{}.dat".format(event_id),
                        "ed_etas_distribution_*.dat",
                        "nB_etas_distribution_*.dat",
                        "nQ_etas_distribution_*.dat",
                        "ecc_ed_*.dat",]

    pre_equilibrium_filelist = [
        'ekt_tIn01_tOut08.music_init_flowNonLinear_pimunuTransverse.txt'
    ]
    hydro_info_filepattern = [
        "eccentricities_evo_*.dat", "momentum_anisotropy_*.dat",
        "meanpT_estimators_*.dat",
        "inverse_Reynolds_number_eta_*.dat",
        "averaged_phase_diagram_trajectory_eta_*.dat",
        "global_conservation_laws.dat", "global_angular_momentum_*.dat",
        "vorticity_evo_*.dat", "FO_nBvseta.dat",
    ]
    photon_filepattern = ['*_Spvn*.dat']
    spin_filepattern = [
        "Smu_dpTdphi_*.dat", "Smu_phi_*.dat", "Smu_pT_*.dat", "Smu_y_*.dat",
        "Smu_Thermal_*.dat", "Rspin_*.dat"
    ]

    hydrofolder = path.join(final_results_folder,
                            "hydro_results_{}".format(event_id))
    spvnfolder = path.join(final_results_folder, results_name)
    photonFolder = path.join(final_results_folder,
                             "photon_results_{}".format(event_id))

    spinfolder = path.join(final_results_folder,
                           "spin_results_{}".format(event_id))

    status = check_an_event_is_good(spvnfolder)
    if status:
        curr_time = time.asctime()
        print("[{}] {} is good, converting results to hdf5".format(
            curr_time, spvnfolder),
              flush=True)

        if para_dict['initial_condition'] == "self" or "fixCentrality":
            # save initial conditions
            if "IPGlasma" in para_dict['initial_type']:
                initial_folder = path.join(
                    final_results_folder,
                    "ipglasma_results_{}".format(event_id))
                for inifilename in initial_state_filelist1:
                    inifile = path.join(initial_folder, inifilename)
                    if path.isfile(inifile):
                        shutil.move(inifile, spvnfolder)
                if para_dict['save_ipglasma']:
                    for inifilename in initial_state_filelist2:
                        inifile = path.join(initial_folder, inifilename)
                        if path.isfile(inifile):
                            shutil.move(inifile, spvnfolder)
            if "3DMCGlauber" in para_dict['initial_type']:
                for iniFilenamePattern in glauber_filelist:
                    iniFileList = glob(path.join(final_results_folder,
                                                 iniFilenamePattern))
                    for iniFile in iniFileList:
                        if path.isfile(iniFile):
                            shutil.move(iniFile, spvnfolder)

            # save pre-equilibrium results
            if (para_dict['initial_type'] == "IPGlasma+KoMPoST"
                    and para_dict['save_kompost']):
                preeq_folder = path.join(final_results_folder,
                                         "kompost_results_{}".format(event_id))
                for prefilename in pre_equilibrium_filelist:
                    prefile = path.join(preeq_folder, prefilename)
                    if path.isfile(prefile):
                        shutil.move(prefile, spvnfolder)

        # save hydro evolution information
        for ipattern in hydro_info_filepattern:
            hydro_info_list = glob(path.join(hydrofolder, ipattern))
            for ihydrofile in hydro_info_list:
                if path.isfile(ihydrofile):
                    shutil.move(ihydrofile, spvnfolder)

        # save photon results
        if para_dict['compute_photons']:
            for ipattern in photon_filepattern:
                photonFileList = glob(path.join(photonFolder, ipattern))
                for photonFile_i in photonFileList:
                    if path.isfile(photonFile_i):
                        shutil.move(photonFile_i, spvnfolder)

        # save spin informaiton
        if para_dict["compute_polarization"]:
            for ipattern in spin_filepattern:
                spin_list = glob(path.join(spinfolder, ipattern))
                for ispinfile in spin_list:
                    if path.isfile(ispinfile):
                        shutil.move(ispinfile, spvnfolder)


        hf = h5py.File("{0}.h5".format(results_name), "w")
        gtemp = hf.create_group("{0}".format(results_name))
        file_list = glob(path.join(spvnfolder, "*"))
        for file_path in file_list:
            file_name = file_path.split("/")[-1]
            if "usedParameters" in file_name:
                parafile = open(file_path)
                for iline, rawline in enumerate(parafile.readlines()):
                    paraline = rawline.strip('\n')
                    gtemp.attrs.create("{0}".format(iline),
                                       np.string_(paraline))
            else:
                dtemp = np.loadtxt(file_path, dtype="float32")
                h5data = gtemp.create_dataset("{0}".format(file_name),
                                              data=dtemp,
                                              compression="gzip",
                                              compression_opts=9)
                # save header
                ftemp = open(file_path, "r")
                header_text = str(ftemp.readline())
                ftemp.close()
                if header_text.startswith("#"):
                    h5data.attrs.create("header", np.string_(header_text))
        hf.close()
        shutil.move("{}.h5".format(results_name), final_results_folder)
        shutil.rmtree(spvnfolder, ignore_errors=True)
    else:
        print("{} is broken, skipped".format(spvnfolder), flush=True)
    return (status)


def remove_unwanted_outputs(final_results_folder, event_id, para_dict):
    """
        This function removes all hydro surface file and UrQMD results
        if they are unwanted to save space

    """
    if not para_dict["save_ipglasma"]:
        ipglasmafolder = path.join(final_results_folder,
                                   "ipglasma_results_{}".format(event_id))
        shutil.rmtree(ipglasmafolder, ignore_errors=True)

    if not para_dict["save_kompost"]:
        kompostfolder = path.join(final_results_folder,
                                  "kompost_results_{}".format(event_id))
        shutil.rmtree(kompostfolder, ignore_errors=True)

    if not para_dict["save_hydro"]:
        hydrofolder = path.join(final_results_folder,
                                "hydro_results_{}".format(event_id))
        shutil.rmtree(hydrofolder, ignore_errors=True)

    if para_dict["compute_polarization"]:
        spinfolder = path.join(final_results_folder,
                               "spin_results_{}".format(event_id))
        shutil.rmtree(spinfolder, ignore_errors=True)

    if not para_dict["save_urqmd"]:
        urqmd_results_name = "particle_list_{}.bin".format(event_id)
        remove(path.join(final_results_folder, urqmd_results_name))

    if para_dict["compute_photons"]:
        photonfolder = path.join(final_results_folder,
                                "photon_results_{}".format(event_id))
        shutil.rmtree(photonfolder, ignore_errors=True)


def checkPoint(startTime, checkPointFileName, finalResultsFolder):
    checkPointTime = time.time()
    if checkPointTime - startTime > 43200:
        # trigger the checkpoint when the simulation runs longer than 12 hours
        if path.exists(checkPointFileName):
            remove(checkPointFileName)
        tar = tarfile.open(checkPointFileName, 'w:gz')
        tar.add(finalResultsFolder)
        tar.close()
        sys.exit(85)


def main(para_dict_):
    """This is the main function"""
    startTime = time.time()
    initial_condition = para_dict_['initial_condition']
    initial_type = para_dict_['initial_type']
    num_threads = para_dict_['num_threads']
    n_urqmd = para_dict_['n_urqmd']
    curr_time = time.asctime()
    print("\U0001F3CE  [{}] Number of threads: {}".format(
        curr_time, num_threads),
          flush=True)

    idx0 = para_dict_['hydro_id0']
    nev = para_dict_['n_hydro']
    exitErrorTrigger = False
    exitErrorTriggerInitial = False
    for iev in range(nev):
        curr_time = time.asctime()

        event_id = str(iev + idx0)
        if (para_dict_['initial_condition'] != "self"
            and para_dict_['initial_condition'] != "fixCentrality"):
            initial_database_name = (
                    initial_condition.split("/")[-1].split(".h5")[0])
            event_id = initial_database_name + "_" + event_id

        final_results_folder = "EVENT_RESULTS_{}".format(event_id)

        # setup OSG checkpoint file
        CHECKPOINT_FILENAME = "{}.tar.gz".format(final_results_folder)
        try:
            tar = tarfile.open("{}".format(CHECKPOINT_FILENAME), 'r:gz')
            tar.extractall()
            tar.close()
            # remove the tar file to save disk space
            remove(CHECKPOINT_FILENAME)
        except FileNotFoundError:
            pass

        if path.exists(final_results_folder):
            print("{} exists ...".format(final_results_folder), flush=True)
            results_file = path.join(final_results_folder,
                                     "spvn_results_{}.h5".format(event_id))
            status = False
            if path.exists(results_file):
                status = True
            else:
                spvnfolder = path.join(final_results_folder,
                                       "spvn_results_{}".format(event_id))
                if path.exists(spvnfolder):
                    status = check_an_event_is_good(spvnfolder)
            if status:
                print(
                    "{} finished properly. No need to rerun.".format(event_id),
                    flush=True)
                continue
            print("Rerun {} ...".format(final_results_folder), flush=True)
        else:
            mkdir(final_results_folder)
        print("[{}] Generate initial condition ... ".format(curr_time),
              flush=True)

        initStauts, ifile = get_initial_condition(initial_condition,
                                                  initial_type, iev,
                                                  idx0 + iev,
                                                  para_dict_['seed_add'],
                                                  final_results_folder,
                                                  para_dict_['time_stamp_str'])
        if not initStauts:
            exitErrorTriggerInitial = True
            continue

        if initial_type == "3DMCGlauber_consttau":
            filename = ifile.split("/")[-1]
            filepath = initial_condition
            shutil.copy(path.join(filepath, filename),
                        "MUSIC/initial/initial_TA.dat")
            shutil.copy(path.join(filepath, re.sub("TA", "TB", filename)),
                        "MUSIC/initial/initial_TB.dat")

        if initial_type == "IPGlasma+KoMPoST":
            kompost_success, kompost_folder_name = run_kompost(
                final_results_folder, event_id)
            hydro_initial_file = "MUSIC/initial/epsilon-u-Hydro.dat"
            if path.islink(hydro_initial_file):
                remove(hydro_initial_file)
            call("ln -s {0:s} {1:s}".format(
                path.join(path.abspath(final_results_folder),
                          kompost_folder_name,
                          ("ekt_tIn01_tOut08"
                           + ".music_init_flowNonLinear_pimunuTransverse.txt")),
                hydro_initial_file),
                 shell=True)

        # first run hydro
        hydro_success, hydro_folder_name = run_hydro_event(
            final_results_folder, event_id)

        if not hydro_success:
            # if hydro didn't finish properly, just skip this event
            print("\U000026D4  {} did not finsh properly, skipped.".format(
                hydro_folder_name),
                  flush=True)
            exitErrorTrigger = True
            continue

        if (initial_type == "3DMCGlauber_dynamical"
                and (initial_condition == "self" or "fixCentrality")):
            # save the initial condition
            shutil.move("MUSIC/initial/strings.dat",
                        path.join(final_results_folder, hydro_folder_name,
                                  "strings_{}.dat".format(event_id)))

        if para_dict_["check_point_flag"]:
            checkPoint(startTime, CHECKPOINT_FILENAME, final_results_folder)

        if para_dict_['compute_photons']:
            # if hydro finishes properly, we continue to do photon radiation
            prepare_evolution_files_for_photon(final_results_folder,
                                               hydro_folder_name)
            photon_success, photon_folder_name = run_photon(
                                final_results_folder, event_id)
            if not photon_success:
                exitErrorTrigger = True
                continue
            if not para_dict["save_hydro"]:
                evoFileName = path.join(final_results_folder,
                                        hydro_folder_name,
                                        "evolution_all_xyeta.dat")
                shutil.rmtree(evoFileName, ignore_errors=True)
            if para_dict_["check_point_flag"]:
                checkPoint(startTime, CHECKPOINT_FILENAME,
                           final_results_folder)

        nUrQMDFolder = n_urqmd
        if para_dict_["compute_polarization"]:
            nUrQMDFolder += 1
        # if hydro finishes properly, we continue to do hadronic transport
        status_success = prepare_surface_files_for_urqmd(final_results_folder,
                                                         hydro_folder_name,
                                                         nUrQMDFolder)
        if not status_success:
            exitErrorTrigger = True
            continue

        # then run UrQMD events in parallel
        urqmd_success, urqmd_file_path = run_urqmd_shell(
            n_urqmd, final_results_folder, event_id, para_dict_,
            startTime, CHECKPOINT_FILENAME)
        if not urqmd_success:
            print("\U000026D4  {} did not finsh properly, skipped.".format(
                urqmd_file_path),
                  flush=True)
            continue

        # finally collect results
        run_spvn_analysis(urqmd_file_path, num_threads, final_results_folder,
                          event_id)

        # zip results into a hdf5 database
        status = zip_results_into_hdf5(final_results_folder, event_id,
                                       para_dict_)

        # remove the unwanted outputs if event is finished properly
        if status:
            remove_unwanted_outputs(final_results_folder, event_id, para_dict_)

    if exitErrorTriggerInitial:
        sys.exit(71)

    if exitErrorTrigger:
        sys.exit(73)


if __name__ == "__main__":
    try:
        INITIAL_CONDITION_TYPE = str(sys.argv[1])
        INITIAL_CONDITION_DATABASE = str(sys.argv[2])
        N_HYDRO_EVENTS = int(sys.argv[3])
        HYDRO_EVENT_ID0 = int(sys.argv[4])
        N_URQMD = int(sys.argv[5])
        N_THREADS = int(sys.argv[6])
        SAVE_IPGLASMA = (sys.argv[7].lower() == "true")
        SAVE_KOMPOST = (sys.argv[8].lower() == "true")
        SAVE_HYDRO = (sys.argv[9].lower() == "true")
        SAVE_URQMD = (sys.argv[10].lower() == "true")
        SEED_ADD = int(sys.argv[11])
        TIME_STAMP = str(sys.argv[12])
        COMP_POLARIZATION = (sys.argv[13].lower() == "true")
        COMP_PHOTONS = (sys.argv[14].lower() == "true")
        CHECK_POINT = (sys.argv[15].lower() == "true")
        AFTERBURNER_TYPE = str(sys.argv[16])
    except IndexError:
        print_usage()
        sys.exit(0)

    known_initial_types = [
        "IPGlasma", "IPGlasma+KoMPoST",
        "3DMCGlauber_dynamical", "3DMCGlauber_participants",
        "3DMCGlauber_consttau",
    ]
    if INITIAL_CONDITION_TYPE not in known_initial_types:
        print("\U0001F6AB  "
              + "Do not recognize the initial condition type: {}".format(
                  INITIAL_CONDITION_TYPE),
              flush=True)
        sys.exit(1)

    para_dict = {
        'initial_condition': INITIAL_CONDITION_DATABASE,
        'initial_type': INITIAL_CONDITION_TYPE,
        'n_hydro': N_HYDRO_EVENTS,
        'hydro_id0': HYDRO_EVENT_ID0,
        'n_urqmd': N_URQMD,
        'num_threads': N_THREADS,
        'save_ipglasma': SAVE_IPGLASMA,
        'save_kompost': SAVE_KOMPOST,
        'save_hydro': SAVE_HYDRO,
        'save_urqmd': SAVE_URQMD,
        'compute_polarization': COMP_POLARIZATION,
        'compute_photons': COMP_PHOTONS,
        'seed_add': SEED_ADD,
        'time_stamp_str': TIME_STAMP,
        'check_point_flag': CHECK_POINT,
        'afterburner_type': AFTERBURNER_TYPE,
    }

    main(para_dict)
