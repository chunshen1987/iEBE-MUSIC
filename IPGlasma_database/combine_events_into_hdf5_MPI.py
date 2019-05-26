#!/usr/bin/env python
"""
     This script combines pre-generated IP-Glasma events into a hdf5 database.
"""

import sys
from os import path, system
from glob import glob
import h5py
from mpi4py import MPI

import numpy as np

def print_help():
    """This function prints out help message"""
    print("{0} results_folder".format(sys.argv[0]))

def collect_one_IPGlasma_event(results_path, event_path, hf):
    """This function collects one IPGlasma event"""
    event_id = event_path.split("/")[-1].split("Parameters")[-1].split(".")[0]
    gtemp = hf.create_group("event-{0}".format(event_id))

    file_name = "usedParameters{0}.dat".format(event_id)
    parafile = open(path.join(results_path, file_name))
    for iline, rawline in enumerate(parafile.readlines()):
        paraline = rawline.strip('\n')
        gtemp.attrs.create("{0}".format(iline), np.string_(paraline))

    file_name = "NcollList{0}.dat".format(event_id)
    filepath = path.join(results_path, file_name)
    if path.isfile(filepath):
        dtemp = np.loadtxt(filepath)
        dtemp = np.nan_to_num(dtemp).reshape(-1, 2)
        dset = gtemp.create_dataset("{0}".format(file_name), data=dtemp,
                                    compression="gzip", compression_opts=9)
    file_name = "NpartList{0}.dat".format(event_id)
    filepath = path.join(results_path, file_name)
    if path.isfile(filepath):
        dtemp = np.loadtxt(path.join(results_path, file_name))
        dtemp = np.nan_to_num(dtemp).reshape(-1, 4)
        dset = gtemp.create_dataset("{0}".format(file_name), data=dtemp,
                                    compression="gzip", compression_opts=9)

    file_name_pattern = "NpartdNdy-t"
    filelist = glob(path.join(results_path, "{0}*-{1}.dat".format(
        file_name_pattern, event_id)))
    for filepath in filelist:
        filename = filepath.split("/")[-1]
        dtemp = np.genfromtxt(filepath, dtype='str')
        data = np.zeros(len(dtemp))
        for idx in range(len(dtemp)):
            if dtemp[idx] != 'N/A':
                data[idx] = float(dtemp[idx])
            else:
                data[idx] = 0.0
        dset = gtemp.create_dataset("{0}".format(filename), data=data,
                                    compression="gzip", compression_opts=9)

    file_name_pattern = "epsilon-u-Hydro-t"
    filelist = glob(path.join(results_path, "{0}*-{1}.dat".format(
        file_name_pattern, event_id)))
    for ifile, filepath in enumerate(filelist):
        filename = filepath.split("/")[-1]
        dtemp    = np.loadtxt(filepath)
        dtemp    = np.nan_to_num(dtemp)
        x_size   = abs(dtemp[0, 1])*2.
        y_size   = abs(dtemp[0, 2])*2.
        data_cut = dtemp[:, 3:]
        dset     = gtemp.create_dataset("{0}".format(filename),
                                        data=data_cut,
                                        compression="gzip", compression_opts=9)
        f = open(filepath)
        header = f.readline().strip('\n')
        dset.attrs.create("header", np.string_(header))
        tmp = header.split()
        dx = float(tmp[12])
        dy = float(tmp[14])
        nx = int(tmp[6])
        ny = int(tmp[8])
        dset.attrs.create("x_size", x_size)
        dset.attrs.create("y_size", y_size)
        dset.attrs.create("dx", dx)
        dset.attrs.create("dy", dy)
        dset.attrs.create("nx", nx)
        dset.attrs.create("ny", ny)

def collect_IPGlasma_events(results_folder):
    """This function collects IPGlasma events in results_folder"""
    mpi_comm = MPI.COMM_WORLD
    mpi_rank = mpi_comm.Get_rank()
    mpi_size = mpi_comm.Get_size()

    if mpi_rank == 0:
        print("MPI using {} threads ...".format(mpi_size))

    results_name = results_folder.split("/")[-1]
    if results_name == "":
        results_name = results_folder.split("/")[-2]
    results_path = path.abspath(path.join(".", results_folder))
    event_list = glob(path.join(results_path, "usedParameters*.dat"))
    nev = len(event_list)

    h5filename = "{0}_rank{1}.h5".format(results_name, mpi_rank)
    print("MPI rank {0}: collect {1} to {2} ... ".format(
        mpi_rank, results_folder, h5filename))
    hf = h5py.File(h5filename, "w")
    mpi_comm.Barrier()
    for ievent in range(nev):
        if mpi_rank == ievent%mpi_size:
            event_path = event_list[ievent]
            print("MPI rank {0:d} processing {1:d}/{2:d} ... ".format(
                mpi_rank, ievent, nev))
            collect_one_IPGlasma_event(results_path, event_path, hf)
    hf.close()
    mpi_comm.Barrier()


    # now combine all the hdf5 files into one using rank 0
    if mpi_rank == 0:
        print("combining to one hdf5 file {}.h5 ...".format(results_name))
        h5_filelist = glob("./{}_rank*.h5".format(results_name))
        for filename in h5_filelist:
            print("processing {0} ... ".format(filename))
            hftemp = h5py.File(filename, "r")
            glist = list(hftemp.keys())
            hftemp.close()
            for gtemp in glist:
                system('h5copy -i {0} -o {1}.h5 -s {2} -d {2}'.format(
                    filename, results_name, gtemp))
            system('rm -fr {0}'.format(filename))


if __name__ == "__main__":
    try:
        RESULTS_FOLDER = str(sys.argv[1])
        collect_IPGlasma_events(RESULTS_FOLDER)
    except IndexError:
        print_help()
        exit(1)
