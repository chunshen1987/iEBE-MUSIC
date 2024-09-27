#! /usr/bin/env python
"""
     This script performs read in a minimum bias database and outputs
     the event id for different centrality bins
"""

from sys import argv, exit
from os import path, mkdir
from glob import glob
from numpy import *
import h5py
import shutil

centrality_cut_list = [
    0., 5., 10., 20., 30., 40., 50., 60., 70., 80., 90., 100.
]
try:
    data_path = path.abspath(argv[1])
    folder_name = data_path.split("/")[-1]
    data_filename = path.join(data_path, "{}.h5".format(folder_name))
    print("input data: {}".format(data_filename))
    hydro_surface_flag = False
    urqmd_flag = False
    hydro_folder = path.join(data_path, "HYDRO_RESULTS")
    urqmd_folder = path.join(data_path, "URQMD_RESULTS")
    if path.exists(hydro_folder):
        print("This run has hydro surface!")
        hydro_surface_flag = True
    if path.exists(urqmd_folder):
        print("This run has UrQMD outputs!")
        urqmd_flag = True
    if not hydro_surface_flag and not urqmd_flag:
        exit(0)
except IndexError:
    print("Usage: {} results_folder".format(argv[0]))
    exit(1)

hf = h5py.File(data_filename, "r")
event_list = list(hf.keys())
dN_dy_mb = []
for ifolder, event_name in enumerate(event_list):
    file_name = "particle_9999_vndata_eta_-0.5_0.5.dat"
    event_group = hf.get(event_name)
    try:
        temp_data = event_group.get(file_name)
        temp_data = nan_to_num(temp_data)
        dN_dy_mb.append(-temp_data[0, 1])
    except:
        continue
dN_dy_mb = -sort(array(dN_dy_mb))
print("total number of events: {}".format(len(dN_dy_mb)))

for icen in range(len(centrality_cut_list) - 1):
    if centrality_cut_list[icen + 1] < centrality_cut_list[icen]:
        continue

    hydro_directory_path = path.join(
        hydro_folder, "C{0:d}-{1:d}".format(int(centrality_cut_list[icen]),
                                            int(centrality_cut_list[icen + 1])))

    urqmd_directory_path = path.join(
        urqmd_folder, "C{0:d}-{1:d}".format(int(centrality_cut_list[icen]),
                                            int(centrality_cut_list[icen + 1])))

    if hydro_surface_flag:
        if path.exists(hydro_directory_path):
            shutil.rmtree(hydro_directory_path)
        mkdir(hydro_directory_path)

    if urqmd_flag:
        if path.exists(urqmd_directory_path):
            shutil.rmtree(urqmd_directory_path)
        mkdir(urqmd_directory_path)

    dN_dy_cut_high = (dN_dy_mb[int(
        len(dN_dy_mb)*centrality_cut_list[icen]/100.)])
    dN_dy_cut_low = dN_dy_mb[min(
        len(dN_dy_mb) - 1,
        int(len(dN_dy_mb)*centrality_cut_list[icen + 1]/100.))]

    selected_events_list = []
    for ifolder, event_name in enumerate(event_list):
        file_name = "particle_9999_vndata_eta_-0.5_0.5.dat"
        event_group = hf.get(event_name)
        temp_data = event_group.get(file_name)
        temp_data = nan_to_num(temp_data)
        if (temp_data[0, 1] > dN_dy_cut_low
                and temp_data[0, 1] <= dN_dy_cut_high):
            selected_events_list.append(event_name)

    nev = len(selected_events_list)
    print("analysis {}%-{}% nev = {}...".format(centrality_cut_list[icen],
                                                centrality_cut_list[icen + 1],
                                                nev))
    for ev_i_name in selected_events_list:
        event_id = ev_i_name.split("_")[-1]
        if hydro_surface_flag:
            hydro_event_name = "hydro_results_{}".format(event_id)
            shutil.move(path.join(hydro_folder, hydro_event_name),
                        hydro_directory_path)
        if urqmd_flag:
            urqmd_event_name = "particle_list_{}.gz".format(event_id)
            shutil.move(path.join(urqmd_folder, urqmd_event_name),
                        urqmd_directory_path)
