#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
@author: Liz Verbeek

This script is part of the TC risk model developed as part of a Master Thesis 
for the Master's Programme Computational Science at the University of Amsterdam, 
see https://github.com/lizverbeek/global_TC_risk_model .

This script reads the IBTrACS data, performs preprocessing to select suitable 
TC tracks and converts these to 2D windfields.

"""

import os
import numpy as np
import xarray as xr
import multiprocessing as mp

import ibtracs_preprocessing
import holland_model as hm
import windfield_functions


# ============================================================================== 
# Get paths for loading and saving data
# ============================================================================== 
current_dir = os.getcwd()
top_dir = os.path.abspath(os.path.join(current_dir, os.pardir))
ibtracs_dir = top_dir + "/data/storm_data/IBTrACS"
preprocessed_dir = ibtracs_dir + "/IBTrACS_preprocessed"
output_dir = ibtracs_dir + "/Historical_Storms"

if not os.path.exists(preprocessed_dir):
    print("Creating directory " + preprocessed_dir)
    os.makedirs(preprocessed_dir)
if not os.path.exists(output_dir):
    print("Creating directory " + output_dir)
    os.makedirs(output_dir)


# ============================================================================== 
# Open IBTrACS dataset
# ============================================================================== 
print("-------------------------------------------------------------")
print("Preprocessing IBTrACS data ...... ")
print("-------------------------------------------------------------")

# Preprocessing
data_path = ibtracs_dir + "/IBTrACS.since1980.v04r00.nc"
data = xr.open_dataset(data_path, decode_times=False)
ibtracs_preprocessing.extract_data(data, preprocessed_dir)
data.close()

# Load preprocessed data
latlist = np.load(preprocessed_dir + "/LATLIST_INTERP.npy", allow_pickle=True).item()
lonlist = np.load(preprocessed_dir + "/LONLIST_INTERP.npy", allow_pickle=True).item()
timelist = np.load(preprocessed_dir + "/TIMELIST_INTERP.npy", allow_pickle=True).item()
windlist = np.load(preprocessed_dir + "/WINDLIST_INTERP.npy", allow_pickle=True).item()
preslist = np.load(preprocessed_dir + "/PRESLIST_INTERP.npy", allow_pickle=True).item()
rmaxlist = np.load(preprocessed_dir + "/RMAXLIST_INTERP.npy", allow_pickle=True).item()
basinlist = np.load(preprocessed_dir + "/BASINLIST_INTERP.npy", allow_pickle=True).item()
yearlist = np.load(preprocessed_dir + "/YEARLIST_INTERP.npy", allow_pickle=True).item()
namelist = np.load(preprocessed_dir + "/NAMELIST_INTERP.npy", allow_pickle=True).item()


# ============================================================================== 
# Run parallel wind field estimation
# ============================================================================== 

# Change structure to list of dicts per track for easier multiprocessing
stormlist = [{'nr': i, 'lat': latlist[i], 'lon': lonlist[i], 'time': timelist[i], 
              'wind': windlist[i], 'pres': preslist[i], 'rmax': rmaxlist[i],
              'basin': basinlist[i][0]} for i in latlist.keys()]

# Run wind field computation
n_cores = mp.cpu_count() - 2
print("-------------------------------------------------------------")
print("Running track to wind field conversion on " + str(n_cores) + " cores")
print("-------------------------------------------------------------")
print()
pool = mp.Pool(n_cores)
results = pool.map(windfield_functions.compute_windfield, stormlist)
pool.close()

print("----------------- Wind field conversion finished ------------------")
print()
