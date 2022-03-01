#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""

@author: Liz Verbeek

This script is part of the TC risk model developed as part of a Master Thesis 
for the Master's Programme Computational Science at the University of Amsterdam, 
see https://github.com/lizverbeek/global_TC_risk_model .

Together with damage_functions.py, this script contains the final damage computation.
From the command line, it should be specified which type of vulnerability function 
should be used ("sigmoid" or "power"); the country of interest (ISO code); and which
climate scenario ("current" or "future") for the TC hazard datasets should 
form the input of the model.

"""

import os
import sys
import fiona
import numpy as np
import pandas as pd
import multiprocessing as mp

from damage_functions import get_region_shape
from damage_functions import get_overlapping_storms
from damage_functions import damage_per_region


# ==============================================================================
vul_func = sys.argv[1]          # "sigmoid" or "power" vulnerability function
country = sys.argv[2]           # Country for damage computation
scenario = sys.argv[3]          # "current" or "future" climate conditions
# ==============================================================================
crs = 4326                      # CRS of input data
level = 1                       # Region level for which to aggregate results
year = "2015"                   # Year of exposure data, used for normalization
                                # of monetary values
n_cores = mp.cpu_count() - 2    # Number of cores to use for parallel computing
# ==============================================================================

# ============================================================================== 
# Get paths for storm data
# ============================================================================== 
current_dir = os.getcwd()
top_dir = os.path.abspath(os.path.join(current_dir, os.pardir))

# Current climate conditions
if scenario == "current":
    storm_dir = (top_dir + "/data/storm_data/Current_climate/GeoTIFF/EPSG_" 
                 + str(crs) + "/World_Res")
    storms = [storm.split(".")[-2] for storm in os.listdir(storm_dir)]

# Future climate conditions (median of future climate projections)
if scenario == "future":
    storm_dir = (top_dir + 
                 "/data/storm_data/Future_climate/ALL_MODELS_MEDIAN/GeoTIFF/EPSG_"
                 + str(crs) + "/World_Res")
    storms = [storm.split(".")[-2] for storm in os.listdir(storm_dir)]


# Exposure, surface area, region and max damage paths
exposure_dir = (top_dir + "/data/GHSL_data/World_" + str(year) 
                + "/EPSG_" + str(crs))
exposure_map = os.listdir(exposure_dir)[0]
exposure_path = exposure_dir + "/" + exposure_map

surface_path = (top_dir + "/data/land_surface_data/World_Res"
               + "/land_surface_km2.tif")
region_file = top_dir + "/data/NUTS_WORLD/gadm36_2.shp"
max_damage_path = top_dir + "/data/max_damage_data/max_damage_GDP_corrected.csv"
reported_damage_path = (top_dir + "/data/storm_data/Storms_with_damage_data/" + 
                        "storm_damage_matches_ALL.csv")

results_dir = top_dir + "/results"


# ==============================================================================
# Read maximum damage values for all countries and given exposure year
# ==============================================================================
max_damage_df = pd.read_csv(max_damage_path, 
                            usecols=["Country", "RES " + str(year), 
                                     "COM " + str(year), "IND " + str(year)])
max_damage_df = max_damage_df.rename(columns={"RES " + str(year): "RES",
                                              "COM " + str(year): "COM",
                                              "IND " + str(year): "IND"})
max_damage_df = max_damage_df.set_index("Country")
max_damage_dict = max_damage_df.to_dict(orient='index')

# ==============================================================================
# Get countries, country shapes and country codes
# ==============================================================================

# Match countries to specified regions (used for calibration)
basin_dict = {"AUS":"OC", "BGD":"NI", "CAN":"NA2", "CHN":"WP3", "CRI":"NA1", 
              "CUB":"NA1", "DOM":"NA1", "FJI":"OC", "GTM":"NA1", "HND":"NA1", 
              "HTI":"NA1", "IDN":"WP1", "IND":"NI", "JPN":"WP4", "KOR":"WP4", 
              "LAO":"WP1", "LKA":"NI", "MAC":"WP4", "MDG":"SI", "MEX":"NA1", 
              "MMR":"NI", "MOZ":"SI", "NIC":"NA1", "OMN":"NI", "PAK":"NI", 
              "PHL":"WP2", "PRK":"WP4", "SLV":"NA1", "SOM":"SI", "THA":"WP1", 
              "USA":"NA2", "VNM":"WP1", "VUT":"OC", "YEM":"NI", "ZWE":"SI"}

# Region information all countries 
basin_dict_new = {"AFG": "NI", "BTN": "NI", "KHM": "WP1", "COL": "NA1", 
                  "DJI": "NI", "ERI": "NI", "ETH": "NI", "GUY": "NA1", 
                  "IRN": "NI", "MWI": "SI", "MYS": "WP1", "MNG": "WP3", 
                  "NPL": "NI", "NZL": "OC", "PAN": "NA1", "PNG": "OC", 
                  "SLB": "OC", "ZAF": "SI", "SUR": "NA1", "SWZ": "SI", 
                  "TZA": "SI", "TLS": "OC", "ARE": "NI", "VEN": "NA1",  
                  "AGO": "SI", "BRN": "WP1", "BWA": "SI", "NAM": "SI", 
                  "RUS": "WP4", "VIR": "NA1", "ZMB": "SI"}

# Region information for uncalibrated regions
basin_dict_other = {"BFA": "NA3", "CIV": "NA3", "DZA": "NA3", "ESP": "NA4", 
                    "FRA": "NA4", "GBR": "NA4", "GHA": "NA3", "GIN": "NA3", 
                    "GMB": "NA3", "GNB": "NA3", "IMN": "NA4", "LBR": "NA3", 
                    "MAR": "NA3", "MLI": "NA3", "MRT": "NA3", "PRT": "NA4", 
                    "SEN": "NA3", "SLE": "NA3"}

# Merge region definitions for full damage computation
basin_dict = {**basin_dict, **basin_dict_new, **basin_dict_other}

# # Run all countries with reported damage
# reported_damage_df = pd.read_csv(reported_damage_path, usecols=["Nr", "ISO"])
# countries = set(reported_damage_df["ISO"]).intersection(max_damage_dict.keys())

# Match country codes to country names
region_shp = fiona.open(region_file, "r")
country_codes = [region["properties"]["GID_0"] for region in region_shp]
country_names = [region["properties"]["NAME_0"] for region in region_shp]
country_code_dict = dict(zip(country_codes, country_names))
country_name_dict = dict(zip(country_names, country_codes))

if level == 1:
    level1_codes = [region["properties"]["GID_1"] for region in region_shp]
    level1_names = [region["properties"]["NAME_1"] for region in region_shp]
    level1_name_dict = dict(zip(level1_names, level1_codes))
elif level == 2:
    level2_codes = [region["properties"]["GID_1"] for region in region_shp]
    level2_names = [region["properties"]["NAME_1"] for region in region_shp]
    level2_name_dict = dict(zip(level2_names, level2_codes))


# ==============================================================================
# Create directories for saving temporary files
# ==============================================================================
total_damage = {}
temp_dir = (top_dir + "/data/temp_cropped_files")
for storm in storms:

    # Create paths for saving temporary data
    storm_path = storm_dir + "/" + storm
    new_dir = temp_dir + "/" + storm_path.split("/")[-1]
    if not os.path.exists(new_dir):
        os.makedirs(new_dir)

    # Create dict for saving damage data
    total_damage[storm] = {}


# ==============================================================================
# Compute damage for given country
# ==============================================================================
print("--------------------------- COUNTRY " + country_code_dict[country] + " --------------------------")
print("------------------------- Running on " + str(n_cores) + " cores ---------------------")
print()

if vul_func == "sigmoid":
    # Run model with regional optimal sigmoidal function vthresh parameters
    if (basin_dict[country] == "OC" or basin_dict[country] == "SI" 
        or basin_dict[country] == "NA3" or basin_dict[country] == "NA4"):
        # Use global vthresh parameter
        param1 = 19.5
    elif basin_dict[country] == "NA1":
        param1 = 11.5
    elif basin_dict[country] == "NA2":
        param1 = 15.0
    elif basin_dict[country] == "NI":
        param1 = 20.5
    elif basin_dict[country] == "WP1":
        param1 = 13.0
    elif basin_dict[country] == "WP2":
        param1 = 19.0
    elif basin_dict[country] == "WP3":
        param1 = 15.5
    elif basin_dict[country] == "WP4":
        param1 = 27.0
    # Vhalf is set to default value
    param2 = 65.7

elif vul_func == "power":
    # Run model with regional optimal power function parameters
    if basin_dict[country] == "NA1":
        param1 = 0.0115
        param2 = 3.5
    elif basin_dict[country] == "NA2":
        param1 = 0.0135
        param2 = 4.5
    elif basin_dict[country] == "NA3" or basin_dict[country] == "NA4":
        # global parameters
        param1 = 0.014
        param2 = 5.5
    elif basin_dict[country] == "NI":
        param1 = 0.011
        param2 = 3.5
    elif basin_dict[country] == "OC":
        param1 = 0.009
        param2 = 2.5
    elif basin_dict[country] == "SI":
        param1 = 0.0095
        param2 = 2.5
    elif basin_dict[country] == "WP1":
        param1 = 0.009
        param2 = 3.5
    elif basin_dict[country] == "WP2":
        param1 = 0.0135
        param2 = 5.0
    elif basin_dict[country] == "WP3":
        param1 = 0.0135
        param2 = 5.0
    elif basin_dict[country] == "WP4":
        param1 = 0.012
        param2 = 6.5
else:
    print("Unknown function type given, valid functions 'sigmoid' or 'power' ")

print("Running damage computation for " + vul_func 
      + " function, parameters: " + str(param1) 
      + ", " + str(param2))
print()

# ==============================================================================
# Get all country subregions for given level
# ==============================================================================
if level == 0:
    regions_level = [country_code_dict[country]]
elif level == 1:
    regions_level = [region["properties"]["GID_1"] for region in region_shp 
                     if region["properties"]["GID_0"] == country]
elif level == 2:
    regions_level = [region["properties"]["GID_2"] for region in region_shp 
                     if region["properties"]["GID_0"] == country]
else:
    print("Region level should be 0, 1 or 2!")
regions_level = set(regions_level)

# ==============================================================================
# First get storms that overlap with country
# ==============================================================================
print("Computing overlap for all storms ...... ")
country_shape = get_region_shape(country, 0, region_file)
storms_overlap = []
for storm in storms:
    storms_overlap = get_overlapping_storms(storms_overlap, storm, storm_dir, 
                                            country, country_shape, temp_dir)

# Set damage to zero for storms with no overlap
storms_no_overlap = set(storms) - set(storms_overlap)
for storm in storms_no_overlap:
    for region in regions_level:
        total_damage[storm][region] = 0.

        
print("---------------------------------------------------------------------")
print("Computing damage for storms: ")
print(storms_overlap)
print("---------------------------------------------------------------------")
print()
print("---------------------------------------------------------------------")
print("Computing damage for regions: ")
print(regions_level)
print("---------------------------------------------------------------------")
print()

# ============================================================================== 
# Crop files for given country and compute damage
# ============================================================================== 
for region in regions_level:

    # Get region shape for cropping
    if level == 0:
        region_shape = country_shape
    else:
        region_shape = get_region_shape(region, level, region_file)

    with mp.Manager() as manager:
        damage_dict = manager.dict()
        processes = []
        for storm in storms_overlap:
            p = mp.Process(target=damage_per_region,
                           args=(damage_dict, storm, storm_dir, country, 
                                 region, region_shape, exposure_path, 
                                 surface_path, max_damage_dict, temp_dir, 
                                 param1, param2, vul_func))
            p.start()
            processes.append(p)
            p.join()
        for storm, damage in damage_dict.items():
            total_damage[storm][region] = damage
    
print()

# ==============================================================================
# Show and save resulting damages
# ==============================================================================
damage_df = pd.DataFrame(total_damage)
damage_df = damage_df.reindex(sorted(damage_df.columns, reverse=True), axis=1)
if level == 0:
    damage_df["ISO"] = [country_name_dict[country] for country in damage_df.index.values]
    damage_df = damage_df.set_index("ISO")

print("----------------------------------------------------------------------")
print("------------------------ Resulting damage: ---------------------------")
print("----------------------------------------------------------------------")
print(damage_df)

damage_df.to_csv(results_dir + "/total_damage_" + vul_func + "_"
                 + str(param1).split(".")[0]
                 + str(param1).split(".")[1] + "_"
                 + str(param2).split(".")[0]
                 + str(param2).split(".")[1] + "_"
                 + country_code_dict[country] + ".csv")
