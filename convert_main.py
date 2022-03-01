#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 5, 2021

@author: Liz Verbeek

This script is part of the TC risk model developed as part of a Master Thesis 
for the Master's Programme Computational Science at the University of Amsterdam, 
see https://github.com/lizverbeek/global_TC_risk_model .

Together with convert_functions.py, this script converts all storm and exposure
files to GeoTIFF files and matches their projection and resolution by upsampling
the files with the lowest resolution.

Please note that the resulting files (after resolution matching) take up around
50 to 100 GB of memory per dataset (historical, current or future climate 
conditions). It is therefore recommended to apply only one of these conversions 
at a time.

"""

import os
import gdal

from convert_functions import netcdf_to_geotiff
from convert_functions import convert_storms_RP
from convert_functions import reproject_GHSL
from convert_functions import match_resolution


# ==============================================================================
nc_variable = "mean"			# NetCDF variable (mean, stdev, conv_5 or conv_95)
crs = 4326						# CRS for conversion
# ==============================================================================


# ==============================================================================
# 								Get paths
# ==============================================================================
current_dir = os.getcwd()
top_dir = os.path.abspath(os.path.join(current_dir, os.pardir))
storm_dir = top_dir + "/data/storm_data"
current_storm_dir = storm_dir + "/Current_climate"
future_storm_dir = (storm_dir + "/Future_climate/ALL_MODELS_MEDIAN/GeoTIFF/EPSG_" 
				    + str(crs) + "/Original_Res")

# ==============================================================================
# Choose directory of exposure data
# ==============================================================================
# exposure_dir = top_dir + "/data/GHSL_data/World_2000"
exposure_dir = top_dir + "/data/GHSL_data/World_2015"

surface_file_dir = top_dir + "/data/land_surface_data/Original_Res"
historical_storm_dir = (storm_dir + "/Storms_with_damage_data/EPSG_" + 
						str(crs) + "/Original_Res")


# ============================================================================
# Reproject GHSL data
# ============================================================================
exposure_dir = reproject_GHSL(exposure_dir, crs)

# ============================================================================
# Convert storm NetCDFs to GeoTIFF with desired crs
# ============================================================================
current_storm_dir = convert_storms_RP(current_storm_dir, nc_variable, crs)


# # ============================================================================
# # Use this for folder structure when skipping conversion and GHSL reprojection
# # ============================================================================
# current_storm_dir = (current_storm_dir + "/GeoTIFF/EPSG_" + str(crs) 
# 						 + "/Original_Res")
# exposure_dir = exposure_dir + "/EPSG_" + str(crs)


# ============================================================================
# Match resolution of storm data to GHSL resolution
# ============================================================================
ref_file = exposure_dir + "/" + os.listdir(exposure_dir)[0]
print("Matching storm files to exposure resolution")

# Historical storms
match_resolution(ref_file, historical_storm_dir, 'average')

# # Return period storms, current climate
match_resolution(ref_file, current_storm_dir, "average")

# Return period storms, future climate, median of models
match_resolution(ref_file, future_storm_dir, "average")


# ============================================================================
# Match resolution of surface area file to GHSL resolution
# ============================================================================
print("Matching surface area files to exposure resolution")
match_resolution(ref_file, surface_file_dir, "max")

