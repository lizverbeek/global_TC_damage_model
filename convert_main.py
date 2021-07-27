#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 5, 2021

@author: Liz Verbeek

Conversion of NetCDF storm data to GeoTIFF files, resolution matching of 
GHSL data, surface area data and storm files. Conversion is included for all 
historical storms and current and future climate NetCDF files.

Please note that the resulting files (after resolution matching) take up around
50 to 100 GB of memory per dataset (historical, current climate conditions or any
of the 4 different future climate dataset). It is therefore recommended to apply
only one conversion at a time.
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

# # Historical storms
# match_resolution(ref_file, historical_storm_dir, 'average')

# # Return period storms, current climate
match_resolution(ref_file, current_storm_dir, "average")

# # Return period storms, future climate, median of models
# match_resolution(ref_file, future_storm_dir, "average")


# ============================================================================
# Match resolution of surface area file to GHSL resolution
# ============================================================================
print("Matching surface area files to exposure resolution")
match_resolution(ref_file, surface_file_dir, "max")

