#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 5, 2021

@author: Liz Verbeek

This script is part of the TC risk model developed as part of a Master Thesis 
for the Master's Programme Computational Science at the University of Amsterdam, 
see https://github.com/lizverbeek/global_TC_risk_model .

This script contains all preprocessing steps for the future climate TC hazard
datasets.
The input should be given as NetCDF with different return period bands. These are
converted to GeoTIFF files per return period of the median of all input datasets.

"""

import os
import gdal
import numpy as np
import netCDF4 as nc

from osgeo import osr


# ==============================================================================
# Get paths
# ==============================================================================
current_dir = os.getcwd()
top_dir = os.path.abspath(os.path.join(current_dir, os.pardir))
storm_dir = top_dir + "/data/storm_data"
future_storm_dir = storm_dir + "/Future_climate"

output_dir = future_storm_dir + "/ALL_MODELS_MEDIAN/GeoTIFF/EPSG_4326/Original_Res"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)


# ==============================================================================
# Get GeoTIFFs for each return period, median of all datasets
# ==============================================================================
basins = ["EP", "NA", "NI", "SI", "SP", "WP"]
crs = 4326
for basin in basins:
    wind_speeds_all_models = []
    for model_dir in os.listdir(future_storm_dir):
        if model_dir.startswith("STORM"):
            for model_dataset in os.listdir(future_storm_dir + "/" 
                                            + model_dir + "/NetCDF"):

                # Get all model results per basin
                if model_dataset.endswith(basin + ".nc"):
                    print(model_dataset)
                    ds = nc.Dataset(future_storm_dir + "/" + model_dir 
                                    + "/NetCDF" + "/" + model_dataset, "r")
                    mean_wind_speeds = np.array(ds.variables["mean"])
                    wind_speeds_all_models.append(mean_wind_speeds)

    # Take median of 4 models
    median_wind_speeds = np.median(wind_speeds_all_models, axis=0)
    median_wind_speeds = np.flipud(median_wind_speeds)

    # Write median wind speeds to GeoTIFF files per return period.
    rows = ds.dimensions["lat"].size
    cols = ds.dimensions["lon"].size
    lat_min = np.min(ds.variables["lat"][:])
    lat_max = np.max(ds.variables["lat"][:])
    lon_min = np.min(ds.variables["lon"][:])
    lon_max = np.max(ds.variables["lon"][:])
    if basin == "NA" or basin == "EP":
        lon_min -= 360
        lon_max -= 360
    lat_resolution = (lat_max - lat_min) / (rows - 1)
    lon_resolution = (lon_max - lon_min) / (cols - 1)

    # Create GeoTIFF file for each return period.
    driver = gdal.GetDriverByName("GTiff")
    source = osr.SpatialReference()
    source.ImportFromEPSG(crs)
    output_proj = source.ExportToPrettyWkt()

    for idx, return_period in enumerate(ds.variables["rp"][:]):
        output_file = (output_dir + "/STORM_RP" + str(int(return_period))
                       + "_" + str(basin) + ".tif")
        output = driver.Create(output_file, cols, rows, 1, gdal.GDT_Float64)
        output.SetProjection(output_proj)
        output.SetGeoTransform((lon_min - 0.5 * lon_resolution, lon_resolution, 0,
                               lat_max + 0.5 * lat_resolution, 0, -lat_resolution))
        output.GetRasterBand(1).WriteArray(median_wind_speeds[:,:,idx])
