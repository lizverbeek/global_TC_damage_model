#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Created on Thu Mar 4, 2021

@author: Liz Verbeek
"""

import os
import re
import gdal
import numpy as np
import netCDF4 as nc

from osgeo import gdal, osr
from tqdm import tqdm

def netcdf_to_geotiff(input_path, nc_variable, crs, output_dir):
    """Convert NetCDF to GeoTIFF.

    Arguments:
        input_path {string}     -- string of data path to NetCDF file.
        nc_variable {string}    -- variable to extract from NetCDF.
        crs {int}               -- EPSG number.
        output_dir {string}     -- GeoTIFF output directory.

    """

    # Get basin from storm file name
    basin = input_path.split("/")[-1].split("_")[-1].split(".")[0]

    # Extract wind speeds from netcdf file
    ds = nc.Dataset(input_path, "r")
    wind_speeds = np.array(ds.variables[nc_variable])

    # Flip (up/down) wind speeds because of NetCDF conventions
    wind_speeds = np.flipud(wind_speeds)

    # Get minimum, maximum and size of latitide and longitude, 
    # longitude coordinates are translated -360 degrees for EP and NA basins.
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
        output.GetRasterBand(1).WriteArray(wind_speeds[:,:,idx])

    return


def convert_storms_RP(storm_dir, nc_variable, crs):
    """Convert storms with fixed return period to GeoTIFF files.

    Arguments:
        storm_dir {string}      -- Input storm files directory.
        nc_variable {string}    -- Variable to extract from NetCDFs.
        crs {int}               -- CRS of original file.

    Returns:
        output_dir {string}     -- Path to output directory

    """

    # Create storm output directory
    output_dir = storm_dir + "/GeoTIFF/EPSG_" + str(crs) + "/Original_Res"
    if not os.path.exists(output_dir):
        print("Creating directory " + output_dir)
        os.makedirs(output_dir)

    # Convert NetCDFs for fixed RP storm data
    input_dir = storm_dir + "/NetCDF"
    print("Creating storm GeoTIFFs")
    print("------------------------------------")
    for storm_file in tqdm(os.listdir(input_dir)):
        if re.search("STORM_FIXED_RETURN_PERIODS_*", storm_file):
            storm_path = input_dir + "/" + storm_file

            netcdf_to_geotiff(storm_path, nc_variable, crs, output_dir)

    print()

    return output_dir


def reproject_GHSL(exposure_dir, out_crs):
    """Reproject exposure GeoTIFF files to different CRS.

    Arguments:
        exposure_dir {string}   -- Exposure files directory.
        out_crs {int}           -- EPSG number for projection.
    
    Returns:
        output_dir {string}     -- Path to output directory

    """
    
    # Create new directory for projected files
    output_dir = exposure_dir + "/EPSG_" + str(out_crs)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Reproject all input files
    input_dir = exposure_dir + "/ESRI_54009"
    projection = "EPSG:" + str(out_crs)
    for in_file in os.listdir(input_dir):
        print("Reprojecting GHSL data to EPSG " + str(out_crs))
        print("------------------------------------")

        in_raster = gdal.Open(input_dir + "/" + in_file)
        out_file = in_file.replace("54009", str(out_crs))
        os.system("gdalwarp -t_srs " + projection + " -srcnodata 'None' " + 
                   input_dir + "/" + in_file + " " + output_dir + "/" + out_file)
        print()

    return output_dir


def match_resolution(ref_file, old_dir, resampleAlg):
    """"Match resolution of storm GeoTIFF to resolution of reference GeoTIFF.

    Arguments:
        ref_file {string}           -- Path to reference GeoTIFF file.
        old_dir {string}            -- Storm input directory.
        resampleAlg {string}        -- Sampling method for interpolation in warping

    """

    
    # Create new directory for high resolution files
    output_dir = (os.path.abspath(os.path.join(old_dir, os.pardir)) + 
                  "/World_Res")
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Get reference resolution
    ref_raster = gdal.Open(ref_file)
    ref_gt = ref_raster.GetGeoTransform()
    xRes = ref_gt[1]
    yRes = -ref_gt[5]

    # Create high resolution files
    print("------------------------------------")
    files = [file for file in os.listdir(old_dir) if file.split(".")[-1] == "tif"]
    for file in tqdm(files):
        old_file_path = os.path.join(old_dir, file)
        os.system("gdalwarp -tr " + str(xRes) + " " + str(yRes) + 
                  " -q -overwrite -r " + resampleAlg + " " +
                   old_file_path + " " + output_dir + "/" + file)

    print()

    return
