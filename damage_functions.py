#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

@author: Liz Verbeek
"""

import os
import fiona
import rasterio
import gdal
import numpy as np
import matplotlib.pyplot as plt

from rasterio.mask import mask
from shapely.geometry import shape
from shapely.ops import unary_union


def sigmoid_vulnerability_function(storm_map, v_thresh, v_half):
    """ Sigmoidal vulnerability function relating wind speeds 
        to building damage ratio.

    Arguments:
        storm_map {2d array}          -- Storm map containing wind speeds
        v_thresh {float}              -- Theoretical damage threshold value (in m/s)
        v_half {float}                -- Theoretical half damage value (in m/s)
    
    Returns:
        damage_ratio {2d array}       -- Ratio of building damage

    """
    zero = np.zeros(storm_map.shape)
    v_n = np.maximum(storm_map - v_thresh, zero) / (v_half - v_thresh)
    damage_ratio = v_n**3/(1 + v_n**3)

    return damage_ratio


def power_vulnerability_function(storm_map, factor, exp):
    """ Exponential vulnerability function.

    Arguments:
        storm_map {numpy ndarray}   -- Storm map containing wind speeds
        factor {float}              -- Multiplication factor (a in a*x^b)
        exp {float}                 -- Exponent (b in a*x^b)
    
    Returns:
        damage_ratio {float}        -- Ratio of building damage

    """
    damage_ratio = (storm_map * factor)**exp
    damage_ratio[np.isnan(damage_ratio)] = 0.
    damage_ratio[damage_ratio > 1.0] = 1.0

    return damage_ratio


def get_max_damage(country, max_damage_dict):
    """ Get maximum damage per km2 for given country, based on maximum damages
        of Huizinga et al., aggregated from building classes using ratio from
        Tiggeloven et al.

        TODO: insert citation to articles!

    Arguments:
        country {string}            -- Country name (as given in max damage file)
        max_damage_dict {string}    -- Dict with all max damage data
        gdp_year {int}              -- Year of exposure data to which GDP should be
                                       corrected

    Returns:
        max_damage {float}          -- Maximum building damage value per grid cell 
                                       for given country
    """

    max_damage = (0.2 * 0.75 * max_damage_dict[country]["RES"] + 
                  0.3 * 0.15 * max_damage_dict[country]["COM"] + 
                  0.3 * 0.10 * max_damage_dict[country]["IND"])

    # Convert m2 to km2
    max_damage *= 1e6

    return max_damage


def get_region_shape(region_code, region_level, region_path):
    """ Get shape of whole region for cropping maps

    Arguments:
        region_code {string}        -- Region to extract from map, must match 
                                       code of NUTS region on any level
        region_level {int}          -- Level (0, 1 or 2) of regions to extract
        region_path {string}        -- Path to file with NUTS regions

    Returns:
        geometry_union {Polygon}    -- Polygon of full region shape

    """

    region_shp = fiona.open(region_path, "r")
    level_name = "GID_" + str(region_level)
    geometries_region = [shape(region["geometry"]) for region in region_shp 
                          if region["properties"][level_name] == region_code]
    geometry_union = unary_union(geometries_region)

    return geometry_union


def crop_region(region_name, region_shape, crop_path, temp_dir):
    """ Crop maps based on polygon shape and saves file to temporary
        directory for further processing.

    Arguments:
        region_name {string}        -- Region name (full name)
        region_shape {polygon}      -- Polygon of region to crop
        crop_path {string}          -- Path to file to crop
        temp_dir {string}           -- Path to temporary directory

    Returns:
        cropped_map {float}         -- Raster of input file cropped to region
                                       of interest
    """

    # Read map to crop from file(s)
    crop_map = rasterio.open(crop_path)
    out_meta = crop_map.meta

    # Extract cropped map based on country shape
    try:
        cropped_map, transform = mask(crop_map, [region_shape], 
                                      all_touched=True, crop=True)

        # Save cropped file as GeoTIFF in temp directory
        cropped_path = (temp_dir + "/" + crop_path.split("/")[-1].split(".")[-2] 
                        + "_region_" + region_name + "_cropped.tif")
        out_meta.update({"height": cropped_map.shape[1],
                         "width": cropped_map.shape[2],
                         "transform": transform})
        with rasterio.open(cropped_path, "w", **out_meta) as dest:
            dest.write(cropped_map)

        return cropped_path
    
    except ValueError:
        # print("No overlap for " + crop_path.split("/")[-1])
        return None


def get_overlapping_storms(storms_overlap, storm, storm_dir, country, country_shape, temp_dir):
    """ Get storms that overlap with country for damage computation.

    Arguments:
        storms_overlap {list}           -- Current list of overlapping storms
        storm {string}                  -- Storm name
        storm_dir {string}              -- Directory with all storms
        country {string}                -- Country code
        country_shape {polygon}         -- Polygon of country shape
        temp_dir {string}               -- Directory for temporary files

    Returns:
        storms_overlap {list}           -- List of storms that overlap with country

    """

    # print("Computing overlap for " + str(storm))
    storm_path = storm_dir + "/" + storm + ".tif"
    save_dir = temp_dir + "/" + storm.split("/")[-1]
    country_storm_cropped = crop_region(country, country_shape, storm_path,
                                        save_dir)
    if country_storm_cropped:
        # print("Overlap for storm " + storm)
        storms_overlap.append(storm)

    for f in os.listdir(save_dir):
        if f.split("_")[0] == "STORM":
            os.remove(save_dir + "/" + f)

    return storms_overlap


def match_rasters(path1, path2):
    """ Get intersection of two raster files. 
        NOTE: raster files should have equal projection and resolution.

    Arguments:
        path1 {string}          -- Path to first raster
        path2 {string}          -- Path to second raster

    Returns:
        array1 {2D array}       -- Raster with first raster data 
                                   for overlapping area
        array2 {2D array}       -- Raster with second raster data 
                                   for overlapping area
    """
    
    # Load raster data, bands and geotransforms
    raster1 = gdal.Open(path1)
    raster2 = gdal.Open(path2)
    band1 = raster1.GetRasterBand(1)
    band2 = raster2.GetRasterBand(1)
    gt1 = raster1.GetGeoTransform()
    gt2 = raster2.GetGeoTransform()
    
    # Get bounding boxes and intersection
    box1 = [gt1[0], gt1[3], gt1[0] + (gt1[1] * raster1.RasterXSize), 
            gt1[3] + (gt1[5] * raster1.RasterYSize)]
    box2 = [gt2[0], gt2[3], gt2[0] + (gt2[1] * raster2.RasterXSize), 
            gt2[3] + (gt2[5] * raster2.RasterYSize)]
    intersect = [max(box1[0], box2[0]), min(box1[1], box2[1]), 
                 min(box1[2], box2[2]), max(box1[3], box2[3])]

    # Check if bounding boxes differ
    if box1 != box2:

        # Check if any overlap exists
        if (intersect[2] < intersect[0]) or (intersect[1] < intersect[3]):
            intersect = None
            return

        # Compute overlap
        else:
            left1 = int(round((intersect[0] - box1[0])/gt1[1]))
            top1 = int(round((intersect[1] - box1[1])/gt1[5]))
            col1 = int(round((intersect[2] - box1[0])/gt1[1])) - left1
            row1 = int(round((intersect[3] - box1[1])/gt1[5])) - top1
            
            left2 = int(round((intersect[0] - box2[0])/gt2[1]))
            top2 = int(round((intersect[1] - box2[1])/gt2[5]))
            col2 = int(round((intersect[2] - box2[0])/gt2[1])) - left2
            row2 = int(round((intersect[3] - box2[1])/gt2[5])) - top2
            
            array1 = band1.ReadAsArray(left1, top1, col1, row1)
            array2 = band2.ReadAsArray(left2, top2, col2, row2)

    else:
        array1 = band1.ReadAsArray()
        array2 = band2.ReadAsArray()
    
    return array1, array2


def compute_damage(exposure, storm, surface_area, max_damage, param1, param2, vul_func):
    """ Compute damage for given exposure and storm maps using vulnerability
        function parameterized by param1 and param2.

    Arguments:
        exposure {string}           -- Path to exposure raster
        storm {string}              -- Path to storm raster
        surface_area {string}       -- Path to land surface ratio file
        max_damage {float}          -- Maximum damage per grid cell

    Parameters:
        param1 {float}            -- Damage threshold parameter for 
                                       vulnerability function
        param2 {float}              -- Half-damage parameter for vulnerability 
                                       function
        vul_func {string}           -- Type of vulnerability function
    
    Returns:
        damage {2D array}           -- Raster of resulting damage
    """

    # Save GeoTIFF meta data
    out_meta = rasterio.open(storm).meta

    # Get overlap for storm and exposure input
    result = match_rasters(exposure, storm)

    # Check if maps overlap and compute damage
    if result and not np.isnan(result[1]).all():
        exposure_area, storm_area = result
        exposure_area[exposure_area == -200] = np.nan

        if vul_func == "sigmoid":
            damage_ratio = sigmoid_vulnerability_function(storm_area, param1, param2)
        elif vul_func == "power":
            damage_ratio = power_vulnerability_function(storm_area, param1, param2)

        damage = (1/100. * exposure_area) * damage_ratio * max_damage

        out_meta.update({"dtype": "float64",
                         "height": damage.shape[0],
                         "width": damage.shape[1]})

        # Save intermediate damage result for surface area correction
        name_list = storm.split("/")[-1].split("_")
        damage_path = (os.path.abspath(os.path.join(exposure, os.pardir))
                       + "/Damage_storm_" + name_list[1] + "_" + name_list[2] 
                       + "_region_" + name_list[-2] + ".tif")
        with rasterio.open(damage_path, "w", **out_meta) as dest:
            dest.write(np.expand_dims(damage, axis=0))

        # Correct for not using equal area projection
        damage, surface = match_rasters(damage_path, surface_area)
        surface[surface == -9999] = np.nan
        damage = damage * surface

        out_meta.update({"dtype": "float64",
                         "height": damage.shape[0],
                         "width": damage.shape[1]})

        # # Save final damage file
        # with rasterio.open(damage_path, "w", **out_meta) as dest:
        #     dest.write(np.expand_dims(damage, axis=0))

    else:
        print("No overlap found, skipping these maps.")
        damage = None

    os.remove(exposure)
    os.remove(storm)
    os.remove(surface_area)

    return damage


def damage_per_region(damage_dict, storm, storm_dir, country, region, region_shape, exposure_path, surface_path, max_damage_dict, temp_dir, param1, param2, vul_func):
    """ Get cropped files per region and, compute damage for given storm and region
        and return updated dict.

    Arguments:
        damage_dict {dict}              -- Dict for saving damage per storm
        storm {string}                  -- Name of storm file
        storm_dir {string}              -- Directory containing all storms
        country {string}                -- Country of region of interest (ISO code)
        region {string}                 -- Region name (full name)
        region_shape {polygon}          -- Polygon of region to crop
        exposure_path {string}          -- Path to exposure file
        surface_path {string}           -- Path to file for surface correction
        max_damage_dict {dict}          -- Dict with maximum damage per country (ISO)
        temp_dir {string}               -- Directory for saving temp files

    Parameters:
        param1 {float}                  -- Damage threshold parameter for 
                                           vulnerability function
        param2 {float}                  -- Half-damage parameter for vulnerability 
                                           function
        vul_func {string}               -- Type of vulnerability function

    Returns:
        damage_dict {dict}              -- New damage dict

    """

    # Create path for saving temporary data
    storm_path = storm_dir + "/" + storm + ".tif"
    save_dir = temp_dir + "/" + storm.split("/")[-1]

    # If storm overlaps with region, compute damage
    storm_cropped_path = crop_region(region, region_shape, storm_path, save_dir)

    if storm_cropped_path:
        print("Computing damage for " + storm_path.split("/")[-1]
              + " and region " + region)
        print("---------------------------------------------------------")
        exposure_cropped_path = crop_region(region, region_shape, 
                                            exposure_path, save_dir)
        surface_cropped_path = crop_region(region, region_shape, 
                                           surface_path, save_dir)

        max_damage = get_max_damage(country, max_damage_dict)
        damage = compute_damage(exposure_cropped_path, storm_cropped_path, 
                                surface_cropped_path, max_damage, 
                                param1, param2, vul_func)
        damage_dict[storm] = np.nansum(damage)

    else:
        # print("No overlap for storm " + storm)
        damage_dict[storm] = 0.

    return damage_dict


