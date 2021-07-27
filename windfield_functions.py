#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

@author: Liz Verbeek

TODO: add adaptation from Bloemendaal et al., + reference to code.
"""

# """
# This script is part of Bloemendaal et al, Estimation of global tropical cyclone wind probabilities using the STORM dataset (in review)
# The script has been developed by Nadia Bloemendaal, Job Dullaart and Sanne Muis. 
# This script is the master program complementary to holland_model.py. The methodology is heavily inspired by 

# Lin, N., and Chavas, D. ( 2012), On hurricane parametric wind and applications in storm surge modeling, 
# J. Geophys. Res., 117, D09120, doi:10.1029/2011JD017126.


# Copyright (C) 2020 Nadia Bloemendaal. All versions released under GNU General Public License v3.0.
# """

import os
import math
import numpy as np

from scipy import spatial
from osgeo import osr


def Basins_WMO(basin, tc_radius=1000.):
    """ 
    Get boundary coordinates for each of the ocean basins. 
        
    Arguments:
        basin {string}      -- Ocean basin

    Returns:
        lat0 {int}          -- Lower bound latitude
        lat1 {int}          -- Upper bound latitude
        lon0 {int}          -- Left bound longitude
        lon1 {int}          -- Right bound longitude

    """
    if basin == 'EP': # Eastern Pacific
        lat0, lat1, lon0, lon1 = 5, 60, -180, -75
    # if basin == 'EP': # Eastern Pacific
    #     lat0, lat1, lon0, lon1 = 5, 60, 180, 285
    elif basin == 'NA': #North Atlantic
        lat0, lat1, lon0, lon1 = 5, 60, -105, -1
    # elif basin == 'NA': #North Atlantic
    #     lat0, lat1, lon0, lon1 = 5, 60, 255, 359
    elif basin == 'NI': #North Indian
        lat0, lat1, lon0, lon1 = 5, 60, 30, 100
    elif basin == 'SI': #South Indian
        lat0, lat1, lon0, lon1 = -60, -5, 10, 135
    elif basin == 'SP': #South Pacific
        lat0, lat1, lon0, lon1 = -60, -5, 135, 240
    elif basin == 'WP': #Western Pacific
        lat0, lat1, lon0, lon1 = 5, 60, 100, 180
    else:
        print("Basin " + basin + " not recognized")
    
    # Add radius in degrees to boundaries so that whole storm fits within basin
    max_distance = tc_radius/111.
    lat0 -= max_distance
    lat1 += max_distance
    lon0 -= max_distance
    lon1 += max_distance

    return lat0, lat1, lon0, lon1


def get_basin_points(basin, res=0.1):
    """
    Get list of (lat, lon), empty wind field with corresponding indices
    and cKDTree of the given points for given basin.

    Arguments:
        basin {string}          -- Ocean basin
        res {float}             -- Resolution for 2D wind fields
        
    Returns:
        points {dict}           -- Dict with (lat,lon) pairs
        wind_field {dict}       -- Empty dict with same indices as points dict
        tree {cKDTree}          -- Search tree of given points 

    """
    lat0, lat1, lon0, lon1 = Basins_WMO(basin)

    if lat0 > 0:
        start = lat0 + 0.5*res
        end = lat1 + 0.5*res
        n_lat = int((end-start)*(1/res))
        latspace = np.linspace(start, end, n_lat)
    else:
        start = lat0 - 0.5*res
        end = lat1 - 0.5*res
        n_lat = int((end-start)*(1/res))
        latspace = np.linspace(start, end, n_lat)

    start = lon0 + 0.5*res
    end = lon1 + 0.5*res
    n_lon = int((end-start)*(1/res))
    lonspace = np.linspace(start, end, n_lon)

    # Create list of (lat, lon) coords and corresponding wind field
    points = [(i,j) for i in latspace for j in lonspace]
    wind_field = {i:[] for i in range(len(points))}
    tree = spatial.cKDTree(points)

    return latspace, lonspace, points, tree


def windfield_to_geotiff(lats, lons, wind_field, crs, file_path, res=0.1):
    """
    Save 2D numpy array of windfield data with specified coordinates
    as a GeoTIFF file at given path.

    Arguments:
        lats {1d array}             -- Raster latitudes
        lons {1d array}             -- Raster longitudes
        wind_field {1d array}       -- Flattened 1D array with wind speed data
        crs {int}                   -- EPSG number for GeoTIFF projection
        file_path {string}          -- Path where data should be saved
        res {float}                 -- Resolution of wind fields
    
    """

    rows = len(lats)
    cols = len(lons)
    lat_min = np.min(lats)
    lat_max = np.max(lats)
    lon_min = np.min(lons)
    lon_max = np.max(lons)
    wind_field = np.reshape(wind_field, (rows, cols))

    # Flip array (needed for correct projection)
    wind_field = np.flipud(wind_field)

    driver = gdal.GetDriverByName("GTiff")
    source = osr.SpatialReference()
    source.ImportFromEPSG(crs)
    output_proj = source.ExportToPrettyWkt()

    output = driver.Create(file_path, cols, rows, 1, gdal.GDT_Float64)
    output.SetProjection(output_proj)
    output.SetGeoTransform((lon_min - 0.5 * res, res, 0,
                            lat_max + 0.5 * res, 0, -res))
    output.GetRasterBand(1).WriteArray(wind_field)

    return


def compute_windfield(storm, alpha=0.55, beta_bg=20., SWRF=0.85, CF=0.915,
                      tc_radius=1000., Patm=101325., n_cols=36, n_rows=1000, 
                      crs=4326, res=0.1):
    """
    Compute the 2D wind field from storm track data using the Holland model.

    Arguments:
        storm {dict}                -- Dict containing all storm variables

    Constants (with reference to literature):
        alpha               -- Deceleration of surface background wind
                               (Lin & Chavas 2012)
        beta_bg             -- Angle of background wind flow
                               (Lin & Chavas 2012)
        SWRF                -- Empirical Surface Wind Reduction Factor
                               (Powell et al 2005)
        CF                  -- Wind conversion factor from 1-min to 10-min average
                               (Harper et al (2012))
        tc_radius           -- Radius of the tropical cyclone, in km
        Patm                -- Ambient pressure
        n_cols              -- Number of gridpoints in angular direction
        n_rows              -- Number of gridpoints in radial direction
        crs                 -- EPSG for wind field projection

    Returns:
        wind_field {2D array}       -- 2D raster with maximum sustained wind speeds
    
    """

    max_distance = tc_radius/111.

    # Get coordinates, search tree and wind field for storm basin
    basin = storm['basin']
    latspace, lonspace, points, tree = get_basin_points(basin)
    wind_field = np.empty((len(latspace)*len(lonspace)))
    wind_field[:] = np.nan

    # Extract values per timestep
    shadowlist = {idx:[] for idx in range(len(points))}   # TODO: copy wind_field? Is the same?
    for j in range(1, len(storm['lat'])):

        lat0, lat1 = storm['lat'][j-1], storm['lat'][j]
        lon0, lon1 = storm['lon'][j-1], storm['lon'][j]
        t0, t1 = storm['time'][j-1], storm['time'][j]
        U10, Rmax, P = storm['wind'][j], storm['rmax'][j], storm['pres'][j]

        # Get time difference in seconds
        # (timesteps in IBTrACS are given as fraction of a day).
        dt = (t1 - t0) * 3600 * 24.

        # Generate the seperate list of coastal points that are in the spyderweb
        distances, indices = tree.query((lat1, lon1), k=len(points), p=2, 
                                         distance_upper_bound=max_distance)
        points_to_save = [points[indices[k]] for k in range(len(distances)) 
                          if distances[k] < max_distance]

        # Spyderweb step 1: Generate the spyderweb mesh      
        rlist, thetalist, xlist, ylist = hm.Generate_Spyderweb_mesh(n_cols, n_rows, 
                                                                    tc_radius, lat0)
        lats, lons = hm.Generate_Spyderweb_lonlat_coordinates(xlist, ylist,
                                                              lat1, lon1)

        # Spyderweb step 2: Calculate the background wind
        [bg, ubg, vbg] = hm.Compute_background_flow(lon0, lat0, lon1, lat1, dt)

        # Spyderweb step 3: Subtract background flow from U10 (TC's 10-meter wind speed)
        # For this, first convert U10 to surface level using the SWRF-constant
        # Next, subtract a fraction alpha of the background flow.
        Usurf = (U10/SWRF) - (bg * alpha)         # 1-min max sustained surface winds
        P_mesh = np.zeros((xlist.shape))
        Pdrop_mesh = np.zeros((xlist.shape))
        up = np.zeros((xlist.shape))
        vp = np.zeros((xlist.shape))

        # Spyderweb step 4: Calculate wind and pressure profile using the Holland model
        for l in range(1, n_rows):
            r = rlist[0][l]
            Vs, Ps = hm.Holland_model(lat1, P, Usurf, Rmax, r)
            Vs = Vs * SWRF                    # Convert back to 10-min wind speed
            P_mesh[:, l].fill(Ps/100.)                  # in Pa
            Pdrop_mesh[:, l].fill((Patm - Ps)/100.)     # in Pa
            beta = hm.Inflowangle(r, Rmax, lat0)
          
            for k in range(0, n_cols):
                ubp = alpha * (ubg * math.cos(math.radians(beta_bg)) - 
                                              np.sign(lat0) * vbg * 
                                              math.sin(math.radians(beta_bg)))
                vbp = alpha * (vbg * math.cos(math.radians(beta_bg)) + 
                                              np.sign(lat0) * ubg * 
                                              math.sin(math.radians(beta_bg)))
                up[k,l] = -Vs * math.sin(thetalist[:,0][k] + beta) + ubp
                vp[k,l] = -Vs * math.cos(thetalist[:,0][k] + beta) + vbp

        u10 = CF * up
        v10 = CF * vp
        windfield = np.sqrt(u10**2. + v10**2.)

        spy_points = []
        wind_points = []
        for k in range(n_cols):
            for l in range(n_rows):
                spy_points.append((lats[k,l], lons[k,l]))
                wind_points.append(windfield[k,l])

        tree2 = spatial.cKDTree(spy_points)

        # Overlay the spyderweb grid with the regular grid
        for (lat, lon), idx in zip(points_to_save, range(len(points_to_save))):
            local_dist, local_ind = tree2.query((lat,lon), k=1, p=2, 
                                                 distance_upper_bound=max_distance)    
            shadowlist[indices[idx]].append(wind_points[local_ind])


    for m in range(len(shadowlist)):
        if len(shadowlist[m]) > 0.:
            # if np.max(shadowlist[m]) >= 18.: (Tropical cyclones only)
            # Create wind field (flattened) raster of maximum wind speeds
            wind_field[m] = np.max(shadowlist[m])


    if not np.all(np.isnan(wind_field)):
        # Save wind field as GeoTIFF
        windfield_to_geotiff(latspace, lonspace, wind_field, crs, output_dir + 
                             "/STORM_WIND_SPEEDS_No" + str(storm['nr']) + "_" + 
                             basin + ".tif")
    else:
        print("WARNING! STORM No " + str(storm['nr']) + " HAS ONLY NaN VALUES!")
    print("STORM No " + str(storm['nr']) + " DONE")
    print("--------------------------------------")

    return wind_field

