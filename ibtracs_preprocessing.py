# -*- coding: utf-8 -*-

"""
Created on Thu Mar 4, 2021

@author: Liz Verbeek

This script is part of the TC risk model developed as part of a Master Thesis 
for the Master's Programme Computational Science at the University of Amsterdam, 
see https://github.com/lizverbeek/global_TC_risk_model .

Code in this script has been adapted from preprocessing.py as used in 
Bloemendaal et al, Estimation of global tropical cyclone wind probabilities 
using the STORM dataset. Sci Data 7, 377 (2020). https://doi.org/10.1038/s41597-020-00720-x

This script contains all preprocessing steps for selecting suitable TC tracks
from the IBTrACS dataset.

"""

import os
import sys
import numpy as np

dir_path = os.path.dirname(os.path.realpath(sys.argv[0]))
__location__ = os.path.realpath(os.path.join(os.getcwd(), 
                                             os.path.dirname(__file__)))


def interpolate(dataset):
    """
    Interpolate dataset to contain values for each 3-hourly timestep.

    Arguments: 
        dataset: the respective dataset

    Returns:
        dataset: the interpolated dataset
    
    """
    
    if (np.any(np.isnan(dataset))==True and 
        len([x for x,v in enumerate(dataset) if np.isnan(v)==False])>1):
        
        # Get indices with values
        ind = [x for x,v in enumerate(dataset) if np.isnan(v)==False]
        # Get indices with no values
        ind1 = [x + ind[0] for x,v in enumerate(dataset[ind[0]:ind[-1]]) 
                if np.isnan(v)==True]
        
        val = [v for v in dataset if np.isnan(v)==False]
        
        # Interpolate between known values
        if len(ind1) > 0:
            interlist=np.interp(ind1,ind,val)
            for ii,jj in zip(ind1,range(len(ind1))):
                dataset[ii]=interlist[jj]
    
    return dataset


def fill_missing_values(dataset):
    """
    Fill missing values by extending first and last known value.

    Input: 
        dataset: the respective dataset
    Output:
        dataset: the dataset extended with filled values
    
    """

    ind = np.where(~np.isnan(dataset))[0]

    # Check if more than one value is available
    if len(ind) > 1:
        missing = np.where(np.isnan(dataset))[0]

        # Check if front has missing values
        if missing[0] == 0:
            for i in missing[np.where(missing < ind[0])]:
                dataset[i] = dataset[ind[0]]

        # Check if end has missing values
        if missing[-1] == (len(dataset) - 1):
            for i in missing[np.where(missing > ind[-1])]:
                dataset[i] = dataset[ind[-1]]
    else:
        # Remove data from further use
        dataset = []

    return dataset


def check_timelist(tlist):
    """
    Check whether the consecutive time steps are 3 hours apart
    Input:
        tlist: list of time steps
    
    """

    for ii in range(1,len(tlist)):
        if tlist[ii]-tlist[ii-1]!=0.125:
            print(tlist)

    return


def convert_wind_speed(wind, agency):
    """
    Convert IBTrACS wind speed to 10-min sustained wind speed.
    
    Input: 
        wind: wind speed 
        agency: name of agency
    Output:
        wind_conv: converted wind    
    """
    
    # Convert for agencies that report 1-minute average
    if (agency == "hurdat_epa" or agency=="hurdat_atl" 
        or agency=="newdelhi" or agency=="atcf"):
        wind = 0.88 * wind
        
    return wind


def extract_data(data, output_path):
    """
    Extract different variables from IBTrACS dataset.
    Input:
        *data*: dataset (IBTrACS)
    Output: 
        *LATLIST_INTERP.npy*: interpolated values of latitude, where each entry in the dictionary stands for one TC
        *LONLIST_INTERP.npy*: interpolated values of longitude (0-360 deg)
        *WINDLIST_INTERP.npy*: interpolated values of wind (m/s)
        *PRESLIST_INTERP.npy*: interpolated values of pressure (hPa)
        *RMAXLIST_INTERP.npy*: interpolated values of Rmax (km)
        *BASINLIST_INTERP.npy*: Basin of TC genesis
        *YEARLIST_INTERP.npy*: Year of TC genesis
    """
    
    basin = data.basin.values
    years = data.season.values
    names = data.name.values
    wind = data.wmo_wind.values
    wind = wind * 0.51444444                            # Convert knots to m/s
    pres = data.wmo_pres.values
    time = data.time.values
    latitude = data.lat.values
    longitude = data.lon.values
    rmax_usa = data.usa_rmw.values * 1.85200            # Convert nm to km
    rmax_reunion = data.reunion_rmw.values * 1.85200    # Convert nm to km
    rmax_bom = data.bom_rmw.values * 1.85200            # Convert nm to km
    wmo_agency = data.wmo_agency.values
    nature = data.nature.values
    landfall = data.landfall.values
    
    """Create a npy list for each of the items"""    
    latlist = {i:[] for i in range(len(years))}
    lonlist = {i:[] for i in range(len(years))}
    timelist = {i:[] for i in range(len(years))}
    windlist = {i:[] for i in range(len(years))}
    preslist = {i:[] for i in range(len(years))}
    basinlist = {i:[] for i in range(len(years))}
    rmaxlist = {i:[] for i in range(len(years))}
    yearlist = {i:[] for i in range(len(years))}
    namelist = {i:[] for i in range(len(years))}
    
    # Loop through all storms
    for i in range(len(years)):

        # Include only storms before 2018 that have made landfall
        if (years[i] < 2018 
            and np.all(np.isnan(landfall[i])) == False 
            and np.nanmin(landfall[i]) == 0):

            # Check if data on agency is available for at least one timestep.
            agency_idx = [x for x,v in enumerate(wmo_agency[i]) if len(v)>1.]
            if len(agency_idx) > 0:
                
                # Convert wind speed to 10-min maximum sustained average.
                agency = wmo_agency[i][agency_idx[0]].decode("utf-8")
                wind_conv = convert_wind_speed(wind[i], agency)

                # Check if wind speed values are available at at least one timestep
                ind = [x for x,v in enumerate(wind_conv)]

                # Exclude SA basin and check availability of wind speed values.
                if (np.all(np.isnan(wind_conv)) == False and len(ind) > 0
                    and basin[i][ind[0]].decode("utf-8")!="SA"):

                    j0 = ind[0]
                    j1 = ind[-1]
                    
                    # Save basin (string) and year data
                    basinlist[i].append(basin[i][ind[0]].decode("utf-8"))
                    yearlist[i].append(years[i])
                    namelist[i].append(names[i].decode("utf-8"))
                    
                    # Get 3-hourly timesteps with wind speed values
                    time_idx = [j0+x for x,v in enumerate(time[i][j0:j1+1]) 
                                if round(v,3)%0.125==0.]
                    new_list = np.intersect1d(ind, time_idx)

                    # Check if more than one timestep has wind speed value
                    if len(new_list) > 1.:
                        # and np.all(np.isnan(rmax[i])) == False:
                        n0 = time_idx.index(new_list[0])
                        n1 = time_idx.index(new_list[-1])
                       
                        new_time = time_idx[n0:n1+1]
                        
                        for j_idx in range(len(new_time)):
                            j = new_time[j_idx]        
                            latlist[i].append(latitude[i][j])

                            if longitude[i][j] >= 180.:
                                longitude[i][j] -= 360.
                            
                            lonlist[i].append(longitude[i][j])
                            timelist[i].append(round(time[i][j],3))
                            windlist[i].append(wind_conv[j])
                            preslist[i].append(pres[i][j])

                            # Check for available Rmax data
                            if np.all(np.isnan(rmax_usa[i])) == False:
                                rmaxlist[i].append(rmax_usa[i][j])
                            elif np.all(np.isnan(rmax_reunion[i])) == False:
                                rmaxlist[i].append(rmax_reunion[i][j])
                            elif np.all(np.isnan(rmax_bom[i])) == False:
                                rmaxlist[i].append(rmax_bom[i][j])

                        check_timelist(timelist[i])
    
    # Fill missing values where possible.
    for i in range(len(latlist)):
        if (len(latlist[i]) > 0 and len(windlist[i]) > 0 and 
            len(preslist[i]) > 0 and len(rmaxlist[i])):

            if np.isnan(windlist[i][0]) or np.isnan(windlist[i][-1]):
                windlist[i] = fill_missing_values(windlist[i])
            if np.isnan(preslist[i][0]) or np.isnan(preslist[i][-1]):
                preslist[i] = fill_missing_values(preslist[i])
            if np.isnan(rmaxlist[i][0]) or np.isnan(rmaxlist[i][-1]):
                rmaxlist[i] = fill_missing_values(rmaxlist[i])
            
            latlist[i] = interpolate(latlist[i])
            lonlist[i] = interpolate(lonlist[i])
            windlist[i] = interpolate(windlist[i])
            preslist[i] = interpolate(preslist[i])
            rmaxlist[i] = interpolate(rmaxlist[i])
    
    # Get indices of storms for which at least one variable is missing
    missing_values = set()
    missing_values.update({index for index, value in latlist.items() if not value})
    missing_values.update({index for index, value in lonlist.items() if not value})
    missing_values.update({index for index, value in timelist.items() if not value})
    missing_values.update({index for index, value in windlist.items() if not value})
    missing_values.update({index for index, value in preslist.items() if not value})
    missing_values.update({index for index, value in rmaxlist.items() if not value})
    missing_values.update({index for index, value in basinlist.items() if not value})
    missing_values.update({index for index, value in yearlist.items() if not value})
    missing_values.update({index for index, value in namelist.items() if not value})

    # Remove storms with missing values for any of the variables
    for i in missing_values:
        del latlist[i]
        del lonlist[i]
        del timelist[i]
        del windlist[i]
        del preslist[i]
        del rmaxlist[i]
        del basinlist[i]
        del yearlist[i]
        del namelist[i]

    print("Storms removed: " + str(len(missing_values)))
    print("Length of resulting dataset: " + str(len(yearlist)))
    print("Saving datasets at " + output_path)

    # Save preprocessed data.
    np.save(os.path.join(output_path, "LATLIST_INTERP.npy"), latlist)
    np.save(os.path.join(output_path, "LONLIST_INTERP.npy"), lonlist)
    np.save(os.path.join(output_path, "TIMELIST_INTERP.npy"), timelist)
    np.save(os.path.join(output_path, "WINDLIST_INTERP.npy"), windlist)
    np.save(os.path.join(output_path, "PRESLIST_INTERP.npy"), preslist)
    np.save(os.path.join(output_path, "RMAXLIST_INTERP.npy"), rmaxlist)
    np.save(os.path.join(output_path, "BASINLIST_INTERP.npy"), basinlist)
    np.save(os.path.join(output_path, "YEARLIST_INTERP.npy"), yearlist)
    np.save(os.path.join(output_path, "NAMELIST_INTERP.npy"), namelist)
