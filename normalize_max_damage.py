#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Created on Thu May 6, 2021

@author: Liz Verbeek

This script is part of the TC risk model developed as part of a Master Thesis 
for the Master's Programme Computational Science at the University of Amsterdam, 
see https://github.com/lizverbeek/global_TC_risk_model .

In this script, maximum damage values per country are normalized to match
the year of the exposure data.

Maximum damage values have been used from 
Huizinga J, de Moel H, Szewczyk W. Global flood depth-damage functions: 
Methodology and the database with guidelines. Joint Research Centre (Seville site); 
2017 Apr. https://doi.org/10.2760/16510

"""

import os
import csv
import fiona
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy import stats

import statsmodels.api as sm
import pylab


# ============================================================================== 
# Get paths 
# ============================================================================== 
current_dir = os.getcwd()
top_dir = os.path.abspath(os.path.join(current_dir, os.pardir))
max_damage_dir = top_dir + "/data/max_damage_data"
region_file = top_dir + "/data/NUTS_WORLD/gadm36_2.shp"

# ============================================================================== 
# Read countrie codes and names from NUTS file
# ============================================================================== 
NUTS_dict = {}
region_shp = fiona.open(region_file, "r")
for region in region_shp:
    code = region["properties"]["GID_0"]
    name = region["properties"]["NAME_0"]
    NUTS_dict[code] = name

# ============================================================================== 
# Read max damage data
# ============================================================================== 
max_damage_df = pd.read_csv(max_damage_dir + "/max_damage_total.csv", 
                            usecols=["Country", "ISO A3 code", 
                                     "RES", "COM", "IND"])
max_damage_df = max_damage_df.rename(columns={"Country": "Country Name"})
max_damage_df = max_damage_df.rename(columns={"ISO A3 code": "Country"})

# ============================================================================== 
# Convert max damage values from Euro to USD
# ============================================================================== 
max_damage_df["RES"] = max_damage_df["RES"] * 1/0.77
max_damage_df["COM"] = max_damage_df["COM"] * 1/0.77
max_damage_df["IND"] = max_damage_df["IND"] * 1/0.77

# ============================================================================== 
# Read GDP data
# ============================================================================== 
current_gdp_df = pd.read_csv(max_damage_dir + "/GDP_current_USD.csv", 
                             usecols=["Country Code", "2000", "2010", "2015"])
current_gdp_df = current_gdp_df.rename(columns={"Country Code": "Country"})

# ============================================================================== 
# For all countries in NUTS file: get max damage data if available, 
# otherwise fill max damage data with estimate.
# ============================================================================== 
gdp_USA_2010 = current_gdp_df[current_gdp_df["Country"] == "USA"]["2010"].values[0]
max_damage_USA = max_damage_df[max_damage_df["Country"] == "USA"]
for country in NUTS_dict.keys():
    if country not in max_damage_df["Country"].values:
        gdp_country_2010 = current_gdp_df[current_gdp_df["Country"] == country]["2010"].values
        if len(gdp_country_2010) > 0 and not np.isnan(gdp_country_2010[0]):
            ratio = gdp_country_2010[0]/gdp_USA_2010
            new_row = {"Country Name": NUTS_dict[country], "Country": country, 
                       "RES": max_damage_USA["RES"].values[0]*ratio, 
                       "COM": max_damage_USA["COM"].values[0]*ratio,
                       "IND": max_damage_USA["IND"].values[0]*ratio}
            max_damage_df = max_damage_df.append(new_row, ignore_index=True)
        else:
            new_row = {"Country Name": NUTS_dict[country], "Country": country, 
                       "RES": 0, "COM": 0, "IND": 0}

# ============================================================================== 
# Normalize maximum damages for years 2000 and 2015
# ============================================================================== 
# Compute ratios and add to new dataframe
current_gdp_df["Ratio GDP 2015"] = current_gdp_df["2015"]/current_gdp_df["2010"]
current_gdp_df["Ratio GDP 2000"] = current_gdp_df["2000"]/current_gdp_df["2010"]
total_df_2000 = pd.merge(max_damage_df, current_gdp_df[["Country" , "Ratio GDP 2000"]])
total_df_2015 = pd.merge(max_damage_df, current_gdp_df[["Country" , "Ratio GDP 2015"]])


# Fill missing values with mean ratio
mean_2000 = total_df_2000["Ratio GDP 2000"].mean()
mean_2015 = total_df_2015["Ratio GDP 2015"].mean()
total_df_2000["Ratio GDP 2000"] = total_df_2000["Ratio GDP 2000"].fillna(mean_2000)
total_df_2015["Ratio GDP 2015"] = total_df_2015["Ratio GDP 2015"].fillna(mean_2015)

# Compute new maximum damage values
max_damage_df = pd.merge(max_damage_df, total_df_2015[["Country", "Ratio GDP 2015"]], 
                         how="outer")
max_damage_df = pd.merge(max_damage_df, total_df_2000[["Country", "Ratio GDP 2000"]], 
                         how="outer")
max_damage_df["RES 2015"] = max_damage_df["RES"] * total_df_2015["Ratio GDP 2015"]
max_damage_df["COM 2015"] = max_damage_df["COM"] * total_df_2015["Ratio GDP 2015"]
max_damage_df["IND 2015"] = max_damage_df["IND"] * total_df_2015["Ratio GDP 2015"]
max_damage_df["RES 2000"] = max_damage_df["RES"] * total_df_2000["Ratio GDP 2000"]
max_damage_df["COM 2000"] = max_damage_df["COM"] * total_df_2000["Ratio GDP 2000"]
max_damage_df["IND 2000"] = max_damage_df["IND"] * total_df_2000["Ratio GDP 2000"]
max_damage_df = max_damage_df.rename(columns={"RES": "RES 2010", "COM": "COM 2010",
                                              "IND": "IND 2010"})

# Save new max damage values
max_damage_df = max_damage_df.drop(["Ratio GDP 2000", "Ratio GDP 2015"],  axis=1)
max_damage_df = max_damage_df.set_index(["Country"])
max_damage_df.to_csv(max_damage_dir + "/max_damage_GDP_corrected.csv")