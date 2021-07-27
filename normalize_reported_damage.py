#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Created on Thu May 12, 2021

@author: Liz Verbeek
"""

import os
import rasterio
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


# ============================================================================== 
# Get paths
# ============================================================================== 
current_dir = os.getcwd()
top_dir = os.path.abspath(os.path.join(current_dir, os.pardir))
exposure_dir = top_dir + "/data/GHSL_data/World/EPSG_4326"

match_dir = top_dir + "/data/storm_data/Storms_with_damage_data"
total_damage_file = match_dir + "/storm_damage_matches_ALL.csv"
gdp_file = top_dir + "/data/max_damage_data/GDP_current_USD.csv"

storm_dir = match_dir + "/EPSG_4326/World_Res"

total_damage_df = pd.read_csv(total_damage_file, keep_default_na=False)
gdp_df = pd.read_csv(gdp_file, header=2)

ratio_2000 = []
ratio_2015 = []
for country, year in zip(total_damage_df["ISO"], total_damage_df["Year"]):
	gdp_old = gdp_df[gdp_df["Country Code"] == country][str(year)].values[0]
	gdp_2000 = gdp_df[gdp_df["Country Code"] == country]["2000"].values[0]
	gdp_2015 = gdp_df[gdp_df["Country Code"] == country]["2015"].values[0]
	if np.isnan(gdp_old):
		ratio_2000.append((gdp_df[str(year)]/gdp_df["2000"]).mean())
		ratio_2015.append((gdp_df[str(year)]/gdp_df["2015"]).mean())
	else:
		ratio_2000.append(gdp_2000/gdp_old)
		ratio_2015.append(gdp_2015/gdp_old)

total_damage_df["Damage year 2000"] = total_damage_df["Damage"] * ratio_2000
total_damage_df["Damage year 2015"] = total_damage_df["Damage"] * ratio_2015
total_damage_df = total_damage_df.set_index("Nr")

total_damage_df.to_csv(match_dir + "/storm_damage_matches_ALL.csv")

# plt.bar(total_damage_df["Country"])
# plt.xticks(rotation="vertical")
# plt.show()

# # Too big
# src = rasterio.open(exposure_dir + 
# 						     "/GHS_BUILT_LDS2000_GLOBE_R2018A_4326_1K_V2_0.tif")
# exposure_map = src.read(1)
# plt.imshow(exposure_map)
# plt.show()