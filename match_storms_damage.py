#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Created on Mon Apr 19, 2021

@author: Liz Verbeek

This script is part of the TC risk model developed as part of a Master Thesis 
for the Master's Programme Computational Science at the University of Amsterdam, 
see https://github.com/lizverbeek/global_TC_risk_model .

In this script, damage events from the EM-DAT dataset are matched to TC events
from the IBTrACS dataset.

"""

import os
import shutil
import datetime
import numpy as np
import pandas as pd


def preprocess_emdat_names(emdat_df):
	""" Removes all parts of storm names in the EM-DAT dataset that are 
		irrelevant for name matching. This function should be checked and
		possibly adapted when new data is available. Also note that the EM-DAT
		dataset storm names contain many spelling mistakes that have to be
		removed by hand as well.

	Arguments:
		emdat_df {dataframe}	-- Pandas dataframe with EM-DAT data

	Returns:
		emdat_df {dataframe}	-- Preprocessed EM-DAT dataframe

	"""

	# Remove irrelevant parts of name
	emdat_df["Event Name"] = emdat_df["Event Name"].str.lower()
	emdat_df["Event Name"] = emdat_df["Event Name"].str.replace("sorm", "storm")
	emdat_df["Event Name"] = emdat_df["Event Name"].str.replace("tropical ", "")
	emdat_df["Event Name"] = emdat_df["Event Name"].str.replace("strom", "storm")
	emdat_df["Event Name"] = emdat_df["Event Name"].str.replace("tropial", "tropical")
	emdat_df["Event Name"] = emdat_df["Event Name"].str.replace("tropcal", "tropical")
	emdat_df["Event Name"] = emdat_df["Event Name"].str.replace("topical", "tropical")
	emdat_df["Event Name"] = emdat_df["Event Name"].str.replace("cylone", "cyclone")
	emdat_df["Event Name"] = emdat_df["Event Name"].str.replace("tyhoon", "typhoon")
	emdat_df["Event Name"] = emdat_df["Event Name"].str.replace("typhhon", "typhoon")

	emdat_df["Event Name"] = emdat_df["Event Name"].str.replace("tropical ", "")
	emdat_df["Event Name"] = emdat_df["Event Name"].str.replace("cyclone ", "")
	emdat_df["Event Name"] = emdat_df["Event Name"].str.replace("cyclones ", "")
	emdat_df["Event Name"] = emdat_df["Event Name"].str.replace("storm ", "")
	emdat_df["Event Name"] = emdat_df["Event Name"].str.replace("hurricane ", "")
	emdat_df["Event Name"] = emdat_df["Event Name"].str.replace("depression ", "")
	emdat_df["Event Name"] = emdat_df["Event Name"].str.replace("depression", "")
	emdat_df["Event Name"] = emdat_df["Event Name"].str.replace("typhoon ", "")
	emdat_df["Event Name"] = emdat_df["Event Name"].str.replace("winter ", "")
	emdat_df["Event Name"] = emdat_df["Event Name"].str.replace("(", "")
	emdat_df["Event Name"] = emdat_df["Event Name"].str.replace(")", "")
	emdat_df["Event Name"] = emdat_df["Event Name"].str.replace("'", "")
	emdat_df["Event Name"] = emdat_df["Event Name"].str.replace('"', "")

	return emdat_df


def match_storm_names(emdat_df, ibtracs_df):
	""" Get storms that match by name and year (some names have been reused)
		and save these storms in new directory for later use.
	
	Arguments:
		emdat_df {dataframe}		-- Pandas dataframe with EM-DAT data
		ibtracs_df {dataframe}		-- Pandas dataframe with (preprocessed) 
									   IBTrACS data

	Returns:
		ibtracs_match_df {dataframe}		-- Dataframe with matching storms
		emdat_match_df {dataframe}			-- Dataframe with matching damage data
		ibtracs_no_match_df {dataframe}		-- Dataframe with rest of storms
		emdat_no_match_df {dataframe}		-- Dataframe with rest of damage data

	"""

	ibtracs_df = ibtracs_df.copy()
	emdat_df = emdat_df.copy()

	# Get names that occur in both datasets	
	ibtracs_names = set(ibtracs_df["Name"])
	emdat_names = set(emdat_df["Event Name"])
	matching_names = ibtracs_names.intersection(emdat_names)

	ibtracs_match_df = ibtracs_df[ibtracs_df["Name"].isin(matching_names)].copy()
	emdat_match_df = emdat_df[emdat_df["Event Name"].isin(matching_names)].copy()

	# Get storms that match on both name and year
	ibtracs_match_df["Year, Name"] = (ibtracs_match_df["Year"].astype(str) 
									  + ", " + ibtracs_match_df["Name"])
	emdat_match_df["Year, Name"] = (emdat_match_df["Year"].astype(str) 
									+ ", " + emdat_match_df["Event Name"])
	matches_ibtracs = ibtracs_match_df["Year, Name"].isin(emdat_match_df["Year, Name"])
	matches_emdat = emdat_match_df["Year, Name"].isin(ibtracs_match_df["Year, Name"])
	ibtracs_match_df = ibtracs_match_df[matches_ibtracs]
	emdat_match_df = emdat_match_df[matches_emdat]
	
	# Return remaining storms for further processing
	ibtracs_no_match_df = ibtracs_df.drop(ibtracs_match_df.index)
	emdat_no_match_df = emdat_df.drop(emdat_match_df.index)

	return ibtracs_match_df, emdat_match_df, ibtracs_no_match_df, emdat_no_match_df


# ============================================================================== 
crs = 4326						# CRS of storm data
# ============================================================================== 

# ============================================================================== 
# Get paths for loading and saving data
# ============================================================================== 
current_dir = os.getcwd()
top_dir = os.path.abspath(os.path.join(current_dir, os.pardir))
ibtracs_dir = top_dir + "/data/storm_data/IBTrACS"
preprocessed_dir = ibtracs_dir + "/IBTrACS_preprocessed"
damage_dir = top_dir + "/data/EM_DAT"
damage_file = damage_dir + "/EM_DAT_1980_2021.xlsx" 

# Load preprocessed IBTrACS data (resulting from ibtracs_preprocessing.py)
latlist = np.load(preprocessed_dir + "/LATLIST_INTERP.npy", allow_pickle=True).item()
lonlist = np.load(preprocessed_dir + "/LONLIST_INTERP.npy", allow_pickle=True).item()
timelist = np.load(preprocessed_dir + "/TIMELIST_INTERP.npy", allow_pickle=True).item()
windlist = np.load(preprocessed_dir + "/WINDLIST_INTERP.npy", allow_pickle=True).item()
preslist = np.load(preprocessed_dir + "/PRESLIST_INTERP.npy", allow_pickle=True).item()
rmaxlist = np.load(preprocessed_dir + "/RMAXLIST_INTERP.npy", allow_pickle=True).item()
basinlist = np.load(preprocessed_dir + "/BASINLIST_INTERP.npy", allow_pickle=True).item()
yearlist = np.load(preprocessed_dir + "/YEARLIST_INTERP.npy", allow_pickle=True).item()
namelist = np.load(preprocessed_dir + "/NAMELIST_INTERP.npy", allow_pickle=True).item()


# Create IBTrACS dataframe from variable lists
ibtracs_df = pd.DataFrame({"Nr": [key for key in namelist.keys()],
						   "Name": [name[0].lower() for name in namelist.values()], 
						   "Year": [int(year[0]) for year in yearlist.values()],
						   "Time (start)": [time[0] for time in timelist.values()],
						   "Time (end)": [time[-1] for time in timelist.values()],
						   "Basin": [basin[0] for basin in basinlist.values()]})

# Load EM-DAT damage data
emdat_df = pd.read_excel(damage_file, header=6, 
						 usecols=["Dis No", "Event Name", "Year", 
						 		  "Total Damages ('000 US$)", "Country", "ISO",
						 		  "Start Year", "Start Month", "Start Day",
						 		  "End Year", "End Month", "End Day"])
emdat_df["Start Month"] = emdat_df["Start Month"].fillna(0).astype(int)
emdat_df["Start Day"] = emdat_df["Start Day"].fillna(0).astype(int)
emdat_df["End Month"] = emdat_df["End Month"].fillna(0).astype(int)
emdat_df["End Day"] = emdat_df["End Day"].fillna(0).astype(int)
emdat_df["Event Name"] = emdat_df["Event Name"].fillna("")

# Use only data where damage data is available
emdat_df = emdat_df[emdat_df["Total Damages ('000 US$)"].notna()]

# Preprocess storm names in EM-DAT data
emdat_df = preprocess_emdat_names(emdat_df)

# Return remaining damage datapoints
ibtracs_match_df, emdat_match_df, _, _ = match_storm_names(emdat_df, ibtracs_df)

# Convert start and end times of storms to datetime
# NOTE: IBTrACS dates are given as days from 1858-11-17 00:00:00
start_date = datetime.datetime(1858, 11, 17)
ibtracs_df["Time (start)"] = start_date + pd.to_timedelta(ibtracs_df["Time (start)"], 
											 			  unit="days")
ibtracs_df["Time (end)"] = start_date + pd.to_timedelta(ibtracs_df["Time (end)"],
														unit="days")
ibtracs_match_df["Time (start)"] = (start_date + 
									pd.to_timedelta(ibtracs_match_df["Time (start)"], 
											 		unit="days"))
ibtracs_match_df["Time (end)"] = (start_date + 
								  pd.to_timedelta(ibtracs_match_df["Time (end)"],
												  unit="days"))

emdat_df["Time (start)"] = (emdat_df["Start Year"].astype(str) + "-"
						    + emdat_df["Start Month"].astype(str) + "-"
						    + emdat_df["Start Day"].astype(str))
emdat_df["Time (end)"] = (emdat_df["End Year"].astype(str) + "-"
						  + emdat_df["End Month"].astype(str) + "-"
						  + emdat_df["End Day"].astype(str))
emdat_df = emdat_df.drop(["Start Year", "Start Month", "Start Day", "End Year",
						  "End Month", "End Day"], axis=1)

emdat_match_df["Time (start)"] = (emdat_match_df["Start Year"].astype(str) + "-"
						   		  + emdat_match_df["Start Month"].astype(str) + "-"
						   		  + emdat_match_df["Start Day"].astype(str))
emdat_match_df["Time (end)"] = (emdat_match_df["End Year"].astype(str) + "-"
						   		+ emdat_match_df["End Month"].astype(str) + "-"
						   		+ emdat_match_df["End Day"].astype(str))
emdat_match_df = emdat_match_df.drop(["Start Year", "Start Month", "Start Day", 
								"End Year", "End Month", "End Day"], axis=1)


# Merge automatic and manual name matching dataframes
total_df1 = pd.merge(emdat_match_df, ibtracs_match_df[["Year, Name", "Nr", "Basin"]])

total_df1 = total_df1[["Nr", "Year", "Event Name", "ISO", "Country", "Basin", 
					   "Total Damages ('000 US$)"]]

match_dir = top_dir + "/data/storm_data/Storms_with_damage_data"
total_df2 = pd.read_csv(match_dir + "/EMDAT_matches.csv")
total_df2 = pd.merge(total_df2, ibtracs_df[["Nr", "Basin"]])
total_df2 = total_df2[["Nr", "Year", "Event Name", "ISO", "Country", "Basin", 
					   "Total Damages ('000 US$)"]]
total_df = total_df1.append(total_df2)
total_df["Damage"] = total_df["Total Damages ('000 US$)"] * 1000
total_df = total_df.drop("Total Damages ('000 US$)", axis=1)
total_df = total_df.set_index("Nr")
total_df.to_csv(match_dir + "/storm_damage_matches_ALL.csv")


# Get matching storms and save them to new directory
match_dir = (top_dir + "/data/storm_data/Storms_with_damage_data/EPSG_" + 
			 str(crs) + "/Original_Res")
old_storm_dir = ibtracs_dir + "/Historical_Storms"

if not os.path.exists(match_dir):
    print("Creating directory " + match_dir)
    os.makedirs(match_dir)

storm_indices = total_df.index
for storm in os.listdir(old_storm_dir):
	if int(storm.split("_")[-2][2:]) in storm_indices:
		shutil.copy(old_storm_dir + "/" + storm, match_dir)

