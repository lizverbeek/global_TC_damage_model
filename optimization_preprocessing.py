#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Created on Sat May 15 2021

@author: Liz Verbeek
"""

import os
import fiona
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from matplotlib import cm


def append_sim_damage(rep_damage_file, sim_damage_dir, param_name, year):
	""" Append simulated damage from files to dataframe with reported damage.

	Arguments:
		rep_damage_file {string}		-- Path to file with reported damage
		sim_damage_dir {string}			-- Directory with simulated damage files
		param_name {string}				-- Name of parameter which is optimized

	Returns:
		new_df {dataframe}				-- New dataframe with reported and
										   simulated damage

	"""
	rep_damage_df = pd.read_csv(rep_damage_file, 
								usecols=["Nr", "ISO", "Damage year " + str(year)])
	for file in os.listdir(sim_damage_dir):
		if param_name == "vthresh":
			param_str = file.split(".")[0].split("_")[-2].strip(param_name)
			param_value = float(param_str[:-1] + "." + param_str[-1])
		elif param_name == "vhalf":
			param_str = file.split(".")[0].split("_")[-1].strip(param_name)
			param_value = float(param_str[:-1] + "." + param_str[-1])
		# elif param_name == "grid":
		# 	vthresh_str = file.split(".")[0].split("_")[-2].strip("vthresh")
		# 	vhalf_str = file.split(".")[0].split("_")[-1].strip("vhalf")
		# 	vthresh_value = float(vthresh_str[:-1] + "." + vthresh_str[-1])
		# 	vhalf_value = float(vhalf_str[:-1] + "." + vhalf_str[-1])
		elif param_name == "exp":
			factor_str = file.split(".")[0].split("_")[-2].strip("factor")
			exp_str = file.split(".")[0].split("_")[-1].strip("exp")
			factor_value = float(factor_str[:1] + "." + factor_str[1:])
			exp_value = float(exp_str[:-1] + "." + exp_str[-1])

		sim_damage_df = pd.read_csv(sim_damage_dir + "/" + file)
		sim_damage_all = []
		for storm, country in zip(rep_damage_df["Nr"].values, 
								  rep_damage_df["ISO"].values):
			sim_damage = (sim_damage_df[sim_damage_df["ISO"] == 
										country][str(storm)].values[0])
			sim_damage_all.append(sim_damage)

		new_df = rep_damage_df
		if param_name == "exp":
			new_df["exp " + str(exp_value) + 
				   " factor " + str(factor_value)] = sim_damage_all
		# elif param_name == "grid":
		# 	new_df["vthresh " + str(vthresh_value) + 
		# 		   " vhalf " + str(vhalf_value)] = sim_damage_all	
		else:
			new_df[param_name + " " + str(param_value)] = sim_damage_all
	
	new_df = new_df.set_index("Nr")
	new_df.to_csv(os.path.abspath(os.path.join(sim_damage_dir, os.pardir)) 
				 + "/total_damage_" + param_name + "_range.csv")

	return new_df


# ============================================================================== 
# Get paths 
# ============================================================================== 
current_dir = os.getcwd()
top_dir = os.path.abspath(os.path.join(current_dir, os.pardir))

# Load reported and simulated damage
rep_damage_file = (top_dir + "/data/storm_data/Storms_with_damage_data" 
				   + "/storm_damage_matches_ALL.csv")
results_dir = top_dir + "/results"
max_damage_dir = top_dir + "/data/max_damage_data"
income_data_file = max_damage_dir + "/METADATA_GDP.csv"

vthresh_dir = results_dir + "/v_thresh_range"
vhalf_dir = results_dir + "/v_half_range"
# grid_dir = results_dir + "/grid_range"
exp_dir = results_dir + "/exp_range"

# Convert to single dataframe per experiment
# damage_df = append_sim_damage(rep_damage_file, vthresh_dir, "vthresh", 2000)
damage_df = append_sim_damage(rep_damage_file, vhalf_dir, "vhalf", 2000)
# damage_df = append_sim_damage(rep_damage_file, grid_dir, "grid")
# damage_df = append_sim_damage(rep_damage_file, exp_dir, "exp", 2000)