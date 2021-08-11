#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Created on Mon June 28 2021

@author: Liz Verbeek
"""

import os
import fiona
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# ============================================================================== 
# Get paths 
# ============================================================================== 
current_dir = os.getcwd()
top_dir = os.path.abspath(os.path.join(current_dir, os.pardir))
results_dir = top_dir + "/results"

# Specify target results directory
target_dir = "/Current_climate_level1regions_sigmoidfunc"
results_files_dir = results_dir + target_dir

AAD_dict = {}
for file in os.listdir(results_files_dir):
	df = pd.read_csv(results_files_dir + "/" + file)
	df = df.rename(columns={"Unnamed: 0": "ISO"})
	return_periods = np.array([int(column[1].lstrip("RP")) 
					  		   for column in df.columns.str.split("_") 
					  		   if column != ["ISO"]])
	annual_damage = df.values[:,1:]/return_periods.T
	df["AAD"] = np.sum(annual_damage, axis=1)
	df = df.set_index("ISO")
	AAD_dict = {**AAD_dict, **df["AAD"].to_dict()}

AAD_df = pd.DataFrame(AAD_dict.items(), columns=["ISO", "Average Annual Damage"])
AAD_df = AAD_df.set_index(["ISO"])
AAD_df.to_csv(results_dir + "/risk_assessment_" + target_dir[1:].lower() + ".csv")