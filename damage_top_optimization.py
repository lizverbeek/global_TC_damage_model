#!/usr/bin/env python
# -*- coding: utf-8 -*-


import os
import numpy as np


# total_vthresh = np.arange(25.5, 30.0, 0.5)
# total_vhalf = np.arange(42.0, 43.0, 1.0)
total_vthresh = [0.0]
# total_vhalf = [65.7]
for vthresh in total_vthresh:
	for vhalf in total_vhalf:
		print("Running damage computation for parameters "
		      + str(vthresh) + "," + str(vhalf))
		os.system("sbatch --export=v_thresh=" + str(vthresh) + ",v_half="
			  + str(vhalf) + " damage_main.sh")

total_factor = np.arange(0.014, 0.015, 0.0005)
total_exp = np.arange(2, 9.5, 0.5)
for factor in total_factor:
	for exp in total_exp:
		os.system("sbatch --export=v_thresh=" + str(factor) + ",v_half=" 
				  + str(exp) + " damage_main.sh")