#!/usr/bin/env python
# -*- coding: utf-8 -*-


import os
import numpy as np


# Match countries to specified regions
basin_dict = {"AUS":"OC", "BGD":"NI", "CAN":"NA2", "CHN":"WP3", "CRI":"NA1", 
              "CUB":"NA1", "DOM":"NA1", "FJI":"OC", "GTM":"NA1", "HND":"NA1", 
              "HTI":"NA1", "IDN":"WP1", "IND":"NI", "JPN":"WP4", "KOR":"WP4", 
              "LAO":"WP1", "LKA":"NI", "MAC":"WP4", "MDG":"SI", "MEX":"NA1", 
              "MMR":"NI", "MOZ":"SI", "NIC":"NA1", "OMN":"NI", "PAK":"NI", 
              "PHL":"WP2", "PRK":"WP4", "SLV":"NA1", "SOM":"SI", "THA":"WP1", 
              "USA":"NA2", "VNM":"WP1", "VUT":"OC", "YEM":"NI", "ZWE":"SI"}

# Region information all countries 
basin_dict_new = {"AFG": "NI", "BTN": "NI", "KHM": "WP1", "COL": "NA1", 
                  "DJI": "NI", "ERI": "NI", "ETH": "NI", "GUY": "NA1", 
                  "IRN": "NI", "MWI": "SI", "MYS": "WP1", "MNG": "WP3", 
                  "NPL": "NI", "NZL": "OC", "PAN": "NA1", "PNG": "OC", 
                  "SLB": "OC", "ZAF": "SI", "SUR": "NA1", "SWZ": "SI", 
                  "TZA": "SI", "TLS": "OC", "ARE": "NI", "VEN": "NA1"}

basin_dict_new_new = {"AGO": "SI", "BRN": "WP1", "BWA": "SI", "NAM": "SI", 
                      "RUS": "WP4", "VIR": "NA1", "ZMB": "SI"}

basin_dict_other = {"BFA": "NA3", "CIV": "NA3", "DZA": "NA3", "ESP": "NA4", 
                    "FRA": "NA4", "GBR": "NA4", "GHA": "NA3", "GIN": "NA3", 
                    "GMB": "NA3", "GNB": "NA3", "IMN": "NA4", "LBR": "NA3", 
                    "MAR": "NA3", "MLI": "NA3", "MRT": "NA3", "PRT": "NA4", 
                    "SEN": "NA3", "SLE": "NA3"}

# Merge region definitions for full damage computation
basin_dict = {**basin_dict, **basin_dict_new}
basin_dict_2 = {**basin_dict_new_new, **basin_dict_other}
basin_dict_all = {**basin_dict, **basin_dict_2}

# # Run all countries
# countries = list(basin_dict_all.keys())
countries = list(basin_dict_2.keys())

vul_func = "sigmoid"
for country in countries:
    os.system("sbatch --export=vul_func=" + vul_func + ",country=" + country 
              + " damage_main.sh")
