import scipy as sci
import pandas as pd
import numpy as np
import scipy.signal
import os

def peak_calling(dir_to_WPS):
	for root, dirs, files in os.walk(dir_to_WPS):
		for file in files:
			print(file)
			path_to_WPS = os.path.join(root, file) # abs path to df that contains WPS information
			df = pd.read_csv(path_to_WPS, sep = "\t")

			idx = sci.signal.find_peaks(df["WPS"])[0]

			peak = np.array([False]*len(df["WPS"])) # an array full of False
			peak[idx] = True

			df["peak"] = peak
			df.to_csv(path_to_WPS, sep = "\t", index = False) # save a new copy of the df with the peak info

peak_calling(dir_to_WPS = "/groups/wyattgrp/users/amunzur/mt-ctDNA/WGS/WPS/data/WPS_files/window_size_15")







