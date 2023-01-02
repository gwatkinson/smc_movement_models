# To compute variances sigma_epsilon and sigma_o

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from statsmodels.tsa.stattools import acovf

df = pd.read_csv("smc_movement_models/data/clean_data.csv")

# To build a time window sequenced vector
# day_str = "18/01/2008"
day_str_min = "18/01/2008"
day_str_max = "18/01/2008"
mask_day = (df["Date"] >= day_str_min) & (df["Date"] <= day_str_max)
day = df.loc[mask_day, ["Dtime", "Velocity", "Depth"]]
day["Dtime"] = pd.to_datetime(day["Dtime"])

window_size = 312
overlap_size = int(window_size / 2)
num_windows = int((len(day) - overlap_size) / overlap_size)

windows = []
ids = []
for i in range(num_windows):
    window = day[i * overlap_size : (i + 2) * overlap_size - 1]["Velocity"].values
    ids.append((i + 1) * overlap_size)
    windows.append(window)

# We need to compute the sigma value for each time window
# Polynomial regression close to the origin
sigma_zs = []
sigma_os = []
sigma_epsilons = []
acvf_list = []

for window in windows:
    acvf = acovf(window)
    acvf_list.append(acvf)

    x = np.linspace(0, len(acvf[1:9]), len(acvf[1:9]))
    coeffs = np.polyfit(x, acvf[1:9], 2)

    sigma_z = coeffs[2]
    sigma_o = acvf[0] - sigma_z
    sigma_epsilon = (1-(acvf[1]/acvf[0])**2)*sigma_z

    sigma_zs.append(sigma_z)
    sigma_os.append(sigma_o)
    sigma_epsilons.append(sigma_epsilon)