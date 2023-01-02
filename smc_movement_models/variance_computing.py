# To compute variances sigma_epsilon and sigma_o

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from statsmodels.tsa.stattools import acovf

df = pd.read_csv(".../data/clean_data.csv")

acvf = acovf(df.Velocity.values)

# Polynomial regression, general framework

x = np.linspace(0, len(acvf), len(acvf))
coeffs = np.polyfit(x, acvf, 2)

y = coeffs[0] + coeffs[1] * x + coeffs[2] * x**2

plt.plot(x,y)
plt.show()

# Polynomial regression close to the origin

x = np.linspace(0, len(acvf[0:9]), len(acvf[0:9]))
coeffs = np.polyfit(x, acvf[0:9], 2)

def regfit2(k):
    """
    k must be an integer (lag)
    The function is equal to the estimate of sigma_z for all k > 0
    """
    return coeffs[0] + coeffs[1] * k + coeffs[2] * k**2

sigma_z0 = regfit2(0)
sigma_o = acvf[0] - sigma_z0

sigma_epsilon = (1-(acvf[1]/acvf[0])**2)*sigma_z0