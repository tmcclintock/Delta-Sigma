"""
This is an example file of how to interface with the python wrapper
that can run the delta sigma code.
"""

import sys
sys.path.insert(0,"src/wrapper/")
import py_Delta_Sigma

import numpy as np

#First load in the power spectrum data
klin = np.genfromtxt("test_data/matter_power_lin/k_h.txt")
Plin = np.genfromtxt("test_data/matter_power_lin/p_k.txt")
knl = np.genfromtxt("test_data/matter_power_nl/k_h.txt")
Pnl = np.genfromtxt("test_data/matter_power_nl/p_k.txt")

cosmo = {"h":0.7,"om":0.286,"ok":0.0,"ode":0.7}
cosmo["ode"]=1.0-cosmo["om"]

input_params = {"Mass": 1.e14,"NR":300,"Rmin":0.01,"Rmax":200.0,"Nbins":15,"R_bin_min":0.01,"R_bin_max":200.0,"delta":200,"Rmis":0.3,"fmis":0.23,"timing":1,"miscentering":1,"averaging":1}
input_params["concentration"] = 4.0*(input_params["Mass"]/5.e14)**-0.1


py_Delta_Sigma.calc_Delta_Sigma(klin,Plin,knl,Pnl,cosmo,input_params)
