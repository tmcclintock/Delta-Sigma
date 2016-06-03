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

Mass = 1.e14
concentration = 4.0*(Mass/5e14)**-0.1
NR = 300
Rmin,Rmax = 0.01,200.0
Nbins = 15
R_bin_min, R_bin_max = 0.01,200.0
delta = 200
Rmis, fmis = 0.3,0.23
timing,miscentering,averaging = 1,1,1
