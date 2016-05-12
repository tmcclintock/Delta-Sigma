"""
This is a script that can be used to visualize outputs.
"""

import numpy as np
import matplotlib.pyplot as plt

testpath = "../test_data/matter_power_%s/"
ps = "nl"
k = np.genfromtxt(testpath%ps+"k_h.txt")
P = np.genfromtxt(testpath%ps+"p_k.txt")
plt.loglog(k,P)
plt.show()

R = np.genfromtxt("R.txt")
xi_mm = np.genfromtxt("xi_mm.txt")

plt.loglog(R,xi_mm)
plt.show()
