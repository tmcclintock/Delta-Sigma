"""
This is a script that can be used to visualize outputs.
"""

import numpy as np
import matplotlib.pyplot as plt

dopowerspec = False
doxinfw = True
doximm = False
dobias = False

testpath = "../test_data/matter_power_%s/"
ps = "nl"
k = np.genfromtxt(testpath%ps+"k_h.txt")
P = np.genfromtxt(testpath%ps+"p_k.txt")
if dopowerspec:
    plt.loglog(k,P)
    plt.show()
    plt.clf()

R = np.genfromtxt("R.txt")
xi_mm = np.genfromtxt("xi_mm.txt").T[0]
if doximm:
    plt.loglog(R,xi_mm)
    plt.show()
    plt.clf()

xi_nfw = np.genfromtxt("xi_nfw.txt")
if doxinfw:
    plt.loglog(R,xi_nfw)
    plt.loglog(R,xi_mm)
    plt.show()
    plt.clf()


M = np.genfromtxt("M.txt")
bias = np.genfromtxt("bias.txt")
nu = np.genfromtxt("nu.txt")
if dobias:
    plt.plot(M,bias)
    plt.xscale('log')
    plt.show()
    plt.clf()
    plt.plot(M,nu)
    plt.xscale('log')
    plt.show()
    plt.clf()
