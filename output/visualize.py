"""
This is a script that can be used to visualize outputs.
"""

import numpy as np
import matplotlib.pyplot as plt

dopowerspec = False
doxinfw = True
doximm = False
dobias = False
dosigmar = False

testpath = "../test_data/matter_power_%s/"
ps = "nl"

k = np.genfromtxt(testpath%ps+"k_h.txt")
P = np.genfromtxt(testpath%ps+"p_k.txt")
R = np.genfromtxt("R.txt")
xi_mm = np.genfromtxt("xi_mm.txt").T[0]
xi_nfw = np.genfromtxt("xi_nfw.txt")
xi_hm = np.genfromtxt("xi_hm.txt")
sigma_r = np.genfromtxt("sigma_r.txt")
delta_sigma = np.genfromtxt("delta_sigma.txt")

M = np.genfromtxt("M.txt")
bias = np.genfromtxt("bias.txt")
nu = np.genfromtxt("nu.txt")

if dopowerspec:
    plt.loglog(k,P)
    plt.show()
    plt.clf()

if doximm:
    plt.loglog(R,xi_mm)
    plt.show()
    plt.clf()

if doxinfw:
    plt.loglog(R,xi_nfw)
    plt.loglog(R,xi_mm)
    plt.loglog(R,xi_hm)
    plt.loglog(R,sigma_r)
    plt.loglog(R,delta_sigma)
    plt.show()
    plt.clf()

if dobias:
    plt.plot(M,bias)
    plt.xscale('log')
    plt.show()
    plt.clf()
    plt.plot(M,nu)
    plt.xscale('log')
    plt.show()
    plt.clf()

if dosigmar:
    plt.loglog(R,sigma_r)
    plt.show()
    plt.clf()
