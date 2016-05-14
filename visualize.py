"""
This is a script that can be used to visualize outputs.
"""

import numpy as np
import matplotlib.pyplot as plt

dopowerspec = False
doxinfw = False
doximm = False
dobias = False
dosigmar = True

testpath = "test_data/matter_power_%s/"
outpath = "output/"
ps = "nl"

k = np.genfromtxt(testpath%ps+"k_h.txt")
P = np.genfromtxt(testpath%ps+"p_k.txt")
R = np.genfromtxt(outpath+"R.txt")
xi_mm = np.genfromtxt(outpath+"xi_mm.txt")
xi_nfw = np.genfromtxt(outpath+"xi_nfw.txt")
xi_2h = np.genfromtxt(outpath+"xi_2halo.txt")
xi_hm = np.genfromtxt(outpath+"xi_hm.txt")
sigma_r = np.genfromtxt(outpath+"sigma_r.txt")
delta_sigma = np.genfromtxt(outpath+"delta_sigma.txt")
mis_sigma_r = np.genfromtxt(outpath+"miscentered_sigma_r.txt")

M = np.genfromtxt(outpath+"M.txt")
bias = np.genfromtxt(outpath+"bias.txt")
nu = np.genfromtxt(outpath+"nu.txt")

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
    plt.loglog(R,xi_mm,ls='--',alpha=0.5)
    plt.loglog(R,xi_hm)
    plt.loglog(R,xi_2h,ls=':',alpha=0.5)
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
    plt.loglog(R,delta_sigma)
    plt.loglog(R,mis_sigma_r,ls='--')
    plt.show()
    plt.clf()
