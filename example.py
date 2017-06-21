"""
This is an example file of how to import
and run the delta sigma code.

1) Set the path to the .so file, which is located in ./src/wrapper
"""
import sys
sys.path.insert(0,"./src/wrapper/")
import py_Delta_Sigma
import matplotlib.pyplot as plt
import numpy as np
plt.rc("font", size=12)

"""
2) Load a linear and mm power spectrum from somewhere.
In this example it's pre-computed, but one could use CAMB or CLASS

3) Create a dictionary with the cosmology.
For reference, fox has h=0.670435, om:0.31834
"""
klin = np.loadtxt("./test_data/matter_power_lin/k_h.txt")
Plin = np.loadtxt("./test_data/matter_power_lin/p_k.txt")
knl  = np.loadtxt("./test_data/matter_power_nl/k_h.txt")
Pnl  = np.loadtxt("./test_data/matter_power_nl/p_k.txt")
cosmo = {"h":0.7,"om":0.3,"ok":0.0}
cosmo["ode"]=1.0-cosmo["om"]

"""
4) Create a dictionary with parameters specifying
the profile.
Mass - mass of the halo Msun/h
concentration - concentration of the halo
delta - overdensity
NR - number of radial points you want to sample
Rmin - minimum R
Rmax - maximum R
Nbins - number of (logarithmic) radial bins to do averaging over
R_bin_min - left edge of the left most bin
R_bin_max - right edge of the right most bin
Rmis - width of the miscentering 2D Gaussian
fmis - fraction of miscentered halos
miscentering - 1 if you want to run miscentering, 0 if not
averaging - 1 if you want to run averaging, 0 if not
single_miscentering - 1 if you want miscentering of a single cluster, 0 if not
"""
input_params = {"Mass": 3*10**14,"NR":300,"Rmin":0.01,
                "Rmax":200.0,"Nbins":15,"R_bin_min":0.01,"R_bin_max":200.0,
                "delta":200,"Rmis":0.25,"fmis":0.25,
                "miscentering":1,"averaging":1,"single_miscentering":1}
input_params["concentration"] = 5.0 #Arbitrary

"""
5) Call the calc_Delta_Sigma function, which returns a dictionary.
"""
return_dict = py_Delta_Sigma.calc_Delta_Sigma(klin,Plin,knl,Pnl,cosmo,input_params)
print "The output dictionary contains:"
for key in return_dict.keys(): print "\t",key

"""
6) Pull out quantities that we want to plot.
Nomenclature:
- sigma and delta_sigma are centered
- ave_ means the quantity is averaged over radial bins
- miscentered_ means the quantity has been convolved with a 2D gaussian of width Rmis
- _mis means the quantity is miscentered by an amount Rmis
- full_ means the quantity is a mean weighted by fmis of the centered and miscentered part
"""
R = return_dict["R"]

"""
sigma = return_dict['sigma']
sigma_miscentered = return_dict['miscentered_sigma']
sigma_m = return_dict['sigma_mis']
plt.loglog(R,sigma,label=r"$\Sigma$")
plt.loglog(R,sigma_miscentered,label=r"$\Sigma_{\rm mis}$",ls='--')
plt.loglog(R,sigma_m,label=r"$\Sigma(R|R_{\rm mis})$")
plt.legend()
plt.xlabel(r"$R\ [{\rm Mpc}/h]$",fontsize=24)
plt.ylabel(r"$\Sigma\ [{\rm M_\odot}\ h/{\rm pc^2}]$",fontsize=24)
plt.subplots_adjust(bottom=0.15, left=0.15)
plt.show()
plt.clf()
"""

delta_sigma = return_dict['delta_sigma']
miscentered_delta_sigma = return_dict['miscentered_delta_sigma']
delta_sigma_mis = return_dict['delta_sigma_mis']
plt.loglog(R,delta_sigma,label=r"$\Delta\Sigma_{\rm centered}$")
plt.loglog(R,miscentered_delta_sigma,label=r"$\Delta\Sigma_{\rm mis}$",ls='--')
plt.loglog(R,delta_sigma_mis,label=r"$\Delta\Sigma(R|R_{\rm mis})$")

plt.legend(fontsize=12)
plt.xlabel(r"$R\ [{\rm Mpc}/h]$",fontsize=24)
plt.ylabel(r"$\Delta\Sigma\ [{\rm M_\odot}\ h/{\rm pc^2}]$",fontsize=24)
plt.subplots_adjust(bottom=0.15)
plt.show()
plt.clf()

Rbins = return_dict['Rbins']
ave_delta_sigma = return_dict['ave_delta_sigma']
ave_miscentered_delta_sigma = return_dict['ave_miscentered_delta_sigma']
ave_delta_sigma_mis = return_dict['ave_delta_sigma_mis']
plt.loglog(Rbins,ave_delta_sigma,label=r"$\bar{\Delta\Sigma}_{\rm centered}$")
plt.loglog(Rbins,ave_miscentered_delta_sigma,label=r"$\bar{\Delta\Sigma}_{\rm mis}$",ls='--')
plt.loglog(Rbins,ave_delta_sigma_mis,label=r"$\bar{\Delta\Sigma}(R|R_{\rm mis})$")
plt.legend()
plt.xlabel(r"$R\ [{\rm Mpc}/h]$",fontsize=24)
plt.ylabel(r"$\Delta\Sigma\ [{\rm M_\odot}\ h/{\rm pc^2}]$",fontsize=24)
plt.subplots_adjust(bottom=0.15)
plt.show()
plt.clf()
