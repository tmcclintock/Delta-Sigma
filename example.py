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
klin = np.loadtxt("./test_data/k.txt")
Plin = np.loadtxt("./test_data/plin_z0_l3.txt")
Pnl = np.loadtxt("./test_data/pnl_z0_l3.txt")
knl = np.copy(klin)
#klin = np.loadtxt("./test_data/matter_power_lin/k_h.txt")
#Plin = np.loadtxt("./test_data/matter_power_lin/p_k.txt")
#knl  = np.loadtxt("./test_data/matter_power_nl/k_h.txt")
#Pnl  = np.loadtxt("./test_data/matter_power_nl/p_k.txt")
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
averaging - 1 if you want to run averaging, 0 if not
"""
input_params = {"Mass": 3*10**14,"NR":300,"Rmin":0.01,
                "Rmax":200.0,"Nbins":15,"R_bin_min":0.01,"R_bin_max":200.0,
                "delta":200,"averaging":1}
input_params["concentration"] = 5.0 #Arbitrary
input_params["concentration"] = 4.916498 #Arbitrary


"""
5) Call the calc_Delta_Sigma function, which returns a dictionary.
"""
params = py_Delta_Sigma.calc_Delta_Sigma(klin,Plin,knl,Pnl,cosmo,input_params)
print "The output dictionary contains:"
for key in params.keys(): print "\t",key

"""
6) Pull out quantities that we want to plot.
Nomenclature:
- sigma and delta_sigma are centered
- ave_ means the quantity is averaged over radial bins
"""
R = params["R"]
sigma = params['sigma']
delta_sigma = params['delta_sigma']
Rbins = params['Rbins']
ave_delta_sigma = params['ave_delta_sigma']
plt.loglog(R,sigma,label=r"$\Sigma$")
plt.loglog(R,delta_sigma,label=r"$\Delta\Sigma$")
plt.loglog(Rbins,ave_delta_sigma,label=r"$\bar{\Delta\Sigma}$")
plt.legend()
plt.xlabel(r"$R\ [{\rm Mpc}/h]$",fontsize=24)
plt.ylabel(r"$[{\rm M_\odot}\ h/{\rm pc^2}]$",fontsize=24)
plt.subplots_adjust(bottom=0.15, left=0.15)
plt.show()
plt.clf()



