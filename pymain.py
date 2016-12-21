"""
This is an example file of how to interface with the python wrapper
that can run the delta sigma code.

STEP 0: compile with:
make SHARED=yes

This should result in a main.so file.
"""

"""
First, we have to give the path to the py_Delta_Sigma.py file,
which is located in src/wrapper/

One day I'll figure out how to write a setup.py script and
have it install in the python directory
"""
import sys
sys.path.insert(0,"src/wrapper/")
import py_Delta_Sigma
import matplotlib.pyplot as plt
import numpy as np

test_path = "test_data/"

#First load in the power spectrum data from somewhere
klin = np.genfromtxt(test_path+"matter_power_lin/k_h.txt")
Plin = np.genfromtxt(test_path+"matter_power_lin/p_k.txt")
knl  = np.genfromtxt(test_path+"matter_power_nl/k_h.txt")
Pnl  = np.genfromtxt(test_path+"matter_power_nl/p_k.txt")

#Create a dictionary with the cosmology
cosmo = {"h":0.7,"om":0.3,"ok":0.0}
cosmo["ode"]=1.0-cosmo["om"]

#This is the buzzard cosmology.
#cosmo = {"h":0.7,"om":0.286,"ok":0.0}
#cosmo["ode"]=1.0-cosmo["om"]

#This is the fox sim cosmology
#cosmo = {"h":0.670435,"om":0.31834,"ok":0.0}
#cosmo["ode"]=1.0-cosmo["om"]

"""
Create a dictionary with all starting params. They are:
Mass - mass of the halo
concentration - concentration of the halo
NR - number of radial points you want to sample
Rmin - minimum R
Rmax - maximum R
Nbins - number of (logarithmic) radial bins to do averaging over
R_bin_min - left edge of the left most bin
R_bin_max - right edge of the right most bin
delta - overdensity
Rmis - width of the miscentering 2D Gaussian
fmis - fraction of miscentered halos
timing - 1 if you want to print timing information, 0 if not
miscentering - 1 if you want to run miscentering, 0 if not
averaging - 1 if you want to run averaging, 0 if not
"""
input_params = {"Mass": 2.56*10**14,"NR":300,"Rmin":0.01,
                "Rmax":200.0,"Nbins":15,"R_bin_min":0.01,"R_bin_max":200.0,
                "delta":200,"Rmis":0.249697,"fmis":0.23374,
                "timing":1,"miscentering":1,"averaging":1}
input_params["concentration"] = 5.0 #Completely arbitrary
#input_params["concentration"] = 4.0*(input_params["Mass"]/5.e14)**-0.1
#Above is an example M-c relation.

#Results come out in a dictionary
return_dict = py_Delta_Sigma.calc_Delta_Sigma(klin,Plin,knl,Pnl,cosmo,input_params)
print "The return dictionary contains:"
for key in return_dict.keys():
    print "\t",key
print ""

R = return_dict["R"]

sigma = return_dict['sigma']
sigma_miscentered = return_dict['miscentered_sigma']
sigma_m = return_dict['sigma_mis']
plt.loglog(R,sigma,label=r"$\Sigma$")
plt.loglog(R,sigma_miscentered,label=r"$\Sigma_{\rm mis}$",ls='--')
plt.loglog(R,sigma_m,label=r"$\Sigma(R|R_{\rm mis})$")
plt.legend()
plt.xlabel(r"$R\ [{\rm Mpc}/h]$",fontsize=24)
plt.ylabel(r"$\Sigma\ [{\rm M_\odot}\ h/{\rm pc^2}]$",fontsize=24)
plt.subplots_adjust(bottom=0.15)
#plt.show()
plt.clf()

delta_sigma = return_dict['delta_sigma']
miscentered_delta_sigma = return_dict['miscentered_delta_sigma']
delta_sigma_mis = return_dict['delta_sigma_mis']
plt.loglog(R,delta_sigma,label=r"$\Delta\Sigma_{\rm centered}$")
plt.loglog(R,miscentered_delta_sigma,label=r"$\Delta\Sigma_{\rm mis}$",ls='--')
plt.loglog(R,delta_sigma_mis,label=r"$\Delta\Sigma(R|R_{\rm mis})$")
plt.legend()
plt.xlabel(r"$R\ [{\rm Mpc}/h]$",fontsize=24)
plt.ylabel(r"$\Delta\Sigma\ [{\rm M_\odot}\ h/{\rm pc^2}]$",fontsize=24)
plt.subplots_adjust(bottom=0.15)
#plt.show()
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
