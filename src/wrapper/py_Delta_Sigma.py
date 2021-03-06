"""
This is the python wrapper that calls the python_wrapper()
c function. This interfaces through c_types so that the user
doesn't have to.
"""
import numpy as np
from ctypes import c_double,c_int,POINTER,cdll
import inspect
import os
library_path = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))+"/Delta_Sigma.so"
dslib = cdll.LoadLibrary(library_path)
interface = dslib.python_interface
interface.restype = c_int
"""
Arguments to the interface are: 
k_lin, P_lin, N_k_lin,
k_nl, P_nl, N_k_nl,
NR,Rmin,Rmax,
h,om,
Mass,concentration,
delta,
averaging, Nbins,
R_bin_min,R_bin_max,

R,xi_1halo,xi_mm,xi_lin
xi_2halo,xi_hm,sigma,

delta_sigma,Rbins,
ave_delta_sigma,
bias,nu
"""
interface.argtypes=[POINTER(c_double),POINTER(c_double),c_int,
                    POINTER(c_double),POINTER(c_double),c_int,
                    c_int,c_double,c_double,
                    c_double,c_double,
                    c_double,c_double,
                    c_int,
                    c_int,c_int,
                    c_double,c_double,
                    POINTER(c_double),POINTER(c_double),POINTER(c_double),
                    POINTER(c_double),POINTER(c_double),POINTER(c_double),
                    POINTER(c_double),POINTER(c_double),
                    POINTER(c_double),
                    POINTER(c_double),POINTER(c_double)]

def calc_Delta_Sigma(k_lin,P_lin,k_nl,P_nl,cosmo_dict,params):
    """Calculates the DeltaSigma profile given some cosmology, matter power spectra, and input parameters (e.g. mass, concentraton, etc.)

    Note: Mass units are Msun/h. Distances are Mpc/h comoving.

    k_lin (array_like): Wavenumbers of input linear matter power spectrum; h/Mpc.
    P_lin (array_like): Linear matter power spectrum; (h/Mpc)^3.
    k_nl (array_like): Wavenumbers of input nonlinear matter power spectrum; h/Mpc.
    P_nl (array_like): Nonlinear matter power spectrum; (h/Mpc)^3.
    cosmo_dict (dictionary): Contains key-value pairs of cosmological parameters. Required parameters: h and om.
    params (dictionary): Contains key-value pairs of halo parameters, including: Mass, delta, concentration, NR, Rmin, Rmax, Nbins, R_bin_min, R_bin_max, averaging.

    Returns:
        output (dictionary): Contains key-value pairs of all possible quantities assosciated with the halo.
    """

    Nk_lin = len(k_lin)
    Nk_nl  = len(k_nl)
    k_lin_in = k_lin.ctypes.data_as(POINTER(c_double))
    P_lin_in = P_lin.ctypes.data_as(POINTER(c_double))
    k_nl_in = k_nl.ctypes.data_as(POINTER(c_double))
    P_nl_in = P_nl.ctypes.data_as(POINTER(c_double))

    Mass,concentration,delta = params["Mass"],params["concentration"],params['delta']
    NR,Rmin,Rmax = params["NR"],params["Rmin"],params["Rmax"]

    #Default values to pass to C
    Nbins = 2
    R_bin_min = Rmin
    R_bin_max = Rmax
    if "averaging" in params: 
        averaging = params['averaging']
        Nbins = params['Nbins']
        R_bin_min = params['R_bin_min']
        R_bin_max = params['R_bin_max']

    h,om = cosmo_dict['h'],cosmo_dict['om']

    R = np.zeros(NR)
    R_in = R.ctypes.data_as(POINTER(c_double))
    xi_1halo = np.zeros(NR)
    xi_1halo_in = xi_1halo.ctypes.data_as(POINTER(c_double))
    xi_mm = np.zeros(NR)
    xi_mm_in = xi_mm.ctypes.data_as(POINTER(c_double))
    xi_lin = np.zeros(NR)
    xi_lin_in = xi_lin.ctypes.data_as(POINTER(c_double))
    xi_2halo = np.zeros(NR)
    xi_2halo_in = xi_2halo.ctypes.data_as(POINTER(c_double))
    xi_hm = np.zeros(NR)
    xi_hm_in = xi_hm.ctypes.data_as(POINTER(c_double))
    sigma = np.zeros(NR)
    sigma_in = sigma.ctypes.data_as(POINTER(c_double))
    delta_sigma = np.zeros(NR)
    delta_sigma_in = delta_sigma.ctypes.data_as(POINTER(c_double))

    Rbins = np.zeros(Nbins)
    Rbins_in = Rbins.ctypes.data_as(POINTER(c_double))
    ave_delta_sigma = np.zeros(Nbins)
    ave_delta_sigma_in = ave_delta_sigma.ctypes.data_as(POINTER(c_double))

    bias = np.zeros(1)
    bias_in = bias.ctypes.data_as(POINTER(c_double))
    nu = np.zeros(1)
    nu_in = nu.ctypes.data_as(POINTER(c_double))
    
    result = interface(k_lin_in,P_lin_in,Nk_lin,
                       k_nl_in,P_nl_in,Nk_nl,
                       NR,Rmin,Rmax,
                       h,om,
                       Mass,concentration,
                       delta,
                       averaging,
                       Nbins,R_bin_min,R_bin_max,
                       R_in,xi_1halo_in,xi_mm_in,xi_lin_in,
                       xi_2halo_in,xi_hm_in,
                       sigma_in,delta_sigma_in,Rbins_in,
                       ave_delta_sigma_in,bias_in,nu_in)

    #Now build a dictionary and return it
    output = {"R":R,"xi_1halo":xi_1halo,"xi_mm":xi_mm,"xi_lin":xi_lin,
              "xi_2halo":xi_2halo,"xi_hm":xi_hm,
              "sigma":sigma,"delta_sigma":delta_sigma,"bias":bias,"nu":nu}
    
    if averaging:
        output["Rbins"] = Rbins
        output["ave_delta_sigma"] = ave_delta_sigma

    return output
