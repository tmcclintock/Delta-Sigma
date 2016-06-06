"""
This is the python wrapper that calls the python_wrapper()
c function. This interfaces through c_types so that the user
doesn't have to.
"""
import numpy as np
import ctypes
from ctypes import c_double,c_int,POINTER,cdll

def calc_Delta_Sigma(klin,Plin,knl,Pnl,cosmo_dict,input_params):
    dslib = cdll.LoadLibrary("main.so")
    interface = dslib.python_interface
    interface.restype = c_int
    """
    Arguments are: 
    k_lin, P_lin, N_k_lin,
    k_nl, P_nl, N_k_nl,
    NR,Rmin,Rmax,
    h,om,ode,ok,
    Mass,concentration,
    Rmis,fmis,delta,
    flow_control,timing,miscentering,
    averaging, Nbins,
    R_bin_min,R_bin_max,

    R,xi_1halo,xi_mm,
    xi_2halo,xi_hm,sigma_r,

    delta_sigma,Rbins,
    ave_delta_sigma,bias,

    nu,miscentered_sigma_r,
    miscentered_delta_sigma,
    miscentered_ave_delta_sigma
    """
    interface.argtypes=[POINTER(c_double),POINTER(c_double),c_int,\
                        POINTER(c_double),POINTER(c_double),c_int,\
                        c_int,c_double,c_double,\
                        c_double,c_double,c_double,c_double,\
                        c_double,c_double,\
                        c_double,c_double,c_int,\
                        POINTER(c_int),c_int,c_int,\
                        c_int,c_int,\
                        c_double,c_double,\
                        POINTER(c_double),POINTER(c_double),POINTER(c_double),\
                        POINTER(c_double),POINTER(c_double),POINTER(c_double),\
                        POINTER(c_double),POINTER(c_double),\
                        POINTER(c_double),POINTER(c_double),\
                        POINTER(c_double),POINTER(c_double),\
                        POINTER(c_double),\
                        POINTER(c_double)]
    klin_in = np.array(klin,dtype=np.double).ctypes.data_as(POINTER(c_double))
    Plin_in = np.array(Plin,dtype=np.double).ctypes.data_as(POINTER(c_double))
    knl_in = np.array(knl,dtype=np.double).ctypes.data_as(POINTER(c_double))
    Pnl_in = np.array(Pnl,dtype=np.double).ctypes.data_as(POINTER(c_double))
    print klin_in

    Mass,concentration,NR,Rmin,Rmax,Nbins,R_bin_min,R_bin_max,\
        delta,Rmis,fmis,timing,miscentering,\
        averaging = input_params["Mass"],input_params["concentration"],\
                    input_params["NR"],input_params["Rmin"],\
                    input_params["Rmax"],input_params["Nbins"],\
                    input_params["R_bin_min"],input_params["R_bin_max"],\
                    input_params["delta"],input_params["Rmis"],\
                    input_params["fmis"],input_params["timing"],\
                    input_params["miscentering"],input_params["averaging"]

    Rarr = np.zeros(NR)
    binarr = np.zeros(Nbins)
    
    

    print Mass,concentration,Rarr.shape,binarr.shape
