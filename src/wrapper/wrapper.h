/*
  This is an interface between the DeltaSigma code and 
  a piece of code that can pass in a cosmology and 
  a power spectrum.
*/
#include "../mc_relation/mc_relation.h"
#include "../delta_sigma/delta_sigma.h"
#include "../sigma/sigma.h"
#include "../xi_hm/xi_hm.h"
#include "../xi_2halo/xi_2halo.h"
#include "../xi_mm/xi_mm.h"
#include "../tinker_bias/tinker_bias.h"
#include "../sigma_mis/sigma_mis.h"
#include "../delta_sigma_mis/delta_sigma_mis.h"
#include "../miscentered_sigma/miscentered_sigma.h"
#include "../miscentered_delta_sigma/miscentered_delta_sigma.h"
#include "../ave_delta_sigma/ave_delta_sigma.h"
#include "../ave_miscentered_delta_sigma/ave_miscentered_delta_sigma.h"
#include "../cosmology/cosmology.h"

#ifndef INTERFACE
#define INTERFACE
typedef struct interface_parameters{
  double Mass;
  double concentration;
  double Rmis;
  int delta;
  int miscentering;
  int averaging;
  int single_miscentering;
  int Nbins;
  double R_bin_min;
  double R_bin_max;
}interface_parameters;
#endif

#ifndef WRAPPER_OUTPUT
#define WRAPPER_OUTPUT
typedef struct wrapper_output{
  double*R;
  double*xi_1halo;
  double*xi_mm;
  double*xi_lin;
  double*xi_2halo;
  double*xi_hm;
  double*sigma;
  double*delta_sigma;
  double*Rbins;
  double*ave_delta_sigma;
  double*bias;
  double*nu;
  double*sigma_mis;
  double*delta_sigma_mis;
  double*miscentered_sigma;
  double*miscentered_delta_sigma;
  double*ave_miscentered_delta_sigma;
  double*ave_delta_sigma_mis;
}wrapper_output;
#endif

int interface(double*k_lin,double*P_lin,int Nk_lin,
	      double*k,double*P,int Nk,
	      int NR,double Rmin,double Rmax,
	      cosmology cosmo, interface_parameters*params,
	      wrapper_output*outputs);

int python_interface(double*k_lin,double*P_lin,int Nk_lin,
		     double*k,double*P,int Nk,
		     int NR,double Rmin,double Rmax,
		     double h,double om,double ode,double ok,
		     double Mass, double concentration,
		     double Rmis, int delta,
		     int miscentering,
		     int averaging, int single_miscentering,
		     int Nbins,
		     double R_bin_min, double R_bin_max,
		     double*R,double*xi_1halo,double*xi_mm,double*xi_lin,
		     double*xi_2halo,double*xi_hm,double*sigma,
		     double*delta_sigma,double*Rbins,
		     double*ave_delta_sigma,double*bias,
		     double*nu,double*sigma_mis,
		     double*delta_sigma_mis,
		     double*miscentered_sigma,
		     double*miscentered_delta_sigma,
		     double*ave_miscentered_delta_sigma,
		     double*ave_delta_sigma_mis);
