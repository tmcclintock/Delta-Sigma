#include "wrapper.h"

int interface(double*k_lin,double*P_lin,int Nk_lin,
	      double*k,double*P,int Nk,
	      int NR,double Rmin,double Rmax,
	      cosmology cosmo, interface_parameters*params,
	      wrapper_output*outputs){
  int i;

  double*R=outputs->R;
  double*xi_1halo=outputs->xi_1halo;
  double*xi_mm=outputs->xi_mm;
  double*xi_lin=outputs->xi_lin;
  double*xi_2halo=outputs->xi_2halo;
  double*xi_hm=outputs->xi_hm;
  double*sigma=outputs->sigma;
  double*delta_sigma=outputs->delta_sigma;
  double*Rbins=outputs->Rbins;
  double*ave_delta_sigma=outputs->ave_delta_sigma;
  double*bias=outputs->bias;
  double*nu=outputs->nu;
  double*sigma_mis=outputs->sigma_mis;
  double*delta_sigma_mis=outputs->delta_sigma_mis;
  double*miscentered_sigma=outputs->miscentered_sigma;
  double*miscentered_delta_sigma=outputs->miscentered_delta_sigma;
  double*ave_miscentered_delta_sigma=outputs->ave_miscentered_delta_sigma;
  double*ave_delta_sigma_mis=outputs->ave_delta_sigma_mis;

  //Used to hold integration errors for some routines
  double*err=(double*)malloc(NR*sizeof(double));

  double Mass=params->Mass;
  double concentration=params->concentration;
  double Rmis=params->Rmis;
  int delta=params->delta;

  int Nbins=params->Nbins;
  double R_bin_min=params->R_bin_min;
  double R_bin_max=params->R_bin_max;

  int miscentering=params->miscentering;
  int averaging=params->averaging;
  int single_miscentering=params->single_miscentering;

  double dlR = (log(Rmax) - log(Rmin))/(float)(NR-1.);
  for(i = 0; i < NR; i++){
    R[i] = exp(log(Rmin) + i*dlR);
  }

  calc_xi_nfw(R,NR,Mass,concentration,delta,xi_1halo,cosmo);
  /*
    These two variables are for the hankel transformation
    They can be changed by hand, depending on how wacky the
    the power spectrum looks. It is difficult to test
    what combination works best for what scenarios.
  */
  int Nevals = 200;
  double h = 0.005;
  calc_xi_mm(R,NR,k,P,Nk,xi_mm,err,Nevals,h);
  calc_xi_mm(R,NR,k_lin,P_lin,Nk_lin,xi_lin,err,Nevals,h);
  calc_tinker_bias(&Mass,1,k_lin,P_lin,Nk_lin,bias,nu,delta,cosmo);
  calc_xi_2halo(NR,xi_mm,*bias,xi_2halo);
  calc_xi_hm(NR,Mass,xi_1halo,xi_2halo,xi_hm);
  calc_sigma(R,Mass,concentration,delta,R,xi_hm,NR,sigma,err,cosmo);
  calc_delta_sigma(R,Mass,concentration,delta,R,sigma,NR,delta_sigma,err,cosmo);
  if(averaging){
    calc_ave_delta_sigma(R,NR,delta_sigma,Nbins,R_bin_min,R_bin_max,Rbins,ave_delta_sigma);
  }

  if(single_miscentering){
    calc_sigma_mis(R,Mass,concentration,delta,Rmis,R,sigma,NR,sigma_mis,err,cosmo);
    calc_delta_sigma_mis(R,Mass,concentration,delta,Rmis,R,sigma,sigma_mis,NR,delta_sigma_mis,err,cosmo);
    if(averaging){
      calc_ave_miscentered_delta_sigma(R,NR,delta_sigma_mis,Nbins,R_bin_min,R_bin_max,Rbins,ave_delta_sigma_mis);
    }
  }

  if(miscentering){
    calc_miscentered_sigma(R,Mass,concentration,delta,Rmis,R,sigma,NR,miscentered_sigma,err,cosmo);
    calc_miscentered_delta_sigma(R,Mass,concentration,delta,Rmis,R,sigma,miscentered_sigma,NR,miscentered_delta_sigma,err,cosmo);
    if(averaging){
      calc_ave_miscentered_delta_sigma(R,NR,miscentered_delta_sigma,Nbins,R_bin_min,R_bin_max,Rbins,ave_miscentered_delta_sigma);
    }
  }

  free(err);
  return 0;
}

int python_interface(double*k_lin,double*P_lin,int Nk_lin,
		     double*k_nl,double*P_nl,int Nk_nl,
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
		     double*ave_delta_sigma_mis){

  cosmology*cosmo = (cosmology*)malloc(sizeof(cosmology));
  cosmo->H0 = h*100.;
  cosmo->h = h;
  cosmo->om = om;
  cosmo->ode = ode;
  cosmo->ok = ok;

  interface_parameters*params=
    (interface_parameters*)malloc(sizeof(interface_parameters));
  params->Mass=Mass;
  params->concentration=concentration;
  params->delta=delta;
  params->Rmis=Rmis;
  params->miscentering=miscentering;//1 is true
  params->averaging=averaging; //1 is true
  params->single_miscentering=single_miscentering; //1 is true
  params->Nbins=Nbins;
  params->R_bin_min=R_bin_min;
  params->R_bin_max=R_bin_max;

  wrapper_output*outputs=(wrapper_output*)malloc(sizeof(wrapper_output));
  outputs->R=R;
  outputs->xi_1halo=xi_1halo;
  outputs->xi_mm=xi_mm;
  outputs->xi_lin=xi_lin;
  outputs->xi_2halo=xi_2halo;
  outputs->xi_hm=xi_hm;
  outputs->sigma=sigma;
  outputs->delta_sigma=delta_sigma;
  outputs->bias=bias;
  outputs->nu=nu;
  outputs->sigma_mis=sigma_mis;
  outputs->delta_sigma_mis=delta_sigma_mis;
  outputs->miscentered_sigma=miscentered_sigma;
  outputs->miscentered_delta_sigma=miscentered_delta_sigma;
  outputs->Rbins=Rbins;
  outputs->ave_delta_sigma=ave_delta_sigma;
  outputs->ave_miscentered_delta_sigma=ave_miscentered_delta_sigma;
  outputs->ave_delta_sigma_mis=ave_delta_sigma_mis;

  interface(k_lin,P_lin,Nk_lin,k_nl,P_nl,Nk_nl,
	    NR,Rmin,Rmax,*cosmo,params,outputs);

  free(cosmo);
  free(params);
  free(outputs);
  return 0;
}
		     
