#include "wrapper.h"

int interface(double*k,double*P,int Nk,int NR,double Rmin,double Rmax,
	      cosmology cosmo, interface_parameters*params,
	      wrapper_output*outputs){
  int i;

  double*R=outputs->R;
  double*xi_1halo=outputs->xi_1halo;
  double*xi_mm=outputs->xi_mm;
  double*xi_2halo=outputs->xi_2halo;
  double*xi_hm=outputs->xi_hm;
  double*sigma_r=outputs->sigma_r;
  double*delta_sigma=outputs->delta_sigma;
  double*bias=outputs->bias;
  double*nu=outputs->nu;
  
  //Used to hold integration errors for some routines
  double*err=(double*)malloc(NR*sizeof(double));

  double Mass=params->Mass;
  double concentration=params->concentration;
  double lnc=params->lnc;
  double fmis=params->fmis;
  int delta=params->delta;

  int*flow_control=params->flow_control;

  for(i = 0; i < NR; i++){
    double dlR = (log(Rmax) - log(Rmin))/(float)(NR-1.);
    R[i] = exp(log(Rmin) + i*dlR);
  }

  calc_xi_nfw(R,NR,Mass,concentration,delta,xi_1halo,cosmo);
  calc_xi_mm(R,NR,k,P,Nk,xi_mm,err);
  calc_tinker_bias(&Mass,1,k,P,Nk,bias,nu,delta,cosmo);
  calc_xi_2halo(NR,xi_mm,*bias,xi_2halo);
  calc_xi_hm(NR,Mass,xi_1halo,xi_2halo,xi_hm);
  calc_sigma_r(R,Mass,concentration,delta,R,xi_hm,NR,sigma_r,err,cosmo);
  calc_delta_sigma(R,Mass,concentration,delta,R,sigma_r,NR,delta_sigma,err,cosmo);

  free(err);
  return 0;
}
