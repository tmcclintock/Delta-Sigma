#include "delta_sigma.h"

int calc_delta_sigma(double*Rp,double Mass,double concentration,
		     int delta,double*R,double*sigma,int NR,
		     double*delta_sigma,double*err,cosmology cosmo){
  int i;
  for(i = 0; i < NR; i++){
    calc_delta_sigma_at_r(Rp[i],Mass,concentration,delta,R,sigma,NR,&delta_sigma[i],&err[i],cosmo);
  }
  return 0;
}
