#include "sigma_mis.h"

int calc_sigma_mis(double*Rp,double Mass,double concentration,
		   int delta,double Rmis,double*R,double*sigma,
		   int NR,double*sigma_mis,
		   double*err,
		   cosmology cosmo){
  int i, status=0;
  for(i = 0; i < NR; i++){
    status |= calc_sigma_mis_at_r(Rp[i],Mass,concentration,delta,Rmis,
				  R,sigma,NR,&sigma_mis[i],&err[i],cosmo);
  }

  return 0;
}
