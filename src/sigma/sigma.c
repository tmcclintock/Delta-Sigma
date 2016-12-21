#include "sigma.h"

int calc_sigma(double*Rp,double Mass,double concentration,
	       int delta,double*R,double*xi,int NR,double*sigma,
	       double*err,cosmology cosmo){
  int i,status=0;

#pragma omp parallel shared(R,xi,NR,sigma,err,status)
#pragma omp for
    for(i = 0; i < NR; i++)
      status |= calc_sigma_at_r(Rp[i],Mass,concentration,delta,R,xi,NR,&sigma[i],&err[i],cosmo);
    
    return status;
}
