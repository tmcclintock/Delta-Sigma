#include "sigma_r.h"

int calc_sigma_r(double*Rp,double Mass,double concentration,int delta,double*R,double*xi,int NR,double*sigma_r,double*err,cosmology cosmo){
  int i,status=0;

#pragma omp parallel shared(R,xi,NR,sigma_r,err,status)
#pragma omp for
    for(i = 0; i < NR; i++)
      status |= calc_sigma_r_at_r(Rp[i],Mass,concentration,delta,R,xi,NR,&sigma_r[i],&err[i],cosmo);
    
    return status;
}
