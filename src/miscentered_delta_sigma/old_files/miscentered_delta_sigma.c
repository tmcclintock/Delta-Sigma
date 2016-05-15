#include "miscentered_delta_sigma.h"

int calc_miscentered_delta_sigma(double*Rp,double Mass,double concentration,
				 int delta,double Rmis,double*R,
				 double*sigma_r,double*miscentered_sigma_r,
				 int NR,double*miscentered_delta_sigma,
				 double*err,cosmology cosmo){
  int i, status=0;
#pragma omp parallel shared(R,sigma_r,NR,miscentered_sigma_r,miscentered_delta_sigma,err,status)
#pragma omp for
  for(i = 0; i < NR; i++){
    status |= calc_miscentered_delta_sigma_at_r(Rp[i],Mass,concentration,delta,
						Rmis,R,sigma_r,
						miscentered_sigma_r,
						NR,&miscentered_delta_sigma[i],
						&err[i],cosmo);
  }
}
