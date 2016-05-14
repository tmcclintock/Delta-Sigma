#include "miscentered_sigma_r.h"

int calc_miscentered_sigma_r(double*Rp,double Mass,double concentration,
			     int delta,double Rmis,double*R,double*xi,
			     int NR,double*sigma_r,double*err,
			     cosmology cosmo){
  int i, status=0;
#pragma omp parallel shared(R,xi,NR,sigma_r,err,status)
#pragma omp for
  for(i = 0; i < NR; i++){
    status |= calc_miscentered_sigma_r_at_r(Rp[i],Mass,concentration,delta,Rmis,R,xi,NR,&sigma_r[i],&err[i],cosmo);
  }

  return 0;
}
