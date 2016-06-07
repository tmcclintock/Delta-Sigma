#include "miscentered_sigma_r.h"

int calc_miscentered_sigma_r(double*Rp,double Mass,double concentration,
			     int delta,double Rmis,double*R,double*sigma_r,
			     int NR,double*mis_sigma_r,double*err,
			     cosmology cosmo){
  int i, status=0;
  double time=omp_get_wtime();

#pragma omp parallel shared(R,sigma_r,NR,mis_sigma_r,err,status)
#pragma omp for
  for(i = 0; i < NR; i++){
    status |= calc_miscentered_sigma_r_at_r(Rp[i],Mass,concentration,delta,Rmis,R,sigma_r,NR,&mis_sigma_r[i],&err[i],cosmo);
  }

  return 0;
}
