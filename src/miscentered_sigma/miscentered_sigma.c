#include "miscentered_sigma.h"

int calc_miscentered_sigma(double*Rp,double Mass,double concentration,
			   int delta,double Rmis,double*R,double*sigma,
			   int NR,double*mis_sigma,
			   double*err,
			   cosmology cosmo){
  int i, status=0;
  double time=omp_get_wtime();

#pragma omp parallel for shared(R,sigma,NR,mis_sigma,err,status)
  for(i = 0; i < NR; i++){
    status |= calc_miscentered_sigma_at_r(Rp[i],Mass,concentration,delta,Rmis,R,sigma,NR,&mis_sigma[i],&err[i],cosmo);
  }

  return 0;
}
