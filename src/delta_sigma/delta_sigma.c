#include "delta_sigma.h"

int calc_delta_sigma(double*Rp,double Mass,double concentration,int delta,double*R,double*sigma_r,int NR,double*delta_sigma,double*err,cosmology cosmo){
    int i;
#pragma omp parallel shared(Rp,R,sigma_r,NR,delta_sigma,err) private(i)
#pragma omp for
    for(i = 0; i < NR; i++){
      calc_delta_sigma_at_r(Rp[i],Mass,concentration,delta,R,sigma_r,NR,&delta_sigma[i],&err[i],cosmo);
    }
    return 0;
}
