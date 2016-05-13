#include "sigma_r.h"

int calc_sigma_r(double*Rp,double Mass,double concentration,int delta,double*R,double*xi,int NR,double*sigma_r,double*err,cosmology cosmo){
    int i;
#pragma omp parallel shared(R,xi,NR,sigma_r,err)
#pragma omp for
    for(i = 0; i < NR; i++){
      calc_sigma_r_at_r(Rp[i],Mass,concentration,delta,R,xi,NR,&sigma_r[i],&err[i],cosmo);
      //printf("%d %e %e %e\n",i,R[i],sigma_r[i],err[i]);
    }
    return 0;
}
