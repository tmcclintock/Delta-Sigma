#include "xi_hm.h"

int calc_xi_hm(int NR,double Mass,double*xi_1h,double*xi_mm,double bias,double*xi_hm){
  int i;
  for(i = 0; i < NR; i++){
    if(xi_1h[i] >= xi_mm[i]*bias)
      xi_hm[i] = xi_1h[i];
    else
      xi_hm[i] = xi_mm[i]*bias;
  }
  return 0;
}
