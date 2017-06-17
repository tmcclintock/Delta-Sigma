#include "xi_hm.h"

int calc_xi_hm(int NR, double*xi_1h, double*xi_2h, double*xi_hm){
  int i;
  for(i = 0; i < NR; i++){
    if(xi_1h[i] >= xi_2h[i])
      xi_hm[i] = xi_1h[i];
    else
      xi_hm[i] = xi_2h[i];
  }
  return 0;
}
