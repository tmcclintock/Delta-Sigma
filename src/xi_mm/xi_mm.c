#include "xi_mm.h"

int calc_xi_mm(double *R,int NR,double*k,double*P,int N,double*xi,double*err){
  int i;
  for(i = 0; i < NR; i++){
    calc_xi_mm_at_r(R[i],k,P,N,xi+i,err+i);
  }
  return 0;
}
