#include "xi_2halo.h"

int calc_xi_2halo(int NR,double*xi_mm,double bias,double*xi_2halo){
  int i;
  for(i = 0; i < NR; i++){
    xi_2halo[i] = bias * xi_mm[i];
  }
}
