#include "tinker_bias.h"

int calc_tinker_bias(double*M,int NM,double*k,double*P,int N,double*b,double*nu,int delta,cosmology cosmo){
  int i;
  for(i = 0; i< NM; i++){
    tinker_bias_at_M(M[i],k,P,N,&b[i],&nu[i],delta,cosmo);
  }
  return 0;
}
