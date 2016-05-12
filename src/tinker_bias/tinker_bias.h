#include "tinker_bias_at_M.h"
#include <omp.h>
#include "../cosmology/cosmology.h"

int calc_tinker_bias(double*M,int NM,double*k,double*P,int N,double*b,double*nu,int delta,cosmology cosmo);
