#include "miscentered_sigma_r_at_r.h"
#include "../cosmology/cosmology.h"
#include <omp.h>

int calc_miscentered_sigma_r(double*Rp,double Mass,double concentration,int delta,double Rmis,double*R,double*sigma_r,int NR,double*mis_sigma_r,double*err,cosmology cosmo);
