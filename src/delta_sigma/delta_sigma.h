#include "delta_sigma_at_r.h"
#include "../cosmology/cosmology.h"
#include "../sigma/sigma.h"
#include <omp.h>

int calc_delta_sigma(double*Rp,double Mass,double concentration,
		     int delta,double*R,double*sigma,int NR,
		     double*delta_sigma,double*err,cosmology cosmo);
