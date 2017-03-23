#include "sigma_at_r.h"
#include "../cosmology/cosmology.h"
#include "../xi_nfw/xi_nfw.h"
#include <omp.h>

int calc_sigma(double*Rp,double Mass,double concentration,
	       int delta,double*R,double*xi,int NR,double*sigma,
	       double*err,cosmology cosmo);