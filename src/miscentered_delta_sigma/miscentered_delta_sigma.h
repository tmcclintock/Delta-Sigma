#include "../miscentered_sigma_r/miscentered_sigma_r_at_r.h"
#include "miscentered_delta_sigma_at_r.h"
#include "../cosmology/cosmology.h"
#include <omp.h>

int calc_miscentered_delta_sigma(double*Rp,double Mass,double concentration,
				 int delta,double Rmis,double*R,
				 double*sigma_r,double*miscentered_sigma_r,
				 int NR,double*miscentered_delta_sigma,
				 double*err,cosmology cosmo);
