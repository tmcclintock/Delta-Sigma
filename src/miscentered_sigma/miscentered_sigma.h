#include "miscentered_sigma_at_r.h"
#include "../cosmology/cosmology.h"
#include <omp.h>

int calc_miscentered_sigma(double*Rp,double Mass,double concentration,
			   int delta,double Rmis,double*R,double*sigma,
			   int NR,double*mis_sigma,
			   double*err,cosmology cosmo);
