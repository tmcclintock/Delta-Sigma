#include "sigma_mis_at_r.h"
#include "../cosmology/cosmology.h"

int calc_sigma_mis(double*Rp,double Mass,double concentration,
		   int delta,double Rmis,double*R,double*sigma,
		   int NR,double*sigma_mis,
		   double*err,cosmology cosmo);
