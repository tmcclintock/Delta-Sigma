#include "../sigma_mis/sigma_mis_at_r.h"
#include "delta_sigma_mis_at_r.h"
#include "../constants/constants.h"
#include "../cosmology/cosmology.h"

int calc_delta_sigma_mis(double*Rp,double Mass,double concentration,
			 int delta,double Rmis,double*R,
			 double*sigma,double*sigma_mis,
			 int NR,double*delta_sigma_mis,
			 double*err,cosmology cosmo);
