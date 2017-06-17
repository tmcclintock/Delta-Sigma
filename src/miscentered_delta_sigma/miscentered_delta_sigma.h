#include "../miscentered_sigma/miscentered_sigma_at_r.h"
#include "miscentered_delta_sigma_at_r.h"
#include "../constants/constants.h"
#include "../cosmology/cosmology.h"

int calc_miscentered_delta_sigma(double*Rp,double Mass,double concentration,
				 int delta,double Rmis,double*R,
				 double*sigma,double*miscentered_sigma,
				 int NR,double*miscentered_delta_sigma,
				 double*err,cosmology cosmo);
