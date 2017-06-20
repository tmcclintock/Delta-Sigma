#include "sigma_at_r.h"
#include "../constants/constants.h"
#include "../cosmology/cosmology.h"
#include "../xi_nfw/xi_nfw.h"

int calc_sigma(double*Rp,double Mass,double concentration,
	       int delta,double*R,double*xi,int NR,double*sigma,
	       double*err,cosmology cosmo);
