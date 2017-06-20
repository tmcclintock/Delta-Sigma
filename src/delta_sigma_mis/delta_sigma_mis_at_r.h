#include "gsl/gsl_spline.h"
#include "gsl/gsl_integration.h"
#include "../constants/constants.h"
#include "../cosmology/cosmology.h"

int calc_delta_sigma_mis_at_r(double Rp,double Mass,double concentration,
			      int delta,double Rmis,double*R,double*sigma,
			      double*sigma_mis,double sigma_less_rmin,
			      int NR,double*delta_sigma_mis,double*err,
			      cosmology cosmo);
