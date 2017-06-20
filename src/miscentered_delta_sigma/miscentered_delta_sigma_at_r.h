#include "gsl/gsl_spline.h"
#include "gsl/gsl_integration.h"
#include "../constants/constants.h"
#include "../cosmology/cosmology.h"

int calc_miscentered_delta_sigma_at_r(double Rp,double Mass,
				      double concentration,int delta,
				      double Rmis,double*R,double*sigma,
				      double*mis_sigma,
				      double sigma_less_rmin,int NR,
				      double*mis_delta_sigma,
				      double*err,cosmology cosmo);
