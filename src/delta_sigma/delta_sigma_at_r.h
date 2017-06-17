#include "gsl/gsl_spline.h"
#include "gsl/gsl_integration.h"
#include "../constants/constants.h"
#include "../cosmology/cosmology.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int calc_delta_sigma_at_r(double Rp,double Mass,double concentration,
			  int delta,double*R,double*sigma,int NR,
			  double*delta_sigma,double*err,cosmology cosmo);
