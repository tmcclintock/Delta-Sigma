#include "gsl/gsl_spline.h"
#include "gsl/gsl_integration.h"
#include "../cosmology/cosmology.h"
#include "../xi_nfw/xi_nfw.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int calc_sigma_at_r(double Rp,double Mass,double concentration,
		    int delta,double*R,double*xi,int NR,double*sigma,
		    double*err,cosmology cosmo);
