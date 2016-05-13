#include "gsl/gsl_spline.h"
#include "gsl/gsl_integration.h"
#include "../cosmology/cosmology.h"
#include "../sigma_r/sigma_r.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int calc_delta_sigma_at_r(double Rp,double Mass,double concentration,int delta,double*R,double*sigma_r,int NR,double*delta_sigma,double*err,cosmology cosmo);
