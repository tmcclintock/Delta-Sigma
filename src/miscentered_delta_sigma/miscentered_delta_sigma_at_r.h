#include "gsl/gsl_spline.h"
#include "gsl/gsl_integration.h"
#include "../cosmology/cosmology.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int calc_miscentered_delta_sigma_at_r(double Rp,double Mass,double concentration,int delta,double Rmis,double*R,double*sigma_r,double*mis_sigma_r,double sigma_less_rmin,int NR,double*mis_delta_sigma,double*err,cosmology cosmo);
