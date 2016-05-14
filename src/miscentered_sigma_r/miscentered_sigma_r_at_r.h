#include "gsl/gsl_spline.h"
#include "gsl/gsl_integration.h"
#include "../cosmology/cosmology.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int calc_miscentered_sigma_r_at_r(double Rp,double Mass,double concentration,int delta,double Rmis,double*R,double*sigma_r,int NR,double*mis_sigma_r,double*err,cosmology cosmo);
