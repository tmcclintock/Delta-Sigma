#include "gsl/gsl_spline.h"
#include "gsl/gsl_integration.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../cosmology/cosmology.h"

int tinker_bias_at_M(double M,double*k,double*P,int N,double*bias,double*nu,int delta,cosmology cosmo);
