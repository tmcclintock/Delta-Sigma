#include "gsl/gsl_spline.h"
#include "gsl/gsl_integration.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int calc_ave_delta_sigma_in_bin(double*R,int NR,double*delta_sigma,
				double lRlow,double lRhigh,
				double*ave_delta_sigma);