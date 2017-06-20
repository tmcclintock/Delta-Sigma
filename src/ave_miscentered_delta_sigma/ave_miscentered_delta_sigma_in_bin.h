#include "gsl/gsl_spline.h"
#include "gsl/gsl_integration.h"
#include "../constants/constants.h"

int calc_ave_miscentered_delta_sigma_in_bin(double*R,int NR,
					    double*miscentered_delta_sigma,
					    double lRlow,double lRhigh,
					    double*ave_miscentered_delta_sigma);
