#include "miscentered_ave_delta_sigma_in_bin.h"
#include <math.h>
#include <omp.h>

int calc_miscentered_ave_delta_sigma(double*R,int NR,
				     double*miscentered_delta_sigma,
				     int Nbins,double R_bin_min,
				     double R_bin_max,
				     double*Rbins,
				     double*ave_delta_sigma);