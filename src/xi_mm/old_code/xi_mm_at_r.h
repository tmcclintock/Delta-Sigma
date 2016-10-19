#include "gsl/gsl_spline.h"
#include "gsl/gsl_integration.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int calc_xi_mm_at_r(double R,double*k,double*P,
		    int N,double*xi,double*err);
