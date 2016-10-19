#include "gsl/gsl_spline.h"
#include "gsl/gsl_integration.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

int calc_xi_mm(double *R,int NR,double*k,double*P,
	       int Nk,double*xi,double*err,int N,double h);
