#include "gsl/gsl_spline.h"
#include "../constants/constants.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int calc_xi_mm(double *R,int NR,double*k,double*P,int Nk,double*xi,double*err,int N,double h);
