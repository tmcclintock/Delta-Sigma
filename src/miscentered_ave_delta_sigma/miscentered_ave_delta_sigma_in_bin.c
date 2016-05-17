#include "miscentered_ave_delta_sigma_in_bin.h"

#define TOL 1e-8
#define workspace_size 8000
#define PI 3.141592653589793

/* These are the parameters passed into the integrand.
   A spline and accelerator for interpolation
   of delta_sigma.*/
typedef struct integrand_params{
  gsl_spline*spline;
  gsl_interp_accel*acc;
  gsl_integration_workspace*workspace;
}integrand_params;

static int do_integral(double*miscentered_ave_delta_sigma,
		       double lRlow,double lRhigh,
		       integrand_params*params);

static double integrand(double lR, void*params);

int calc_miscentered_ave_delta_sigma_in_bin(double*R,int NR,
					    double*miscentered_delta_sigma,
					    double lRlow,double lRhigh,
					    double*miscentered_ave_delta_sigma){
  int status = 0;

  gsl_spline*spline = gsl_spline_alloc(gsl_interp_cspline,NR);
  gsl_spline_init(spline,R,miscentered_delta_sigma,NR);
  gsl_interp_accel*acc= gsl_interp_accel_alloc();
  gsl_integration_workspace * workspace
    = gsl_integration_workspace_alloc(workspace_size);

  integrand_params*params=malloc(sizeof(integrand_params));
  params->acc=acc;
  params->spline=spline;
  params->workspace=workspace;

  double Rlow=exp(lRlow),Rhigh=exp(lRhigh);

  do_integral(miscentered_ave_delta_sigma,lRlow,lRhigh,params);
  *miscentered_ave_delta_sigma *= 2./(Rhigh*Rhigh-Rlow*Rlow);

  gsl_spline_free(spline),gsl_interp_accel_free(acc);
  gsl_integration_workspace_free(workspace);
  free(params);
  return 0;
}

int do_integral(double*miscentered_ave_delta_sigma,double lRlow,double lRhigh,
		integrand_params*params){
  int status = 0;

  gsl_function F;
  F.function=&integrand;
  F.params=params;

  gsl_integration_workspace*workspace=params->workspace;

  double result,abserr;

  status |= gsl_integration_qag(&F,lRlow,lRhigh,TOL,TOL/10.,workspace_size,6,
				workspace,&result,&abserr);

  *miscentered_ave_delta_sigma = result;

  return status;
}

double integrand(double lR, void*params){
  double R = exp(lR);
  integrand_params pars=*(integrand_params *)params;
  gsl_spline*spline = pars.spline;//Mis_Delta_Sigma(R) spline
  gsl_interp_accel*acc = pars.acc;
  return R*R*gsl_spline_eval(spline,R,acc);
}
