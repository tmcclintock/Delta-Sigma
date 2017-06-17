#include "miscentered_delta_sigma_at_r.h"

/* These are the parameters passed into the integrand.
   A spline and accelerator for interpolation and the radius 
   we are evaluating mis_sigma at.*/
typedef struct integrand_params{
  gsl_spline*spline;
  gsl_interp_accel*acc;
  gsl_integration_workspace*workspace;
  double lrperp;
  double lrmin;
  double lrmax;
}integrand_params;

static int do_integral(double*miscentered_delta_sigma,double*err,integrand_params*params);

static double integrand(double lR, void*params){
  double R = exp(lR);
  integrand_params pars=*(integrand_params *)params;
  gsl_spline*spline = pars.spline;//Sigma(R) spline
  gsl_interp_accel*acc = pars.acc;
  return R*R*gsl_spline_eval(spline,R,acc);
}

int calc_miscentered_delta_sigma_at_r(double Rp,double Mass,
				      double concentration,int delta,
				      double Rmis,double*R,double*sigma,
				      double*miscentered_sigma,double inner_result,
				      int NR,double*miscentered_delta_sigma,
				      double*err,cosmology cosmo){
  int status = 0;

  gsl_spline*spline = gsl_spline_alloc(gsl_interp_cspline,NR);
  gsl_spline_init(spline,R,miscentered_sigma,NR);
  gsl_interp_accel*acc= gsl_interp_accel_alloc();
  gsl_integration_workspace * workspace
    = gsl_integration_workspace_alloc(workspace_size);
  integrand_params*params=malloc(sizeof(integrand_params));
  params->acc=acc;
  params->spline=spline;
  params->workspace=workspace;
  params->lrperp=log(Rp);
  params->lrmin=log(R[0]);
  params->lrmax=log(R[NR-1]);

  status |= do_integral(miscentered_delta_sigma,err,params);
  *miscentered_delta_sigma += inner_result;
  *miscentered_delta_sigma *= 2./Rp/Rp;
  *miscentered_delta_sigma -= gsl_spline_eval(spline,Rp,acc);
  *err *= 2./Rp/Rp;

  gsl_spline_free(spline),gsl_interp_accel_free(acc);
  gsl_integration_workspace_free(workspace);
  free(params);
  return status;
}

int do_integral(double*miscentered_delta_sigma,double*err,integrand_params*params){
  int status = 0;

  gsl_function F;
  F.params=params;

  double lrmin = params->lrmin;
  double lrmax = params->lrmax;
  double lRp = params->lrperp;
  gsl_integration_workspace*workspace=params->workspace;

  double result=0,abserr=0;

  F.function=&integrand;
  status |= gsl_integration_qag(&F,lrmin,lRp,TOL,TOL/10.,workspace_size,6,workspace,&result,&abserr);

  *miscentered_delta_sigma = result;
  *err = abserr;
  
  return status;
}

