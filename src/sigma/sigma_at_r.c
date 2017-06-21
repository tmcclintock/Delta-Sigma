#include "sigma_at_r.h"

/* These are the parameters passed into the integrand.
   A spline and accelerator for interpolation and the radius 
   we are evaluating xi at.*/
typedef struct integrand_params{
  gsl_spline*spline;
  gsl_interp_accel*acc;
  gsl_integration_workspace * workspace;
  double rperp;
  double lrmin;
  double lrmax;
  cosmology cosmo;
  double Mass;
  double concentration;
  int delta;
}integrand_params;

static int do_integral(double*sigmar_r,double*err,integrand_params*params);
static double integrand_small_scales(double lrz, void*params);
static double integrand_medium_scales(double lrz, void*params);

int calc_sigma_at_r(double Rp,double Mass,double concentration,
		    int delta,double*R,double*xi,int NR,
		    double*sigma,double*err,cosmology cosmo){

  double h = cosmo.h, om = cosmo.om;
  double H0 = h*100.;
  double rhom = om*rhomconst*1e-12; //SM h^2/pc^2/Mpc
  //note, this isn't ACTUALLY rhom, but it contains the integration constants
  gsl_spline*spline = gsl_spline_alloc(gsl_interp_cspline,NR);
  gsl_spline_init(spline,R,xi,NR);
  gsl_interp_accel*acc= gsl_interp_accel_alloc();
  gsl_integration_workspace * workspace
    = gsl_integration_workspace_alloc(workspace_size);
  
  integrand_params*params=malloc(sizeof(integrand_params));
  params->acc=acc;
  params->spline=spline;
  params->workspace=workspace;
  params->rperp=Rp;
  params->lrmin=log(R[0]);
  params->lrmax=log(R[NR-1]);
  params->cosmo=cosmo;
  params->Mass=Mass;
  params->concentration=concentration;
  params->delta=delta;
  
  do_integral(sigma,err,params);
  *sigma *= rhom*2;
  *err *= rhom*2;

  gsl_spline_free(spline),gsl_interp_accel_free(acc);
  gsl_integration_workspace_free(workspace);
  free(params);
  return 0;
}

int do_integral(double*sigma,double*err,integrand_params*params){
  gsl_function F;
  F.params=params;
  
  double lrmin = params->lrmin;
  double lrmax = params->lrmax;
  double rmax = exp(lrmax);
  double rp = params->rperp;
  double lrz_max = log(sqrt(rmax*rmax - rp*rp));
  gsl_integration_workspace*workspace=params->workspace;

  double result1=0,abserr1=0;
  double result2=0,abserr2=0;
  double result3=0,abserr3=0;

  F.function=&integrand_small_scales;
  int status = gsl_integration_qag(&F,lrmin-10,lrmin,TOL,TOL/10.,workspace_size,6,workspace,&result1,&abserr1);

  F.function=&integrand_medium_scales;
  if (lrz_max > 0)
    status |= gsl_integration_qag(&F,lrmin,lrz_max,TOL,TOL/10.,workspace_size,6,workspace,&result2,&abserr2);
  
  *sigma = result1+result2;
  *err= abserr1+abserr2;
  return status;
}

double integrand_small_scales(double lrz, void*params){
  double rz = exp(lrz);
  integrand_params pars=*(integrand_params *)params;
  double rp = pars.rperp;
  cosmology cosmo = pars.cosmo;
  double arg = sqrt(rz*rz+rp*rp);
  double Mass = pars.Mass;
  double concentration = pars.concentration;
  int delta = pars.delta;
  double answer = 0;
  calc_xi_nfw(&arg,1,Mass,concentration,delta,&answer,cosmo);
  return rz*answer;
}

double integrand_medium_scales(double lrz, void*params){
  double rz = exp(lrz);
  integrand_params pars=*(integrand_params*)params;
  gsl_spline*spline = pars.spline;//Xi_hm(R) spline
  gsl_interp_accel*acc = pars.acc;
  double rp = pars.rperp;
  return rz*gsl_spline_eval(spline,sqrt(rz*rz+rp*rp),acc);
}
