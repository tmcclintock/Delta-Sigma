#include "delta_sigma_at_r.h"

#define TOL 1e-8
#define workspace_size 8000
#define PI 3.141592653589793

//These are physical constants
#define G 4.517e-48//Newton's G in Mpc^3/s^2/Solar Mass
#define Mpcperkm 3.241e-20//Mpc/km used to convert H0 to per seconds

/* These are the parameters passed into the integrand.
   A spline and accelerator for interpolation and the radius 
   we are evaluating sigma_r at.*/
typedef struct integrand_params{
  gsl_spline*spline;
  gsl_interp_accel*acc;
  gsl_integration_workspace * workspace;
  double lrperp;
  double lrmin;
  double lrmax;
  cosmology cosmo;
  double Mass;
  double concentration;
  int delta;
}integrand_params;

static int do_integral(double*sigmar_r,double*err,integrand_params*params);

static double integrand(double r_z, void*params);

static double sigma_r_1halo_analytic(double R,double Mass,double concentration,double om,double H0,double delta);

int calc_delta_sigma_at_r(double Rp,double Mass,double concentration
		      ,int delta,double*R,double*sigma_r,int NR,
		      double*delta_sigma,double*err,cosmology cosmo){

  gsl_spline*spline = gsl_spline_alloc(gsl_interp_cspline,NR);
  gsl_spline_init(spline,R,sigma_r,NR);
  gsl_interp_accel*acc= gsl_interp_accel_alloc();
  gsl_integration_workspace * workspace
    = gsl_integration_workspace_alloc(workspace_size);
  
  integrand_params*params=malloc(sizeof(integrand_params));
  params->acc=acc;
  params->spline=spline;
  params->workspace=workspace;
  params->rperp=log(Rp);
  params->lrmin=log(R[0]);
  params->lrmax=log(R[NR-1]);
  params->cosmo=cosmo;
  params->Mass=Mass;
  params->concentration=concentration;
  params->delta=delta;
  
  do_integral(delta_sigma,err,params);
  *delta_sigma *= 2./Rp/Rp;
  *delta_sigma -= gsl_spline_eval(spline,Rp,acc);
  *err *= 2./Rp/Rp;

  gsl_spline_free(spline),gsl_interp_accel_free(acc);
  gsl_integration_workspace_free(workspace);
  free(params);
  return 0;
}

static int do_integral(double*delta_sigma,double*err,integrand_params*params){
  gsl_function F;
  F.function=&integrand;
  F.params=params;
  
  double lrmin = params->lrmin;
  double lrmax = params->lrmax;
  double lRp = params->lrperp;
  gsl_integration_workspace*workspace=params->workspace;

  double result,abserr;
  int status = gsl_integration_qag(&F,lrmin-10,lRp,TOL,TOL/10.,
				   workspace_size,6,workspace,&result,&abserr);
  *delta_sigma = result;
  *err= abserr;
  return status;
}

static double integrand(double lR, void*params){
  double R = exp(lR);
  integrand_params pars=*(integrand_params *)params;
  gsl_spline*spline = pars.spline;//Sigma_R(R) spline
  gsl_interp_accel*acc = pars.acc;
  double rmin = exp(pars.lrmin);
  double answer = 0;
  cosmology cosmo = pars.cosmo;
  double Mass = pars.Mass;
  double concentration = pars.concentration;
  int delta = pars.delta;


  if(R<rmin){//fix
    calc_sigma_r(&R,1,Mass,concentration,delta,&answer,cosmo);
    return R*R*answer;
  }else
    return R*R*gsl_spline_eval(spline,arg,acc);
}

static double sigma_r_1halo_analytic(double R,double rscale,double Mass,double om,double H0,double delta){

}
