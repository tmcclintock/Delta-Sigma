#include "sigma_r_at_r.h"

#define TOL 1e-8
#define workspace_size 8000
#define PI 3.141592653589793

//These are physical constants
#define G 4.517e-48//Newton's G in Mpc^3/s^2/Solar Mass
#define Mpcperkm 3.241e-20//Mpc/km used to convert H0 to per seconds

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

static double integrand(double r_z, void*params);

int calc_sigma_r_at_r(double Rp,double Mass,double concentration
		      ,int delta,double*R,double*xi,int NR,
		      double*sigma_r,double*err,cosmology cosmo){
  int i;

  double h = cosmo.h, om = cosmo.om;
  double H0 = h*100.;
  double rhom = om*3.*(H0*Mpcperkm*H0*Mpcperkm)/(8.*PI*G)
    /(h*h*1e12);//SM h^2/pc^2/Mpc

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
  
  do_integral(sigma_r,err,params);
  *sigma_r *= rhom*2;
  *err *= rhom*2;

  gsl_spline_free(spline),gsl_interp_accel_free(acc);
  gsl_integration_workspace_free(workspace);
  free(params);
  return 0;
}


static int do_integral(double*sigma_r,double*err,integrand_params*params){
  gsl_function F;
  F.function=&integrand;
  F.params=params;
  
  double lrmin = params->lrmin;
  double lrmax = params->lrmax;
  gsl_integration_workspace*workspace=params->workspace;

  double result,abserr;
  int status = gsl_integration_qag(&F,lrmin-10,lrmax,TOL,TOL/10.,
				   workspace_size,6,workspace,&result,&abserr);
  *sigma_r = result;
  *err= abserr;
  return status;
}

static double integrand(double rz, void*params){
  rz = exp(rz);
  integrand_params pars=*(integrand_params *)params;
  gsl_spline*spline = pars.spline;//Xi_hm(R) spline
  gsl_interp_accel*acc = pars.acc;
  double rp = pars.rperp;
  cosmology cosmo = pars.cosmo;
  double rmin = exp(pars.lrmin),rmax = exp(pars.lrmax);
  double arg = sqrt(rz*rz+rp*rp);
  double Mass = pars.Mass;
  double concentration = pars.concentration;
  int delta = pars.delta;
  double answer = 0;

  if(arg>rmax)
    return answer;
  else if(arg<rmin){
    calc_xi_nfw(&arg,1,Mass,concentration,delta,&answer,cosmo);
    return rz*answer;
  }else
    return rz*gsl_spline_eval(spline,arg,acc);
}
