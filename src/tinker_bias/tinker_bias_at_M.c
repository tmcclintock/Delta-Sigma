#include "tinker_bias_at_M.h"

typedef struct integrand_params{
  gsl_spline *spline;
  gsl_interp_accel *acc;
  gsl_integration_workspace * workspace;
  double r;
  double lkmin;
  double lkmax;
}integrand_params;

static int do_integral(double*nu,integrand_params*params);

static double integrand(double lk,void*params);

int tinker_bias_at_M(double M,double*k,double*P,int N,double*bias,double*nu,int delta,cosmology cosmo){

  double om = cosmo.om;
  double rhom = om*rhomconst;//SM h^2/Mpc^3
  double R=pow(M/(1.33333333333*PI*rhom),0.3333333333);//Lagrangian radius Mpc/h

  double y = log10(delta);
  double xp = exp(-1.0*pow(4./y,4.));
  double A = 1.+0.24*y*xp, a = 0.44*y-0.88;
  double B = 0.183, b = 1.5;
  double C = 0.019+0.107*y+0.19*xp, c = 2.4;

  //Prepare the splines and the arguments to the integrands
  gsl_spline*spline=gsl_spline_alloc(gsl_interp_cspline,N);
  gsl_spline_init(spline,k,P,N);
  gsl_interp_accel*acc=gsl_interp_accel_alloc();
  gsl_integration_workspace*workspace
    =gsl_integration_workspace_alloc(workspace_size);
  integrand_params *params=malloc(sizeof(integrand_params));
  params->spline=spline;
  params->acc=acc;
  params->workspace=workspace;
  params->r=R;
  params->lkmin=log(k[0]);
  params->lkmax=log(k[N-1]);

  //Calculate nu
  do_integral(nu,params);
  //And the actual bias
  *bias = 1 
      - A*pow(*nu,a) / (pow(*nu,a)+pow(delta_c,a))
      + B*pow(*nu,b) 
      + C*pow(*nu,c);

  gsl_spline_free(spline),gsl_interp_accel_free(acc);
  gsl_integration_workspace_free(workspace);
  free(params);
  return 0;
}

int do_integral(double*nu,integrand_params*params){
  gsl_function F;
  F.function=&integrand;
  F.params=params;

  double lkmin = params->lkmin;
  double lkmax = params->lkmax;
  gsl_integration_workspace*workspace=params->workspace;

  double result,abserr;
  int status = gsl_integration_qag(&F,lkmin,lkmax,BIAS_TOL,BIAS_TOL/10.,workspace_size,6,workspace,&result,&abserr);
  *nu = delta_c/sqrt(result/(2.*PI*PI));
  return status;
}

double integrand(double lk, void*params){
  integrand_params pars = *(integrand_params*)params;
  gsl_spline*spline = pars.spline;
  gsl_interp_accel*acc = pars.acc;
  double R = pars.r;
  
  double k = exp(lk);
  double x  = k*R;
  double P = gsl_spline_eval(spline,k,acc);
  double w = 3.0/x/x/x*(sin(x)-x*cos(x)); //Window function
  return k*k*k*P*w*w;
}
